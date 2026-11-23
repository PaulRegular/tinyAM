#' Extract observation residuals
#'
#' Return keyed catch and survey rows for plotting or further checks. This is a
#' convenience extractor. Standardized residuals use fitted states and marginal
#' observation SDs. Shrinkage affects these residuals, and standardization
#' leaves any across-age correlation in place.
#' @param object A fitted `tam_fit`.
#' @param type `"standardized"` (default), `"log"` (log observation minus log
#'   median), `"response"` (observation minus prediction), `"quantile"`, or
#'   `"osa"`.
#' @param method For `type = "quantile"`, use `"mle-mvn"` for a Gaussian draw of
#'   the random effects, `"mle-eb"` for the fitted states, or `"auto"` to choose
#'   by model. With `NULL`, use the retained method or the automatic choice.
#'   For `type = "osa"`, select an algorithm from [add_osa_residuals()].
#'   Supplying a method explicitly requests a fresh calculation.
#' @param ... Computation controls passed to [RTMB::oneStepPredict()], such as
#'   `trace = FALSE` or `parallel = FALSE`.
#'   `subset` and `conditional` index included observations in likelihood order;
#'   unrequested rows retain missing residuals.
#'   Explicit controls request a fresh calculation; other residual types reject them.
#' @param residual_seed Optional nonnegative integer for `type = "quantile"`.
#'   Controls the Gaussian state draw and randomization of recorded zeros.
#'   Defaults to the stored seed, or 123456. Supplying a seed requests a fresh calculation.
#' @details `type = "osa"` reuses cached OSA columns when no computation controls
#' or method are supplied. Otherwise both catch and index residuals are calculated
#' together, using one algorithm at fitted parameter values. It can be expensive,
#' but does not rerun [fit_tam()] or re-estimate fixed parameters. Temporary
#' objectives protect the fitted object's state. For Gamma and generalized gamma, standardized residuals
#' use the mean and SD of the log of positive observations. Log residuals use
#' the positive median. Delta lognormal standardized residuals also use the
#' positive component. Tweedie standardized residuals are response-scale
#' Pearson residuals. Delta and Tweedie response residuals use the mean
#' including zeros.
#'
#' For `"mle-mvn"`, fixed effects stay at their estimates and all fitted random
#' effects are sampled together. The Gaussian approximation uses their
#' conditional mode and precision. Catch and survey distributions are evaluated
#' at this draw, using the approach described in the
#' [sdmTMB residual documentation](https://sdmtmb.github.io/sdmTMB/reference/residuals.sdmTMB.html).
#' The fitted-state method `"mle-eb"` evaluates the same distributions at the
#' mode. For lognormal observations at the mode, quantile and standardized
#' residuals agree. Both quantile methods randomize uniformly across the
#' cumulative probability step at each sampled zero. Missing observations,
#' projections and structural zeros have missing quantile residuals.
#' Marginal quantiles retain any across-age observation correlation.
#'
#' The automatic choice uses `"mle-mvn"` when there are fitted random effects and
#' the conditional draw is usable. Otherwise it uses `"mle-eb"`. If a draw
#' fails, the returned metadata records the reason. An explicit request for
#' `"mle-mvn"` in that case gives an error.
#' OSA is available for all-lognormal fits. The supplied fit is not modified.
#' To retain results for later extraction and dashboard use, assign
#' `fit_osa <- add_osa_residuals(fit, method = "oneStepGeneric", trace = FALSE)`.
#' @return A data frame with year, age, type (catch/index), survey, observations,
#'   predictions, residuals, and available unit/fleet columns. OSA output also
#'   identifies `osa_method`. The `residual_type` attribute names the requested
#'   scale. Quantile output identifies the method and any fallback in the table.
#'   The `quantile` attribute records the requested and used methods, status and
#'   seed. The `osa` attribute retains available method, settings and coverage.
#' @examples
#' \dontrun{
#' residuals(fit)
#' residuals(fit, type = "log")
#' residuals(fit, type = "quantile")
#' residuals(fit, type = "quantile", method = "mle-eb")
#' residuals(fit, type = "osa", method = "oneStepGeneric", trace = FALSE)
#' }
#' @seealso [tidy_obs_pred()], [plot_bubbles()]
#' @importFrom stats residuals
#' @export
residuals.tam_fit <- function(object, type = c("standardized", "log", "response", "osa", "quantile"),
                             method = NULL, ..., residual_seed = NULL) {
  object <- .require_tam_fit(object, arg = "object")
  type <- match.arg(type)
  controls <- list(...)
  .validate_residual_controls(type, method, controls)
  if (type != "quantile" && !is.null(residual_seed)) {
    cli::cli_abort("{.arg residual_seed} is only used with type = 'quantile'.")
  }
  if (type == "osa") {
    .check_osa_family(object)
    explicit <- !is.null(method) || length(controls) > 0L
    if (is.null(method)) method <- "oneStepGaussianOffMode"
    .osa_method(method, controls)
  }
  tables <- .observation_tables(object)
  if (type == "quantile") {
    tables <- .observation_quantiles(object, tables, method, residual_seed)
    quantile <- attr(tables, "quantile", exact = TRUE)
    if (quantile$status == "fallback") {
      cli::cli_inform("Using fitted-state quantile residuals: {quantile$detail}")
    }
  }
  metadata <- attr(tables, "osa", exact = TRUE)
  if (is.null(metadata)) metadata <- object$osa
  missing_osa <- vapply(tables, function(d) is.data.frame(d) && nrow(d) > 0L &&
    !"osa_res" %in% names(d) && any(is.finite(d$obs) & d$obs > 0 & if (is.null(d$is_proj)) TRUE else !d$is_proj), logical(1))
  if (type == "osa" && (explicit || any(missing_osa))) {
    computed <- .compute_osa_residuals(object, method, controls)
    metadata <- computed$metadata
    tables <- .attach_osa_tables(tables, computed$tables, metadata)
  }

  out <- .bind_rows(lapply(names(tables), function(name) {
    d <- tables[[name]]
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    zero_capable <- .is_zero_family(d$family)
    residual <- switch(type, standardized = d$std_res, osa = d$osa_res,
      quantile = d$quantile_res,
      response = d$obs - d$pred, log = {
        if (!is.null(d$log_res)) {
          d$log_res
        } else {
          value <- suppressWarnings(log(d$obs) - log(d$pred))
          if (!is.null(d$sd) && !is.null(d$std_res)) {
            underflow <- which(d$pred == 0 & d$obs > 0 & is.finite(d$std_res) & is.finite(d$sd))
            value[underflow] <- d$std_res[underflow] * d$sd[underflow]
          }
          value
        }
      })
    if (is.null(residual)) residual <- rep(NA_real_, nrow(d))
    keep <- is.finite(d$obs) & d$obs > 0 & is.finite(residual)
    if (type %in% c("quantile", "response", "standardized")) {
      zero <- which(d$obs == 0 & if (type == "standardized") d$family %in% "tweedie" else zero_capable)
      keep[zero] <- is.finite(residual[zero])
    }
    if ("is_proj" %in% names(d)) keep <- keep & !d$is_proj
    residual[!keep] <- NA_real_
    if (type == "osa" && !"osa_method" %in% names(d)) {
      d$osa_method <- rep(if (is.null(metadata$method)) NA_character_ else metadata$method, nrow(d))
    }
    columns <- c("year", "age", "survey", "fleet", "obs", "pred", "unit", "is_proj",
      "family", "delta_type", "pred_scale")
    if (type == "osa") columns <- c(columns, "osa_method")
    if (type == "quantile") columns <- c(columns, "quantile_method", "quantile_status", "quantile_detail")
    d <- d[intersect(columns, names(d))]
    d$type <- name; d$residual <- residual
    if (!"survey" %in% names(d)) d$survey <- NA_character_
    d
  }))
  attr(out, "residual_type") <- type
  if (type == "osa" && !is.null(metadata)) attr(out, "osa") <- metadata
  if (type == "quantile") {
    attr(out, "residual_seed") <- quantile$seed
    attr(out, "quantile") <- quantile
  }
  out
}

#' Calculate and retain observation residuals
#'
#' Store observation diagnostics in a returned copy of a fitted model for use
#' with [residuals()] and [vis_tam()]. One-step-ahead residuals are available
#' for all-lognormal fits. Standardized residuals are already in the fit.
#'
#' @param fit A fitted `tam_fit`.
#' @param type `"quantile"` (default) or `"osa"`.
#' @param ... For `type = "osa"`, arguments to [add_osa_residuals()].
#' @param method Residual method, as in [residuals.tam_fit()].
#' @param residual_seed Optional nonnegative integer for the Gaussian draw and
#'   randomization of observed zeros. Defaults to the stored seed or 123456.
#' @return A `tam_fit` with retained `quantile_res` or `osa_res` columns in
#'   `$obs_pred`. The supplied fit and the caller's random-number state are
#'   unchanged by the quantile calculation.
#' @seealso [residuals.tam_fit()], [add_osa_residuals()], [tidy_obs_pred()]
#' @export
add_residuals <- function(fit, type = c("quantile", "osa"), ..., residual_seed = NULL, method = NULL) {
  fit <- .require_tam_fit(fit, arg = "fit")
  type <- match.arg(type)
  .validate_residual_controls(type, method, list(...))
  if (type == "osa") {
    if (!is.null(residual_seed)) cli::cli_abort("{.arg residual_seed} is only used for quantile residuals.")
    if (is.null(method)) return(add_osa_residuals(fit, ...))
    return(add_osa_residuals(fit, method = method, ...))
  }
  fit$obs_pred <- .observation_quantiles(fit, .observation_tables(fit), method, residual_seed)
  fit
}

.validate_residual_controls <- function(type, method, controls) {
  if ("methods" %in% names(controls)) cli::cli_abort("Use {.arg method}, not {.arg methods}.")
  if (type != "osa" && length(controls)) {
    cli::cli_abort("Additional computation controls are only used with type = 'osa'.")
  }
  if (!is.null(method) && !type %in% c("osa", "quantile")) {
    cli::cli_abort("{.arg method} is only used with type = 'osa' or type = 'quantile'.")
  }
}

.quantile_method <- function(method) {
  choices <- c("auto", "mle-mvn", "mle-eb")
  if (!is.character(method) || length(method) != 1L || !method %in% choices) {
    cli::cli_abort("Quantile method must be one of {.val {choices}}.")
  }
  method
}

.fitted_random_indices <- function(fit) {
  random <- fit$obj$env$lrandom
  if (is.function(random)) which(random()) else integer()
}

.sampled_observation_quantiles <- function(fit, seed) {
  random <- .fitted_random_indices(fit)
  if (!length(random)) cli::cli_abort("There are no fitted random effects to sample.")
  precision <- .diagnostic_joint_precision(fit)
  mode <- fit$parameter_values
  map <- fit$dat$obs_map
  .with_diagnostic_seed(seed, {
    zero <- which(map$is_observed & .is_zero_family(map$family) & fit$dat$obs_values == 0)
    uniforms <- stats::runif(length(zero))
    # At fixed MLEs the random/random block of the joint precision is conditional.
    mode[random] <- mode[random] + .precision_perturbation(precision[random, random, drop = FALSE])
    if (any(!is.finite(mode[random]))) cli::cli_abort("The sampled random effects are non-finite.")
    sampled <- list(dat = fit$dat, rep = fit$obj$report(mode))
    tables <- .add_quantile_residuals(.tidy_obs_predictions(sampled, positive_medians = FALSE),
      seed, uniforms = uniforms)
    for (name in names(tables)) {
      included <- map$is_observed[map$type == name]
      if (any(!is.finite(tables[[name]]$quantile_res[included]))) {
        cli::cli_abort("The sampled predictions did not yield finite residuals for all likelihood observations.")
      }
      tables[[name]]$quantile_res[!included] <- NA_real_
    }
    tables
  })
}

.copy_observation_quantiles <- function(tables, source, columns) {
  for (name in names(tables)) {
    rows <- .observation_match(tables[[name]], source[[name]])
    tables[[name]][columns] <- source[[name]][rows, columns, drop = FALSE]
  }
  tables
}

.observation_quantiles <- function(fit, tables, method = NULL, seed = NULL) {
  cached <- fit$obs_pred
  retained <- attr(cached, "quantile", exact = TRUE)
  if (is.null(method) && is.null(seed) && !is.null(retained)) {
    tables <- .copy_observation_quantiles(tables, cached,
      c("quantile_res", "quantile_method", "quantile_status", "quantile_detail"))
    attr(tables, "quantile") <- retained
    attr(tables, "residual_seed") <- retained$seed
    return(tables)
  }
  if (is.null(seed)) seed <- attr(cached, "residual_seed", exact = TRUE)
  seed <- .observation_residual_seed(seed)
  if (is.null(method)) method <- if (is.null(retained)) "auto" else retained$requested
  requested <- .quantile_method(method)
  method <- requested
  if (requested == "auto") {
    method <- if (length(.fitted_random_indices(fit))) "mle-mvn" else "mle-eb"
  }
  detail <- NA_character_
  if (method == "mle-mvn") {
    sampled <- tryCatch(.sampled_observation_quantiles(fit, seed), error = identity)
    if (inherits(sampled, "error")) {
      detail <- conditionMessage(sampled)
      if (requested == "mle-mvn") cli::cli_abort(c("Sampled-state quantile residuals are unavailable.",
        "x" = "{detail}", "i" = "Use method = 'mle-eb' to evaluate quantiles at the fitted states."))
      method <- "mle-eb"
    } else {
      tables <- .copy_observation_quantiles(tables, sampled, "quantile_res")
    }
  }
  if (method == "mle-eb") tables <- .add_quantile_residuals(tables, seed)
  status <- if (is.na(detail)) "ok" else "fallback"
  for (name in names(tables)) {
    n <- nrow(tables[[name]])
    tables[[name]]$quantile_method <- rep(method, n)
    tables[[name]]$quantile_status <- rep(status, n)
    tables[[name]]$quantile_detail <- rep(detail, n)
  }
  attr(tables, "quantile") <- list(method = method, requested = requested, status = status, detail = detail, seed = seed)
  attr(tables, "residual_seed") <- seed
  tables
}

.observation_residual_seed <- function(seed = NULL) {
  if (is.null(seed)) return(123456L)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    cli::cli_abort("{.arg residual_seed} must be a nonnegative integer.")
  }
  as.integer(seed)
}

.residual_uniform <- function(n, seed) {
  if (!n) return(numeric())
  .with_diagnostic_seed(seed, stats::runif(n))
}

.add_quantile_residuals <- function(tables, seed = NULL, uniforms = NULL) {
  if (is.null(seed)) seed <- attr(tables, "residual_seed", exact = TRUE)
  seed <- .observation_residual_seed(seed)
  index <- lapply(tables, function(d) {
    if (is.null(d$family)) return(integer())
    which(.is_zero_family(d$family) & d$is_observed & !d$is_structural & d$obs == 0)
  })
  if (is.null(uniforms)) uniforms <- .residual_uniform(sum(lengths(index)), seed)
  offset <- 0L
  for (name in names(tables)) {
    d <- tables[[name]]
    values <- if (is.null(d$std_res)) rep(NA_real_, nrow(d)) else d$std_res
    positive <- if (is.null(d$family)) integer() else
      which(d$family != "lognormal" & d$is_observed & !d$is_structural & d$obs > 0)
    for (family in unique(d$family[positive])) {
      rows <- positive[d$family[positive] == family]
      if (family == "tweedie") {
        log_mean <- d$log_obs_positive_mean[rows] + d$log_p_positive[rows]
        log_obs <- log(d$obs[rows])
        log_upper <- .tweedie_log_cdf(log_obs, log_mean, d$sd[rows], d$tweedie_power[rows],
          lower.tail = FALSE)
        values[rows] <- stats::qnorm(log_upper, lower.tail = FALSE, log.p = TRUE)
        use_lower <- which(log_upper > log(.5) | is.na(log_upper))
        if (length(use_lower)) {
          i <- rows[use_lower]
          log_lower <- .tweedie_log_cdf(log_obs[use_lower], log_mean[use_lower],
            d$sd[i], d$tweedie_power[i], lower.tail = TRUE)
          values[i] <- stats::qnorm(log_lower, log.p = TRUE)
        }
        next
      }
      tail <- if (family == "delta_lognormal") {
        function(lower) stats::pnorm(values[rows], lower.tail = lower, log.p = TRUE)
      } else {
        function(lower) .positive_log_cdf(family, log(d$obs[rows]),
          d$log_obs_positive_mean[rows], d$sd[rows], d$gengamma_Q[rows], lower.tail = lower)
      }
      log_lower <- tail(TRUE)
      log_upper <- d$log_p_positive[rows] + tail(FALSE)
      values[rows] <- stats::qnorm(log_upper, lower.tail = FALSE, log.p = TRUE)
      use_lower <- which(log_upper > log(.5))
      if (length(use_lower)) {
        i <- rows[use_lower]
        lower_component <- d$log_p_positive[i] + log_lower[use_lower]
        zero_component <- d$log_p_zero[i]
        largest <- pmax(zero_component, lower_component)
        lower <- largest + log1p(exp(pmin(zero_component, lower_component) - largest))
        lower[which(largest == -Inf)] <- -Inf
        values[i] <- stats::qnorm(lower, log.p = TRUE)
      }
    }
    zero <- index[[name]]
    if (length(zero)) {
      u <- uniforms[offset + seq_along(zero)]
      offset <- offset + length(zero)
      # A zero probability can underflow on the ordinary probability scale.
      values[zero] <- stats::qnorm(d$log_p_zero[zero] + log(u), log.p = TRUE)
    }
    values[!is.finite(values)] <- NA_real_
    if ("is_proj" %in% names(d)) values[which(d$is_proj)] <- NA_real_
    d$quantile_res <- values
    tables[[name]] <- d
  }
  attr(tables, "residual_seed") <- seed
  tables
}

