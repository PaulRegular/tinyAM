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
#' @param method OSA algorithm, as in [add_osa_residuals()]. Only used with
#'   `type = "osa"`. Explicitly supplying it requests a fresh calculation.
#' @param ... Computation controls passed to [RTMB::oneStepPredict()], such as
#'   `trace = FALSE` or `parallel = FALSE`.
#'   `subset` and `conditional` index included observations in likelihood order;
#'   unrequested rows retain missing residuals.
#'   Explicit controls request a fresh calculation; other residual types reject them.
#' @param residual_seed Optional nonnegative integer for `type = "quantile"`.
#'   The default is 123456. Supplying a seed requests a fresh calculation.
#' @details `type = "osa"` reuses cached OSA columns when no computation controls
#' or method are supplied. Otherwise both catch and index residuals are calculated
#' together, using one algorithm at fitted parameter values. It can be expensive,
#' but does not rerun [fit_tam()] or re-estimate fixed parameters. Temporary
#' objectives protect the fitted object's state. Other residual types remain
#' simple extractors. Poisson-link delta lognormal fits use the positive
#' component for log and standardized residuals and the full mean for response
#' residuals. Quantile residuals use the fitted cumulative distribution,
#' drawing uniformly across its probability step at each observed zero. These
#' residuals condition on fitted states. For lognormal observations they equal
#' marginal standardized residuals. Missing and projected observations have
#' missing residuals, and deterministic zeros have no quantile residuals.
#' OSA is available for all-lognormal fits. The supplied fit is not modified.
#' To retain results for later extraction and dashboard use, assign
#' `fit_osa <- add_osa_residuals(fit, method = "oneStepGeneric", trace = FALSE)`.
#' @return A data frame with year, age, type (catch/index), survey, observations,
#'   predictions, residuals, and available unit/fleet columns. OSA output also
#'   identifies `osa_method`. The `residual_type` attribute names the requested
#'   scale, `residual_seed` records the quantile seed, and `osa` retains available
#'   method/settings/coverage metadata.
#' @examples
#' \dontrun{
#' residuals(fit)
#' residuals(fit, type = "log")
#' residuals(fit, type = "quantile")
#' residuals(fit, type = "osa", method = "oneStepGeneric", trace = FALSE)
#' }
#' @seealso [tidy_obs_pred()], [plot_bubbles()]
#' @importFrom stats residuals
#' @export
residuals.tam_fit <- function(object, type = c("standardized", "log", "response", "osa", "quantile"),
                             method = "oneStepGaussianOffMode", ..., residual_seed = NULL) {
  object <- .require_tam_fit(object, arg = "object")
  type <- match.arg(type)
  controls <- list(...)
  explicit <- !missing(method) || length(controls) > 0L
  .osa_method(method, controls)
  if (type != "osa" && explicit) cli::cli_abort("Additional computation controls are only used with type = 'osa'.")
  if (type != "quantile" && !is.null(residual_seed)) {
    cli::cli_abort("{.arg residual_seed} is only used with type = 'quantile'.")
  }
  if (type == "osa") .check_osa_family(object)
  tables <- .observation_tables(object)
  if (type == "quantile") {
    if (!is.null(residual_seed) || any(vapply(tables, function(d) is.null(d$quantile_res), logical(1)))) {
      tables <- .add_quantile_residuals(tables, residual_seed)
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
    pld <- which(d$family == "pld_lognormal")
    residual <- switch(type, standardized = d$std_res, osa = d$osa_res,
      quantile = if (is.null(d$quantile_res)) d$std_res else d$quantile_res,
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
    if (type %in% c("quantile", "response") && length(pld)) {
      keep[pld] <- is.finite(d$obs[pld]) & d$obs[pld] >= 0 & is.finite(residual[pld])
    }
    if ("is_proj" %in% names(d)) keep <- keep & !d$is_proj
    residual[!keep] <- NA_real_
    if (type == "osa" && !"osa_method" %in% names(d)) {
      d$osa_method <- rep(if (is.null(metadata$method)) NA_character_ else metadata$method, nrow(d))
    }
    columns <- c("year", "age", "survey", "fleet", "obs", "pred", "unit", "is_proj",
      "family", "pred_scale")
    if (type == "osa") columns <- c(columns, "osa_method")
    d <- d[intersect(columns, names(d))]
    d$type <- name; d$residual <- residual
    if (!"survey" %in% names(d)) d$survey <- NA_character_
    d
  }))
  attr(out, "residual_type") <- type
  if (type == "osa" && !is.null(metadata)) attr(out, "osa") <- metadata
  if (type == "quantile") attr(out, "residual_seed") <- attr(tables, "residual_seed", exact = TRUE)
  out
}

#' Calculate and retain observation residuals
#'
#' Store observation diagnostics in a returned copy of a fitted model for use
#' with [residuals()] and [vis_tam()]. Quantile residuals condition on the fitted
#' states. One-step-ahead residuals are available for all-lognormal fits.
#'
#' @param fit A fitted `tam_fit`.
#' @param type `"quantile"` (default) or `"osa"`.
#' @param ... For `type = "osa"`, arguments to [add_osa_residuals()].
#' @param residual_seed Optional nonnegative integer for randomizing exact zeros
#'   in quantile residuals. Defaults to the stored seed, or 123456 for a new fit.
#' @return A `tam_fit` with retained `quantile_res` or `osa_res` columns in
#'   `$obs_pred`. The supplied fit and the caller's random-number state are
#'   unchanged by the quantile calculation.
#' @seealso [residuals.tam_fit()], [add_osa_residuals()], [tidy_obs_pred()]
#' @export
add_residuals <- function(fit, type = c("quantile", "osa"), ..., residual_seed = NULL) {
  fit <- .require_tam_fit(fit, arg = "fit")
  type <- match.arg(type)
  if (type == "osa") {
    if (!is.null(residual_seed)) cli::cli_abort("{.arg residual_seed} is only used for quantile residuals.")
    return(add_osa_residuals(fit, ...))
  }
  if (length(list(...))) cli::cli_abort("Additional computation controls are only used with type = 'osa'.")
  fit$obs_pred <- .add_quantile_residuals(.observation_tables(fit), residual_seed)
  fit
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
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  stats::runif(n)
}

.add_quantile_residuals <- function(tables, seed = NULL) {
  if (is.null(seed)) seed <- attr(tables, "residual_seed", exact = TRUE)
  seed <- .observation_residual_seed(seed)
  index <- lapply(tables, function(d) which(d$family == "pld_lognormal" &
    d$is_observed & !d$is_structural & d$obs == 0))
  uniforms <- .residual_uniform(sum(lengths(index)), seed)
  offset <- 0L
  for (name in names(tables)) {
    d <- tables[[name]]
    values <- if (is.null(d$std_res)) rep(NA_real_, nrow(d)) else d$std_res
    positive <- which(d$family == "pld_lognormal" & d$is_observed &
      !d$is_structural & d$obs > 0)
    if (length(positive)) {
      z <- values[positive]
      p_zero <- d$log_p_zero[positive]
      p_positive <- d$log_p_positive[positive]
      upper <- p_positive + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
      use_lower <- which(upper > log(.5))
      values[positive] <- stats::qnorm(upper, lower.tail = FALSE, log.p = TRUE)
      if (length(use_lower)) {
        lower_component <- p_positive[use_lower] + stats::pnorm(z[use_lower], log.p = TRUE)
        zero_component <- p_zero[use_lower]
        largest <- pmax(zero_component, lower_component)
        lower <- largest + log1p(exp(pmin(zero_component, lower_component) - largest))
        lower[largest == -Inf] <- -Inf
        values[positive[use_lower]] <- stats::qnorm(lower, log.p = TRUE)
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
    d$quantile_res <- values
    tables[[name]] <- d
  }
  attr(tables, "residual_seed") <- seed
  tables
}

