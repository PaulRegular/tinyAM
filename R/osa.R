#' Calculate and retain one-step-ahead residuals
#'
#' Calculate observation diagnostics at fitted parameter values and return a fit
#' containing the results. Use the returned object in [vis_tam()] or [residuals()].
#' No model is refitted, and the supplied fit is not modified.
#'
#' @param fit A fitted `tam_fit`.
#' @param method One algorithm supported by [RTMB::oneStepPredict()]:
#'   `"oneStepGaussianOffMode"` (default), `"oneStepGaussian"`,
#'   `"oneStepGeneric"`, `"fullGaussian"`, or `"cdf"`.
#'   Use one of the first four algorithms for models with correlated
#'   observation errors. RTMB's `"cdf"` method does not support its multivariate
#'   normal density.
#' @param ... Named computation controls for [RTMB::oneStepPredict()], such as
#'   `trace`, `parallel`, `subset`, `conditional`, or `splineApprox`.
#'   `subset` and `conditional` index included log observations in likelihood
#'   order. Their unrequested rows retain missing residuals.
#' @details The algorithms use different predictive approximations and can
#' differ numerically. `oneStepGeneric` uses numerical integration and can be
#' expensive. The existing off-mode Gaussian calculation remains the default.
#' Log observations have a continuous Gaussian likelihood, so `discrete = FALSE`
#' is supplied unless explicitly overridden. Method-dependent upstream defaults,
#' including `reverse`, are retained. When SAM is also loaded, parallel OSA
#' calculations use the existing sequential safeguard.
#'
#' Temporary RTMB objectives protect the fitted parameters, uncertainty and
#' objective state. Only OSA columns and metadata are added to a copy of the
#' existing observation tables; row order, confidence intervals, units and custom
#' columns are preserved. Metadata in `fit$osa` records the method, requested and
#' effective settings, and observation coverage. Failed requested calculations
#' remain missing and are recorded separately from unrequested observations.
#' A calculation error returns no partially updated fit.
#' @return A `tam_fit` with `osa_res` and `osa_method` columns in `$obs_pred`,
#'   and calculation metadata in `$osa` and the observation-table attributes.
#' @examples
#' \dontrun{
#' fit_osa <- add_osa_residuals(fit, method = "oneStepGeneric", trace = FALSE)
#' residuals(fit_osa, type = "osa")
#' vis_tam(fit_osa)
#' }
#' @seealso [residuals.tam_fit()], [tidy_obs_pred()]
#' @export
add_osa_residuals <- function(fit, method = "oneStepGaussianOffMode", ...) {
  fit <- .require_tam_fit(fit, arg = "fit")
  result <- .compute_osa_residuals(fit, method, list(...))
  tables <- .observation_tables(fit)
  tables <- .attach_osa_tables(tables, result$tables, result$metadata)
  out <- fit
  out$obs_pred <- tables
  out$osa <- result$metadata
  out
}

.osa_method <- function(method, controls = list()) {
  if ("methods" %in% names(controls)) {
    cli::cli_abort("Use {.arg method}, not {.arg methods}, to select one OSA algorithm.")
  }
  supported <- c("oneStepGaussianOffMode", "oneStepGaussian", "oneStepGeneric",
                 "fullGaussian", "cdf")
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% supported) {
    cli::cli_abort("{.arg method} must be one of: {paste(supported, collapse = ', ')}.")
  }
  if (length(controls) && (is.null(names(controls)) || anyNA(names(controls)) ||
      any(!nzchar(names(controls))) || anyDuplicated(names(controls)))) {
    cli::cli_abort("OSA computation controls must have unique, nonempty names.")
  }
  reserved <- intersect(names(controls), c("obj", "method", "observation.name", "data.term.indicator"))
  if (length(reserved)) {
    cli::cli_abort("These OSA arguments are supplied by tinyAM: {paste(reserved, collapse = ', ')}.")
  }
  method
}

.osa_indices <- function(x, n, name) {
  if (is.null(x)) return(integer())
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x < 1L | x > n) || any(x != floor(x)) || anyDuplicated(x)) {
    cli::cli_abort("OSA {.arg {name}} must contain distinct observation indices between 1 and {n}.")
  }
  as.integer(x)
}

.compute_osa_residuals <- function(fit, method, controls = list()) {
  method <- .osa_method(method, controls)
  if (is.null(fit$obj) || is.null(fit$parameter_values)) {
    cli::cli_abort("OSA calculation requires the fitted objective and stored parameter values.")
  }
  if (is.null(fit$dat$obs_map) || is.null(fit$dat$observed)) {
    cli::cli_abort("OSA calculation requires the stored observations and likelihood mapping.")
  }
  if (method == "cdf" && .obs_has_correlation(fit$dat)) {
    cli::cli_abort(c(
      "OSA method 'cdf' does not support correlated observation errors.",
      "i" = "Use method = 'oneStepGaussianOffMode', 'oneStepGaussian', 'oneStepGeneric' or 'fullGaussian'."
    ))
  }
  requested_settings <- controls
  n <- length(fit$dat$observed)
  conditional <- .osa_indices(controls$conditional, n, "conditional")
  selected <- if (is.null(controls$subset)) setdiff(seq_len(n), conditional) else
    .osa_indices(controls$subset, n, "subset")
  if (length(intersect(selected, conditional))) {
    cli::cli_abort("OSA {.arg subset} and {.arg conditional} must not overlap.")
  }
  if (is.null(controls$discrete)) controls$discrete <- FALSE
  if (isTRUE(controls$parallel) && "stockassessment" %in% loadedNamespaces()) {
    controls$parallel <- FALSE
    cli::cli_inform("Calculating one-step residuals sequentially while SAM and RTMB are loaded.")
  }
  values <- rep(NA_real_, n)
  if (length(selected)) {
    result <- tryCatch({
      obj <- .fitted_objective(fit)
      do.call(RTMB::oneStepPredict, c(list(obj = obj, method = method), controls))
    }, error = function(e) cli::cli_abort(
      "Could not calculate OSA residuals with {method}: {conditionMessage(e)}", parent = e))
    if (!is.numeric(result$residual) || length(result$residual) != length(selected)) {
      cli::cli_abort("OSA results do not match the requested observation subset.")
    }
    values[selected] <- result$residual
  }
  failed <- selected[!is.finite(values[selected])]
  if (length(failed)) {
    values[failed] <- NA_real_
    cli::cli_warn("OSA method {method} returned unavailable residuals for {length(failed)} requested observations.")
  }
  map <- fit$dat$obs_map
  observed <- map$is_observed
  coverage <- map[observed, intersect(c("type", "year", "age", "survey", "fleet"), names(map)), drop = FALSE]
  coverage$observation <- seq_len(n)
  coverage$requested <- seq_len(n) %in% selected
  coverage$conditional <- seq_len(n) %in% conditional
  coverage$available <- is.finite(values)
  coverage$status <- ifelse(coverage$conditional, "conditional",
    ifelse(!coverage$requested, "not_requested", ifelse(coverage$available, "available", "failed")))
  rownames(coverage) <- NULL
  metadata <- list(method = method, requested_settings = requested_settings,
    settings = c(list(method = method), controls), coverage = coverage,
    n_observations = n, n_requested = length(selected), n_available = sum(is.finite(values)),
    n_failed = length(failed), n_excluded = sum(!observed))
  all_values <- rep(NA_real_, nrow(map))
  all_values[observed] <- values
  tables <- fit$dat$obs[c("catch", "index")]
  for (name in names(tables)) {
    tables[[name]]$osa_res <- all_values[map$type == name]
    tables[[name]]$osa_method <- rep(method, nrow(tables[[name]]))
    attr(tables[[name]], "osa") <- metadata
  }
  attr(tables, "osa") <- metadata
  list(tables = tables, metadata = metadata)
}

.observation_tables <- function(fit) {
  tables <- fit$obs_pred
  absent <- vapply(c("catch", "index"), function(name) !is.data.frame(tables[[name]]), logical(1))
  if (any(absent)) {
    interval <- attr(tables, "interval", exact = TRUE)
    if (is.null(interval)) interval <- attr(fit$pop, "interval", exact = TRUE)
    if (is.null(interval)) interval <- .95
    fallback <- tidy_obs_pred(fit, interval = interval)
    for (name in names(absent)[absent]) tables[[name]] <- fallback[[name]]
    if (is.null(attr(tables, "interval", exact = TRUE))) attr(tables, "interval") <- interval
  }
  tables
}

.observation_match <- function(target, source) {
  keys <- intersect(c("year", "age", "survey", "fleet"), names(source))
  if (!all(keys %in% names(target)) || !length(keys)) {
    cli::cli_abort("OSA observation rows must retain the fitted year, age, survey and fleet keys.")
  }
  key <- function(d) do.call(paste, c(d[keys], sep = "\r"))
  source_key <- key(source)
  rows <- match(key(target), source_key)
  if (anyDuplicated(source_key) || anyNA(rows)) {
    cli::cli_abort("OSA observation rows do not match the fitted inputs uniquely.")
  }
  rows
}

.attach_osa_tables <- function(tables, source, metadata = attr(source, "osa", exact = TRUE)) {
  for (name in intersect(c("catch", "index"), names(tables))) {
    if (!is.data.frame(tables[[name]]) || !is.data.frame(source[[name]]) ||
        !"osa_res" %in% names(source[[name]])) next
    rows <- .observation_match(tables[[name]], source[[name]])
    tables[[name]]$osa_res <- source[[name]]$osa_res[rows]
    tables[[name]]$osa_method <- if ("osa_method" %in% names(source[[name]]))
      as.character(source[[name]]$osa_method[rows]) else rep(
        if (is.null(metadata$method)) NA_character_ else metadata$method, nrow(tables[[name]]))
    if (!is.null(metadata)) attr(tables[[name]], "osa") <- metadata
  }
  if (!is.null(metadata)) attr(tables, "osa") <- metadata
  tables
}

.preserve_observation_extras <- function(tables, cached, metadata = attr(cached, "osa", exact = TRUE)) {
  for (name in intersect(names(tables), names(cached))) {
    extras <- setdiff(names(cached[[name]]), names(tables[[name]]))
    if (length(extras)) {
      rows <- .observation_match(tables[[name]], cached[[name]])
      tables[[name]][extras] <- cached[[name]][rows, extras, drop = FALSE]
    }
  }
  .attach_osa_tables(tables, cached, metadata)
}
