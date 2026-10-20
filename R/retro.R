
# Resolve targets using raw inputs, before any modeled plus-group pooling.
.hindcast_controls <- function(fit, years, catch, F_mult) {
  supplied <- if (is.null(F_mult)) catch else F_mult
  kind <- if (is.null(catch)) "F_mult" else "catch"
  if (!is.null(catch) && !is.null(F_mult)) cli::cli_abort("Specify catch or F_mult, not both, for hindcasts.")
  observed <- identical(catch, "observed")
  if (is.null(supplied)) supplied <- 1
  if (!observed && (!is.numeric(supplied) || !length(supplied) || any(!is.finite(supplied)) || any(supplied < 0))) {
    cli::cli_abort("Hindcast catch must be 'observed', NULL, or non-negative annual values; F_mult must be non-negative numeric values.")
  }
  if (!observed && !is.null(names(supplied))) {
    named_years <- suppressWarnings(as.numeric(names(supplied)))
    if (anyDuplicated(names(supplied)) || any(!is.finite(named_years)) || any(named_years != trunc(named_years))) {
      cli::cli_abort("Name hindcast controls with unique integer forecast years.")
    }
  } else if (!observed && length(supplied) > 1L) {
    cli::cli_abort("Name a multi-year hindcast catch or F_mult vector by forecast year.")
  }
  out <- data.frame(year = years, control = kind, target = NA_real_, target_source = NA_character_, error = NA_character_)
  source <- fit$refit_args$data
  ages <- sort(unique(c(fit$dat$ages, fit$dat$plus_ages)))
  for (i in seq_along(years)) {
    value <- if (observed) NULL else if (is.null(names(supplied))) supplied else supplied[as.character(years[i])]
    if (length(value) && !is.na(value)) {
      out$target[i] <- value
      out$target_source[i] <- if (kind == "catch") "supplied catch" else "F multiplier"
    } else if (kind == "catch") {
      d <- source$catch[source$catch$year == years[i], ]
      w <- source$weight[source$weight$year == years[i], ]
      counts <- d$obs[match(ages, d$age)]
      weights <- w$obs[match(ages, w$age)]
      if (any(!is.finite(counts)) || any(!is.finite(weights))) {
        out$error[i] <- "Observed catch biomass unavailable: incomplete catch-at-age or weights; supply an annual catch override."
      } else {
        out$target[i] <- sum(counts * weights)
        out$target_source[i] <- "observed catch-at-age x supplied weight"
      }
    } else {
      out$error[i] <- "No F multiplier supplied for this forecast year."
    }
  }
  out
}

# Pair forecasts with the original held-out indices, including when a reference
# fit failed numerical checks. Future observations never enter the fold fit.
.hindcast_index_data <- function(predictions, observations) {
  keys <- c("year", "age", "survey")
  p <- predictions[predictions$is_proj, c(keys, "pred"), drop = FALSE]
  o <- observations[!observations$is_proj, c(keys, "obs"), drop = FALSE]
  p$obs <- rep(NA_real_, nrow(p)); p$is_proj <- rep(TRUE, nrow(p))
  o$pred <- rep(NA_real_, nrow(o)); o$is_proj <- rep(FALSE, nrow(o))
  rbind(p, o[, names(p), drop = FALSE])
}


#' Run a retrospective (peel) analysis
#'
#' @title Retrospective fits and catch-conditioned hindcasts
#'
#' @description
#' Creates a sequence of terminal-year peels and refits the model on each
#' truncated year range.
#'
#' @details
#' Ordinary retrospectives fit historical years only, removing inherited projections.
#' `folds` peels plus the full historical reference are returned. Existing
#' projection years never become terminal assessment years.
#'
#' [fit_hindcast()] additionally projects each peeled fit one year ahead and
#' scores survey indices on the log scale. The default conditions on observed
#' aggregate catch biomass in that held-out year. This is a conditional assessment
#' of survey forecast skill, not a prediction of future removals or a forecast
#' made with only information available at the original assessment date.
#'
#' Catch biomass is derived from complete raw catch-at-age times supplied weights,
#' including retained ages within the model plus group. Stock weights may differ
#' from fishery weights: supply authoritative annual catch totals when appropriate.
#' Missing catch cells are not zero. Only the aggregate target is used; held-out
#' catch composition, survey observations and biological inputs do not enter the
#' fit. Forecast biology is based on training data, with `n_mean = 1`.
#'
#' The full historical reference has no forecast, so unknown catch beyond the
#' observed period is not required. Failed fits and unavailable targets are
#' excluded from stacked predictions but retained in `fold_diagnostics`.
#' Mohn's rho is unavailable if the full historical reference failed.
#' Scores describe retained observations; compare models only over identical
#' forecast observations. Each survey observation receives equal weight.
#'
#' Refits use [furrr::future_map()]; set your own [future::plan()] for parallel
#' execution. The supplied fit is unchanged.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param folds Integer; number of terminal peels (default `2`).
#' @param start_from_fit Logical; use terminal parameter estimates as starting
#'    values for retrospective fits?
#' @param hindcast Logical; add one-year projections to peeled fits? Prefer
#'   [fit_hindcast()] for this workflow.
#' @param catch Hindcast catch control: `"observed"` (default), a nonnegative
#'   scalar applied to all forecast years, a numeric vector named by forecast
#'   year, or `NULL` for status-quo F. Named values override observed catches;
#'   remaining years use derived observed catches. Values use native biomass
#'   units (catch numbers times individual weight), without automatic conversion.
#' @param F_mult Optional hindcast multiplier, scalar or named by forecast year.
#'   An explicitly supplied multiplier replaces the implicit observed-catch
#'   default. Explicit non-NULL `catch` and `F_mult` together are an error.
#'   These controls apply only when `hindcast = TRUE`.
#' @param grad_tol Numeric tolerance for `max|grad|`. If `NULL` (default),
#'   the tolerance stored on `fit` (from [fit_tam()]) is used; otherwise the
#'   supplied value is passed to [check_tam()]. Folds failing numerical checks
#'   are dropped; structural and advisory findings remain separate.
#' @param progress Logical; show progress bar using [progressr::with_progress()].
#' @param globals Character vector naming global objects to supply to the workers.
#'
#' @return
#' A list whose elements are:
#'   - `obs_pred` — a named list of stacked data frames (e.g., `catch`, `index`);
#'   - `pop` — a named list of stacked data frames (e.g., `ssb`, `N`, `M`,
#'   `mu_M`, `F`, `mu_F`, `Z`, `ssb_mat`);
#'   - `mohns_rho` - data frame of Mohn's rho (measure of retrospective bias) values
#'     for each data frame within `pop` (calculated using [compute_mohns_rho()]
#'     function);
#'   - `hindcast_rmse` (when `hindcast = TRUE`) - root mean squared error for
#'     use as an overall measure of the forecast skill of the model (calculated
#'     from survey index observations only using the [compute_hindcast_rmse()]
#'     function); and,
#'   - `fits` - retained fitted objects, named by terminal year; and
#'   - `fold_diagnostics` - controls, targets, convergence, warnings, failure
#'     reasons and scored-index counts for every requested fold.
#' The `obs_pred` and `pop` objects are created using [tidy_tam()]. A `fold` column
#' is added to each data.frame within these objects which specifies the terminal
#' year.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' \donttest{
#' if (interactive()) {
#'   # Choose your parallel plan (set once per session)
#'   future::plan(future::multisession, workers = 4)
#'   retros <- fit_retro(fit, folds = 5, progress = TRUE)
#'   head(retros$pop$ssb)
#'   head(retros$mohns_rho)
#'
#'   vis_tam(retros$fits)
#'   hc <- fit_hindcast(fit, folds = 5, progress = TRUE)
#'   hc$hindcast_rmse  # survey index forecasts conditional on observed removals
#'   hc$fold_diagnostics
#'   # Compare with a forecast of status-quo fishing:
#'   hc_F <- fit_hindcast(fit, folds = 5, F_mult = 1, progress = TRUE)
#'   # A named annual catch vector overrides derived catches for those years.
#'   # fit_hindcast(fit, folds = 5, catch = c(`2020` = 1000))
#' }
#' }
#'
#' @seealso
#' [fit_tam()], [sim_tam()]
#'
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr with_progress progressor
#' @importFrom stats update
#'
#' @export
fit_retro <- function(
    fit,
    folds = 2,
    start_from_fit = FALSE,
    hindcast = FALSE,
    grad_tol = NULL,
    progress = TRUE,
    globals = NULL,
    catch = "observed",
    F_mult = NULL
) {
  .require_tam_fit(fit, arg = "fit")
  if (missing(catch) && !is.null(F_mult)) catch <- NULL
  if (hindcast && !is.null(catch) && !is.null(F_mult)) cli::cli_abort("Specify catch or F_mult, not both, for hindcasts.")
  if (!hindcast && (!missing(catch) || !is.null(F_mult))) {
    cli::cli_abort("catch and F_mult controls apply only to hindcasts; use proj_settings for ordinary model projections.")
  }
  if (is.null(grad_tol)) grad_tol <- if (!is.null(fit$grad_tol)) fit$grad_tol else 1e-2
  min_year <- min(fit$dat$years[!fit$dat$is_proj])
  max_year <- .terminal_year(fit)
  if (!is.numeric(folds) || length(folds) != 1L || !is.finite(folds) ||
      folds != trunc(folds) || folds < 0 || folds > max_year - min_year - 1) {
    cli::cli_abort("folds must be a non-negative integer leaving at least two historical years in every fit.")
  }
  retro_years <- seq.int(max_year - folds, max_year)
  controls <- if (hindcast && folds > 0) .hindcast_controls(fit, retro_years[-length(retro_years)] + 1, catch, F_mult) else NULL
  start_par <- if (start_from_fit) as.list(fit$sdrep, "Estimate") else NULL
  actual <- fit$obs_pred$index
  results <- progressr::with_progress({
    tick <- progressr::progressor(steps = length(retro_years))
    furrr::future_map(seq_along(retro_years), function(i) {
      reference <- retro_years[i] == max_year
      control <- if (hindcast && !reference) controls[i, ] else NULL
      diagnostic <- data.frame(fold = retro_years[i],
                               forecast_year = if (hindcast && !reference) retro_years[i] + 1 else NA_real_,
                               is_reference = reference, control = if (is.null(control)) "none" else control$control,
                               target = if (is.null(control)) NA_real_ else control$target,
                               target_source = if (is.null(control)) NA_character_ else control$target_source,
                               target_unit = if (is.null(control)) NA_character_ else if (control$control == "catch") .population_unit(fit$dat, "total_yield") else "multiplier",
                               is_converged = FALSE, n_scored_index = 0L, error = NA_character_, warnings = "")
      warnings <- character()
      r <- withCallingHandlers(tryCatch({
        if (!is.null(control) && !is.na(control$error)) cli::cli_abort(control$error)
        if (reference && !any(fit$dat$is_proj)) {
          r <- fit
        } else {
          projection <- if (is.null(control)) NULL else {
            x <- list(n_proj = 1, n_mean = 1)
            x[[control$control]] <- control$target
            x
          }
          r <- update(fit, years = min_year:retro_years[i], ages = fit$dat$ages,
                      proj_settings = projection, start_par = start_par, silent = TRUE)
        }
        r$grad_tol <- grad_tol
        r$diagnostics <- check_tam(r, grad_tol = grad_tol)
        r$is_converged <- r$diagnostics$is_converged
        diagnostic$is_converged <- isTRUE(r$is_converged)
        if (!diagnostic$is_converged) diagnostic$error <- "Numerical convergence checks did not pass."
        if (hindcast && !reference && diagnostic$is_converged) {
          d <- r$obs_pred$index
          pairs <- merge(d[d$is_proj, c("year", "age", "survey", "pred")],
                         actual[!actual$is_proj, c("year", "age", "survey", "obs")], by = c("year", "age", "survey"))
          diagnostic$n_scored_index <- sum(is.finite(pairs$pred) & pairs$pred > 0 & is.finite(pairs$obs) & pairs$obs > 0)
        }
        r
      }, error = function(e) {
        diagnostic$error <<- conditionMessage(e)
        NULL
      }), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      diagnostic$warnings <- paste(unique(warnings), collapse = "\n")
      tick()
      list(fit = r, diagnostic = diagnostic)
    }, .options = furrr::furrr_options(seed = 1, packages = "tinyAM", globals = globals))
  }, enable = progress)
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostic"))
  if (any(!diagnostics$is_converged)) {
    cli::cli_warn("Model may not have converged or a projection was unavailable for folds: {paste(diagnostics$fold[!diagnostics$is_converged], collapse = ', ')}. See fold_diagnostics.")
  }
  fits <- lapply(results[diagnostics$is_converged], `[[`, "fit")
  names(fits) <- retro_years[diagnostics$is_converged]
  if (!length(fits)) cli::cli_abort("All folds failed convergence checks or projection validation.")
  out <- c(tidy_tam(model_list = fits, label = "fold"), list(fits = fits, fold_diagnostics = diagnostics))
  rhos <- lapply(names(out$pop), function(nm) {
    d <- out$pop[[nm]]
    if ("age" %in% names(d)) {
      values <- sapply(split(d, d$age), compute_mohns_rho)
      data.frame(metric = nm, age = as.numeric(names(values)), rho = values)
    } else {
      data.frame(metric = nm, age = NA_real_, rho = compute_mohns_rho(d))
    }
  })
  out$mohns_rho <- do.call(rbind, rhos)
  if (!as.character(max_year) %in% names(fits)) out$mohns_rho$rho <- NA_real_
  if (hindcast) out$hindcast_rmse <- compute_hindcast_rmse(.hindcast_index_data(out$obs_pred$index, actual))
  out
}

#' @rdname fit_retro
#' @param ... Further arguments to [fit_retro()], including `folds`,
#'   `start_from_fit`, `grad_tol`, `progress` and `globals`.
#' @export
fit_hindcast <- function(fit, ..., catch = "observed", F_mult = NULL) {
  .require_tam_fit(fit, arg = "fit")
  if (missing(catch) && !is.null(F_mult)) catch <- NULL
  fit_retro(fit, hindcast = TRUE, catch = catch, F_mult = F_mult, ...)
}


