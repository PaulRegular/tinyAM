
#' Merge user-supplied starting values into parameter list
#'
#' @description
#' Internal helper that overlays user-provided starting values onto the default
#' parameter list from `make_par(dat)`. Only matching names are updated.
#'
#' @details
#' Scalars are replaced directly, vectors align by names (if present) or by
#' position, and matrices align by row/column names when available, otherwise
#' by shape. Non-matching elements remain unchanged. This allows reusing
#' estimated parameters from previous fits as starting values for new runs or
#' retros, ensuring consistent initialization without manual trimming.
#'
#' @param par0 Base parameter list (usually from `make_par(dat)`).
#' @param start List of user-supplied starting values with matching structure.
#'
#' @return A parameter list like `par0`, updated with matching values from
#'   `start`.
#'
#' @keywords internal
#' @noRd

.merge_start_par <- function(par0, start) {
  for (nm in intersect(names(par0), names(start))) {
    x0 <- par0[[nm]]
    xs <- start[[nm]]

    # Scalars
    if (is.numeric(x0) && length(dim(x0)) == 0 && length(x0) == 1L) {
      if (is.numeric(xs) && length(xs) == 1L) par0[[nm]] <- xs
      next
    }

    # Vectors with names
    if (is.null(dim(x0))) {
      if (!is.null(names(x0)) && !is.null(names(xs))) {
        ii <- intersect(names(x0), names(xs))
        x0[ii] <- xs[ii]
      } else if (length(xs) == length(x0)) {
        x0[] <- xs
      }
      par0[[nm]] <- x0
      next
    }

    # Matrices: align by dimnames when possible, else recycle shape if equal
    if (is.matrix(x0)) {
      rnm <- rownames(x0); cnm <- colnames(x0)
      if (!is.null(dimnames(x0)) && !is.null(dimnames(xs))) {
        rr <- intersect(rnm, rownames(xs))
        cc <- intersect(cnm, colnames(xs))
        if (length(rr) && length(cc)) x0[rr, cc] <- xs[rr, cc]
      } else if (all(dim(xs) == dim(x0))) {
        x0[,] <- xs
      }
      par0[[nm]] <- x0
      next
    }
  }
  par0
}


#' Fit a Tiny Assessment Model (TAM)
#'
#' @description
#' Builds data with [make_dat()], initializes parameters with [make_par()],
#' constructs the RTMB objective, optimizes it, and returns a fitted object with
#' reports and standard errors.
#'
#' @details
#' Random-effect blocks are chosen automatically from the model settings:
#'
#' - Always includes `log_f` and `log_r`.
#' - Includes `missing` if there are missing observations.
#' - Includes `log_n` if `N_settings$process != "off"`.
#' - Includes `log_m` if `M_settings$process != "off"`.
#'
#' A warning is issued if the number of random effects exceeds 1.5 times the
#' number of observed data points (rough identifiability check).
#'
#' @param obs A named list of tidy observation tables (e.g., `catch`, `index`,
#'   `weight`, `maturity`). See [cod_obs] for an example.
#' @param interval Level in `(0, 1)` to use to generate confidence intervals,
#'   where applicable; default `0.95.`
#' @param silent Logical; if `TRUE`, disables tracing information.
#' @param add_osa_res Logical; add one-step-ahead residuals? Hard wired to
#'   apply the `"oneStepGaussianOffMode"` method. See [RTMB::oneStepPredict()]
#'   for details.
#' @param start_par Optional list of starting parameter values, matching the
#'   structure returned by [make_par()]. Useful for warm starts or retrospective
#'   runs where parameter estimates from a previous fit are reused to improve
#'   convergence and speed. Non-matching or missing entries are ignored.
#' @param grad_tol Numeric tolerance passed to [check_convergence()] when
#'   evaluating the fitted object's gradients. Defaults to `1e-2`.
#' @inheritDotParams make_dat
#'
#' @return
#' A list with components:
#'
#' - **call**: matched call.
#' - **dat**: data list returned by [make_dat()].
#' - **obj**: RTMB `ADFun` object.
#' - **opt**: `[stats::nlminb()]` optimization result.
#' - **rep**: list from `obj$report()`.
#' - **sdrep**: [RTMB::sdreport()] result.
#' - **fixed_par**: fixed parameter estimates in a tidy format (see [tidy_par()]).
#' - **random_par**: list of random parameter estimates in a tidy format (see
#'                   [tidy_par()]).
#' - **is_converged**: Did the model converge? (see [check_convergence()])
#' - **obs_pred**: `obs$catch` and `obs$index` data augmented with
#'                  predicted values, parameter estimates, and
#'                  standardized residuals (see [tidy_obs_pred()]).
#' - **pop**: A collection of population summaries in tidy format (see
#'            [tidy_pop()]).
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' fit$sdrep
#'
#' ## Fit with projections (status quo F)
#' fit2 <- update(fit,
#'   proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)
#' )
#'
#' if (interactive()) {
#'   fits <- list("without_proj" = fit, "with_proj" = fit2)
#'   vis_tam(fits)
#' }
#'
#' @importFrom stats nlminb rnorm
#'
#' @seealso
#' [make_dat()], [make_par()], [sim_tam()], [fit_retro()],
#' [RTMB::MakeADFun()], [RTMB::sdreport()]
#' @export
fit_tam <- function(
    obs,
    interval = 0.95,
    add_osa_res = FALSE,
    silent = FALSE,
    start_par = NULL,
    grad_tol = 1e-2,
    ...
) {

  call <- match.call()

  dat <- make_dat(obs, ...)
  par <- make_par(dat)
  if (!is.null(start_par)) {
    par <- .merge_start_par(par, start_par)
  }

  ran <- c("log_f", "log_r")
  if (dat$any_fill_missing) {
    ran <- c(ran, "missing")
  }
  if (dat$N_settings$process != "off") {
    ran <- c(ran, "log_n")
  }
  if (dat$M_settings$process != "off") {
    ran <- c(ran, "log_m")
  }
  map <- list()
  if (dat$M_settings$process == "ar1" && nlevels(dat$M_settings$age_blocks) == 1) {
    par$logit_phi_m[1] <- qlogis(0)
    map$logit_phi_m <- factor(c(NA, 1)) # phi_age moot when only one age block
  }

  make_nll_fun <- function(f, d) function(p) f(p, d) # use closure to avoid global assignment of data
  obj <- RTMB::MakeADFun(
    make_nll_fun(nll_fun, dat),
    par,
    random = ran,
    map = map,
    silent = silent
  )

  opt <- try(stats::nlminb(
    obj$par, obj$fn, obj$gr,
    control = list(eval.max = 1000, iter.max = 1000)
  ))
  rep <- obj$report()
  sdrep <- RTMB::sdreport(obj, getJointPrecision = TRUE)

  out <- list(
    call = call,
    dat = dat,
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep
  )

  class(out) <- c("tam_fit", "list")

  par_tabs <- tidy_par(out, interval = interval)
  out$fixed_par <- par_tabs$fixed
  out$random_par <- par_tabs$random
  out$obs_pred <- tidy_obs_pred(out, add_osa_res = add_osa_res, trace = !silent)
  out$pop <- tidy_pop(out, interval = interval)
  out$is_converged <- check_convergence(out, grad_tol = grad_tol, quiet = TRUE)
  out$grad_tol <- grad_tol

  .new_tam_fit(out)

}


#' Run a retrospective (peel) analysis
#'
#' @title Retrospective fits for TAM
#'
#' @description
#' Creates a sequence of terminal-year peels and refits the model on each
#' truncated year range.
#'
#' @details
#' Peel years are `(max_year - folds) : max_year`.
#' Each refit is attempted with `try()` so individual failures do not stop the sequence.
#' Refits are generated via `update(fit, years = ...)`; ensure your `fit` object
#' supports `update()` with a `years` argument. This function uses
#' [furrr::future_map()] to run the retros in parallel. Remember to plan your session (e.g.,
#' `future::plan(multisession, workers = 4)`).
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param folds Integer; number of terminal peels (default `2`).
#' @param start_from_fit Logical; use terminal parameter estimates as starting
#'    values for retrospective fits?
#' @param hindcast Logical; fit one-step-ahead projections? If `TRUE`,
#'    `proj_settings` will be set to `list(n_proj = 1, n_mean = 1, F_mult = 1)`
#'    to generate a one year status-quo F projection.
#' @param grad_tol Numeric tolerance for `max|grad|`. If `NULL` (default),
#'   the tolerance stored on `fit` (from [fit_tam()]) is used; otherwise the
#'   supplied value is passed to [check_convergence()]. Output from retro fits
#'   that exceed this tolerance are dropped (see [check_convergence()]).
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
#'     from catch and index observations using the [compute_hindcast_rmse()]
#'     function); and,
#'   - `fits` - refitted objects for each peel (named by terminal year).
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
    globals = NULL
) {
  .require_tam_fit(fit, arg = "fit")

  if (is.null(grad_tol)) {
    grad_tol <- if (!is.null(fit$grad_tol)) fit$grad_tol else 1e-3
  }

  min_year <- min(fit$dat$years)
  max_year <- max(fit$dat$years)
  retro_years <- seq(max_year - folds, max_year)
  if (start_from_fit) {
    start_par <- as.list(fit$sdrep, "Estimate")
  } else {
    start_par <- NULL
  }

  if (hindcast) {
    fit$call$proj_settings <- list(n_proj = 1, n_mean = 1, F_mult = 1)
  }

  progressr::with_progress({
    update_progress <- progressr::progressor(steps = length(retro_years))
    retro <- furrr::future_map(seq_along(retro_years), function(i) {
      if (retro_years[i] == max_year) {
        r <- fit # need not re-run terminal year
      } else {
        r <- suppressWarnings(
          try(
            stats::update(
              fit,
              years = min_year:retro_years[i],
              start_par = start_par,
              silent = TRUE
            ),
            silent = TRUE
          )
        )
      }
      if (inherits(r, "try-error")) {
        update_progress()
        return(r)
      }
      r$is_converged <- check_convergence(r, grad_tol = grad_tol)
      update_progress()
      r
    }, .options = furrr::furrr_options(seed = 1, packages = "tinyAM", globals = globals))
  }, enable = progress)
  names(retro) <- retro_years

  converged <- sapply(retro, function(x) !inherits(x, "try-error") && isTRUE(x$is_converged))

  if (any(!converged)) {
    cli::cli_warn(
      paste0("Model may not have converged for the following folds: ",
             paste(retro_years[!converged], collapse = ", "))
    )
  }

  fits <- retro[converged]
  if (length(fits) == 0L) cli::cli_abort("All folds failed convergence checks")

  out <- c(tidy_tam(model_list = fits, label = "fold"),
           list(fits = fits))

  rhos <- lapply(names(out$pop), function(nm) {
    d <- out$pop[[nm]]
    if ("age" %in% names(d)) {
      split_d <- split(d, d$age)
      rhos <- sapply(split_d, compute_mohns_rho)
      data.frame(metric = nm, age = as.numeric(names(rhos)), rho = rhos)
    } else {
      rho <- compute_mohns_rho(d)
      data.frame(metric = nm, age = NA, rho = rho)
    }
  })
  rhos <- do.call(rbind, rhos)

  out$mohns_rho <- rhos

  if (hindcast) {
    cols <- c("year", "age", "obs", "pred", "fold", "is_proj")
    d <- rbind(out$obs_pred$catch[, cols],
               out$obs_pred$index[, cols])
    out$hindcast_rmse <- compute_hindcast_rmse(d)
  }

  out
}


#' Run one-step-ahead hindcasts (convenience wrapper)
#'
#' @description
#' Convenience alias for [fit_retro()] with `hindcast = TRUE`. This performs
#' retrospective peels and, for each peel, fits a **one-step-ahead projection**
#' (status-quo F: `n_proj = 1`, `n_mean = 1`, `F_mult = 1`), and returns the
#' same structure as [fit_retro()] including `hindcast_rmse`.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param ... Passed through to [fit_retro()] (e.g., `folds`, `grad_tol`,
#'   `progress`, `globals`).
#'
#' @return
#' Exactly the return value of [fit_retro()] with `hindcast = TRUE`—i.e.,
#' a list containing `obs_pred`, `pop`, optional `mohns_rho`, `hindcast_rmse`,
#' and `fits` (see [fit_retro()] for full details).
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' \donttest{
#' if (interactive()) {
#'   # Choose your parallel plan (set once per session)
#'   future::plan(future::multisession, workers = 4)
#'   hc  <- fit_hindcast(fit, folds = 5, progress = TRUE)
#'   hc$hindcast_rmse
#'
#'   vis_tam(hc$fits)
#' }
#' }
#'
#' @seealso [fit_retro()]
#' @export
fit_hindcast <- function(fit, ...) {
  .require_tam_fit(fit, arg = "fit")

  fit_retro(fit, hindcast = TRUE, ...)
}



#' Quick convergence check for a TAM fit
#'
#' @description
#' Checks two basics and returns `TRUE` only if all pass:
#' (1) maximum absolute gradient from `sdreport`,
#' (2) Hessian positive-definite flag.
#'
#' If all pass, a short success message is printed unless `quiet = TRUE`.
#' If any check fails, a warning is emitted (not suppressed by `quiet`).
#'
#' @param fit A fitted TAM object containing `$sdrep`.
#' @param grad_tol Numeric tolerance for `max|grad|`. Default `1e-3`.
#' @param quiet Logical; if `TRUE` (default) suppresses the success message.
#'
#' @return Logical: `TRUE` if all checks pass, otherwise `FALSE`.
#' @importFrom cli cli_inform cli_warn format_warning
#' @export
check_convergence <- function(fit, grad_tol = 1e-3, quiet = TRUE) {
  sdrep <- NULL

  if (.is_tam_fit(fit)) {
    sdrep <- fit$sdrep
  } else if (inherits(fit, "sdreport")) {
    sdrep <- fit
  } else if (is.list(fit) && !is.null(fit$sdrep)) {
    sdrep <- fit$sdrep
  }

  if (is.null(sdrep)) {
    cli::cli_abort(c(
      "`{.arg fit}` must be either a {.cls tam_fit}, an {.cls sdreport},",
      "or a list containing an {.field sdrep} element."
    ))
  }

  if (inherits(sdrep, "sdreport")) {
    grad <- sdrep$gradient.fixed
    pd_hess <- sdrep$pdHess
  } else if (is.list(sdrep)) {
    grad <- sdrep$gradient.fixed
    pd_hess <- sdrep$pdHess
  } else {
    cli::cli_abort(c(
      "`{.arg fit}$sdrep` must be either an {.cls sdreport} or a list",
      "with components {.field gradient.fixed} and {.field pdHess}."
    ))
  }

  if (is.null(grad)) {
    cli::cli_abort("`{.arg fit}` must provide a gradient via `sdrep$gradient.fixed`.")
  }

  if (is.null(pd_hess)) {
    cli::cli_abort("`{.arg fit}` must provide a Hessian flag via `sdrep$pdHess`.")
  }

  max_grad <- max(abs(grad))
  grad_ok  <- is.finite(max_grad) && max_grad <= grad_tol
  hess_ok  <- isTRUE(pd_hess)
  ok       <- grad_ok && hess_ok

  main_text <- if (ok) "{.strong Model converged}" else "{.strong Model may not have converged}"
  grad_text <- sprintf("Maximum gradient [%s] %s tolerance [%s])",
                       signif(max_grad, 1),
                       if (grad_ok) "<=" else ">",
                       grad_tol)
  hess_text <- sprintf("Hessian %s positive definite",
                       if (hess_ok) "was" else "was not")
  grad_bullet <- if (grad_ok) "v" else "x"
  hess_bullet <- if (hess_ok) "v" else "x"
  bullets <- c(main_text, grad_text, hess_text)
  names(bullets) <- c("", grad_bullet, hess_bullet)

  if (ok) {
    if (!quiet) {
      cli::cli_inform(message = bullets)
    }
  } else {
    cli::cli_warn(message = bullets)
  }

  ok
}



