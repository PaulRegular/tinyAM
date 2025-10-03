
#' Fit a Tiny Assessment Model (TAM)
#'
#' @title Fit TAM
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
#' @param interval Level in `(0, 1)` to use to generate confidence
#'                 intervals, where applicable; default `0.95 (see [tidy_tam()]).`
#' @param silent Logical; if `TRUE`, disables tracing information.
#' @inheritDotParams make_dat
#'
#' @return
#' A list with components:
#'
#' - **call**: matched call.
#' - **dat**: data list returned by [make_dat()].
#' - **obj**: RTMB `ADFun` object.
#' - **opt**: `nlminb` optimization result.
#' - **rep**: list from `obj$report()`.
#' - **sdrep**: [RTMB::sdreport()] result.
#' - **is_converged**: Did the model converge? (see [check_convergence()])
#' - **obs_pred**: `obs$catch` and `obs$index` data augmented with
#'                  predicted values, parameter estimates, and
#'                  standardized residuals (see [tidy_obs_pred()]).
#' - **pop**: A collection of population summaries in tidy format (see
#'            [tidy_pop()]).
#'
#' @examples
#' fit <- fit_tam(
#'   cod_obs,
#'   years = 1983:2024, ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw"),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' fit$sdrep
#'
#' ## Fit with projections (status quo F)
#' fit2 <- update(fit,
#'   proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)
#' )
#'
#' @seealso
#' [make_dat()], [make_par()], [sim_tam()], [fit_retro()],
#' [RTMB::MakeADFun()], [RTMB::sdreport()]
#' @export
fit_tam <- function(obs, interval = 0.95, silent = FALSE, ...) {

  call <- match.call()

  dat <- make_dat(obs, ...)
  par <- make_par(dat)

  ran <- c("log_f", "log_r")
  if (length(par$missing) > 0) {
    ran <- c(ran, "missing")
  }
  if (dat$N_settings$process != "off") {
    ran <- c(ran, "log_n")
  }
  if (dat$M_settings$process != "off") {
    ran <- c(ran, "log_m")
  }

  make_nll_fun <- function(f, d) function(p) f(p, d) # use closure to avoid global assignment of data
  obj <- RTMB::MakeADFun(
    make_nll_fun(nll_fun, dat),
    par,
    random = ran,
    silent = silent
  )

  opt <- try(nlminb(
    obj$par, obj$fn, obj$gr,
    control = list(eval.max = 1000, iter.max = 1000)
  ))
  rep <- obj$report()
  sdrep <- RTMB::sdreport(obj)

  out <- list(
    call = call,
    dat = dat,
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep
  )

  out$obs_pred <- tidy_obs_pred(out)
  out$pop <- tidy_pop(out)
  out$is_converged <- check_convergence(out, quiet = TRUE)

  out

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
#' @param grad_tol Numeric tolerance for `max|grad|`. Default `1e-3`. Output
#'                 from retro fits that exceed this tolerance are dropped (see
#'                 [check_convergence()]).
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
#'     function); and,
#'   - `fits` - refitted objects for each peel (named by terminal year).
#' The `obs_pred` and `pop` objects are created using [tidy_tam()]. A `fold` column
#' is added to each data.frame within these objects which specifies the terminal
#' year.
#'
#' @examples
#' \dontrun{
#' # Choose your parallel plan (set once per session)
#' future::plan(future::multisession, workers = 4)
#'
#' fit <- fit_tam(
#'   cod_obs,
#'   years = 1983:2024,
#'   ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~M_assumption),
#'   obs_settings = list(q_form = ~ q_block, sd_form = ~ sd_obs_block)
#' )
#' retros <- fit_retro(fit, folds = 5, progress = TRUE)
#' head(retros$pop$ssb)
#' head(retros$mohns_rho)
#' }
#'
#' @seealso
#' [fit_tam()], [sim_tam()]
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr with_progress progressor
#' @export
fit_retro <- function(fit, folds = 2, grad_tol = 1e-3, progress = TRUE, globals = NULL) {
  min_year <- min(fit$dat$years)
  max_year <- max(fit$dat$years)
  retro_years <- seq(max_year - folds, max_year)

  progressr::with_progress({
    update_progress <- progressr::progressor(steps = length(retro_years))
    retro <- furrr::future_map(seq_along(retro_years), function(i) {
      r <- suppressWarnings(
        try(update(fit, years = min_year:retro_years[i], silent = TRUE), silent = TRUE)
      )
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
      paste0("Model may not have converged for the following retro folds: ",
             paste(retro_years[!converged], collapse = ", "))
    )
  }

  fits <- retro[converged]
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

  out
}


#' Run one-step-ahead hindcasts
#'
#' @description
#' Like [fit_retro], this function creates a sequence of terminal-year peels and
#' refits the model on each truncated year range, except this function forces
#' the projection of one year for comparison of one-step-ahead projections relative
#' to actual observations to be used for hindcast validation.
#'
#' @inheritParams fit_retro
#' @inheritDotParams fit_retro
#'
#' @details
#' This is a convenience wrapper around [fit_retro()] that always performs
#' a **one-step-ahead projection** with `F_mult = 1` (status-quo F),
#' `n_proj = 1`, and `n_mean = 1`.
#'
#' @return
#' A list whose elements are:
#'   - `obs_pred` — a named list of stacked data frames (e.g., `catch`, `index`);
#'   - `pop` — a named list of stacked data frames (e.g., `ssb`, `N`, `M`,
#'   `mu_M`, `F`, `mu_F`, `Z`, `ssb_mat`);
#'   - `rmse` - root mean squared error for use as an overall measure of the forecast
#'     skill of the model (calculated from catch and index observations using the
#'     [compute_hindcast_rmse()] function); and,
#'   - `fits` - refitted objects for each peel (named by terminal year).
#' The `obs_pred` and `pop` objects are created using [tidy_tam()]. A `fold` column
#' is added to each data.frame within these objects which specifies the terminal
#' year.
#'
#' @seealso [fit_retro()]
#'
#' @examples
#' \dontrun{
#' # Choose your parallel plan (set once per session)
#' future::plan(future::multisession, workers = 4)
#'
#' fit1 <- fit_tam(
#'   cod_obs,
#'   years = 1983:2024,
#'   ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~M_assumption),
#'   obs_settings = list(q_form = ~ q_block, sd_form = ~ sd_obs_block)
#' )
#' fit2 <- update(fit1, N_settings = list(process = "ar1", init_N0 = FALSE))
#' hindcasts1 <- fit_hindcast(fit1, folds = 5, progress = TRUE)
#' hindcasts2 <- fit_hindcast(fit2, folds = 5, progress = TRUE)
#' hindcasts1$rmse
#' hindcasts2$rmse
#' }
#'
#' @export
fit_hindcast <- function(fit, ...) {
  fit$call$proj_settings <- list(n_proj = 1, n_mean = 1, F_mult = 1)
  out <- fit_retro(fit, ...)
  out$mohns_rho <- NULL
  cols <- c("year", "age", "obs", "pred", "fold", "is_proj")
  d <- rbind(out$obs_pred$catch[, cols],
             out$obs_pred$index[, cols])
  out$rmse <- compute_hindcast_rmse(d)
  out
}



#' Simulate from a fitted TAM
#'
#' @title Simulation using a fitted TAM
#'
#' @description
#' Simulates random effects and observations from a fitted model by calling the
#' likelihood with `simulate = TRUE`, then regenerates the REPORTed
#' quantities. Simulation is performed purely in R (not via RTMB `simref`).
#'
#' @details
#' The function:
#'
#' 1. Extracts MLEs as a parameter list: `as.list(fit$sdrep, "Estimate")`.
#' 2. Calls `nll_fun(p, simulate = TRUE)` to draw new random effects
#'    and `log_obs`.
#' 3. If `obs_only = FALSE` (default), injects the simulated random
#'    effects back into `p` (`p[obj$env$.random]`), calls
#'    `nll_fun` again (to update any derived quantities), and runs
#'    `obj$report(unlist(p))`.
#' 4. Returns the REPORT list with `log_obs` replaced by the simulated
#'    values; if `obs_only = FALSE`, all derived objects reflect the
#'    simulated random effects.
#'
#' @param fit A fitted TAM object returned by [fit_tam()].
#' @param obs_only Logical; if `TRUE`, only new `log_obs` values are
#'   simulated and inserted into the final report.
#'   If `FALSE` (default), simulated random effects are also injected and all REPORTed quantities are recomputed.
#'
#' @return
#' A list of REPORTed objects (same structure as `fit$rep`) with
#' `log_obs` replaced by simulated values.
#' If `obs_only = FALSE`, the report also reflects simulated random effects.
#'
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' rep_sim <- sim_tam(fit, obs_only = FALSE)
#' str(rep_sim$log_obs)
#'
#' @seealso
#' [fit_tam()], [fit_retro()]
#'
#' @export
sim_tam <- function(fit, obs_only = FALSE) {

  obj <- fit$obj
  dat <- fit$dat
  p <- as.list(fit$sdrep, "Estimate")
  sims <- nll_fun(p, dat, simulate = TRUE)
  rep <- obj$report()
  if (!obs_only) {
    p[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(p, dat, simulate = TRUE)
    rep <- obj$report(unlist(p))
  }
  rep$log_obs <- sims$log_obs
  rep

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
  max_grad <- max(abs(fit$sdrep$gradient.fixed))
  grad_ok  <- is.finite(max_grad) && max_grad <= grad_tol
  hess_ok  <- isTRUE(fit$sdrep$pdHess)
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
