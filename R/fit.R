
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
#' This function assigns the constructed data list to a global (`dat <<- ...`)
#' to satisfy closures inside the likelihood (`nll_fun`).
#'
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
#' @param inputs A named list of tidy observation tables (e.g., `catch`, `index`,
#'   `weight`, `maturity`). See [northern_cod_data] for an example.
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
#' - **obs_pred**: `inputs$catch` and `inputs$index` data augmented with
#'                  predicted values, parameter estimates, and
#'                  standardized residuals (see [tidy_obs_pred()]).
#' - **pop**: A collection of population summaries in tidy format (see
#'            [tidy_pop()]).
#'
#' @examples
#' fit <- fit_tam(
#'   northern_cod_data,
#'   years = 1983:2024, ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw"),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' fit$sdrep
#'
#' @seealso
#' [make_dat()], [make_par()], [sim_tam()], [fit_retro()],
#' [RTMB::MakeADFun()], [RTMB::sdreport()]
#' @export
fit_tam <- function(inputs, silent = FALSE, ...) {

  call <- match.call()

  dat <<- make_dat(inputs, ...) # dat not found without global assignment
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

  n_ran <- length(unlist(par[ran]))
  n_obs <- sum(!is.na(dat$log_obs))
  if (n_ran > (n_obs * 1.5)) {
    warning(sprintf("Number of random effects (%d) exceed 1.5 times the number of observations (%d). Consider simplifying your model.", n_ran, n_obs))
  }

  obj <- MakeADFun(
    nll_fun,
    par,
    random = ran,
    silent = silent
  )

  opt <- nlminb(
    obj$par, obj$fn, obj$gr,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  opt$objective
  rep <- obj$report()
  sdrep <- sdreport(obj)

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
#' truncated year range, attaching the list of refits to `fit$retro`.
#'
#' @details
#' Peel years are `(max_year - folds) : max_year`.
#' Each refit is attempted with `try()` so individual failures do not stop the sequence.
#' Refits are generated via `update(fit, years = ...)`; ensure your `fit` object
#' supports `update()` with a `years` argument.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param folds Integer; number of terminal peels (default `2`).
#'
#' @return
#' The input `fit` with an added `$retro` list whose elements are the
#' refitted objects for each peel (named by terminal year).
#'
#' @examples
#' \dontrun{
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' retro_fit <- fit_retro(fit, folds = 5)
#' lapply(retro_fit$retro, function(x) x$rep$ssb)
#' }
#'
#' @seealso
#' [fit_tam()], [sim_tam()]
#' @export
fit_retro <- function(fit, folds = 2) {

  min_year <- min(fit$dat$years)
  max_year <- max(fit$dat$years)
  retro_years <- (max_year - folds):(max_year)

  retro <- vector("list", length(retro_years))
  names(retro) <- retro_years
  for (i in seq_along(retro_years)) {
    r <- try(update(fit, years = min_year:retro_years[i]))
    if (inherits(r, "try-error") || !r$is_converged) {
      retro[[i]] <- r
    } else {
      retro[[i]] <- NULL
      warning("Model failed when terminal year was ", retro_years[i])
    }
  }

  fit$retro <- retro
  fit

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
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' rep_sim <- sim_tam(fit, obs_only = FALSE)
#' str(rep_sim$log_obs)
#'
#' @seealso
#' [fit_tam()], [fit_retro()]
#'
#' @export
sim_tam <- function(fit, obs_only = FALSE) {

  obj <- fit$obj
  p <- as.list(fit$sdrep, "Estimate")
  sims <- nll_fun(p, simulate = TRUE)
  rep <- obj$report()
  if (!obs_only) {
    p[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(p, simulate = TRUE)
    rep <- obj$report(unlist(p))
  }
  rep$log_obs <- sims$log_obs
  rep

}


#' Quick convergence check for a TAM fit
#'
#' @description
#' Checks three basics and returns `TRUE` only if all pass:
#' (1) maximum absolute gradient from `sdreport`,
#' (2) Hessian positive-definite flag, and
#' (3) optimizer convergence code (`0` for `nlminb()`).
#'
#' If all pass, a short success message is printed unless `quiet = TRUE`.
#' If any check fails, a warning is emitted (not suppressed by `quiet`).
#'
#' @param fit A fitted TAM object containing `$sdrep` and `$opt`.
#' @param grad_tol Numeric tolerance for `max|grad|`. Default `1e-3`.
#' @param quiet Logical; if `TRUE` (default) suppresses the success message.
#'
#' @return Logical: `TRUE` if all checks pass, otherwise `FALSE`.
#' @export
check_convergence <- function(fit, grad_tol = 1e-3, quiet = TRUE) {
  fmt_num <- function(x) if (is.finite(x)) signif(x, 3) else "NA"

  g        <- tryCatch(fit$sdrep$gradient.fixed, error = function(e) NULL)
  max_grad <- if (length(g)) max(abs(g), na.rm = TRUE) else NA_real_
  pdHess   <- tryCatch(fit$sdrep$pdHess,        error = function(e) NA)
  conv     <- tryCatch(fit$opt$convergence,     error = function(e) NA_integer_)

  grad_ok  <- is.finite(max_grad) && max_grad <= grad_tol
  hess_ok  <- isTRUE(pdHess)
  optim_ok <- isTRUE(conv == 0)
  ok       <- grad_ok && hess_ok && optim_ok

  issues <- c(
    if (!grad_ok) sprintf("max|grad|=%s > tol=%g", fmt_num(max_grad), grad_tol),
    if (!hess_ok) "Hessian not PD",
    if (!optim_ok) sprintf("opt code=%s != 0", as.character(conv))
  )
  msg <- sprintf(
    "Convergence %s: max|grad|=%s (tol=%g), pdHess=%s, opt code=%s%s",
    if (ok) "OK" else "FAIL",
    fmt_num(max_grad), grad_tol, as.character(pdHess), as.character(conv),
    if (length(issues) && !ok) paste0(" | issues: ", paste(issues, collapse = "; ")) else ""
  )

  if (ok) {
    if (!quiet) message(msg)
  } else {
    warning(msg, call. = FALSE)
  }
  ok
}
