
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

