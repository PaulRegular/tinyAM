

.one_sim <- function(fit, obs_only = FALSE, sim_number = 1) {

  ## Simulate observations | par
  obj <- fit$obj
  dat <- sim_dat <- fit$dat
  par <- sim_par <- as.list(fit$sdrep, "Estimate")
  sims <- nll_fun(par, dat, simulate = TRUE)

  ## Use simulated random effects?
  if (!obs_only) {
    sim_par[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(sim_par, dat, simulate = TRUE)
  }

  ## Update report using sim_dat and sim_par
  sim_dat$log_obs <- sims$log_obs
  make_nll_fun <- function(f, d) function(p) f(p, d)
  sim_obj <- RTMB::MakeADFun(make_nll_fun(nll_fun, sim_dat), sim_par)
  sim_rep <- sim_obj$report()

  ## Tidy, label, and return
  sim_pop <- tidy_rep(list(dat = sim_dat, rep = sim_rep))
  for (nm in names(sim_pop)) sim_pop[[nm]]$sim <- sim_number
  sim_obs <- fit$dat$obs[c("catch", "index")]
  split_obs <- split(exp(sims$log_obs), fit$dat$obs_map$type)
  sim_obs$catch$obs <- split_obs$catch
  sim_obs$catch$sim <- sim_number
  sim_obs$index$obs <- split_obs$index
  sim_obs$index$sim <- sim_number

  list(sim_obs = sim_obs, sim_pop = sim_pop)

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
#' @inheritParams fit_retro
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
sim_tam <- function(fit, n = 10, obs_only = FALSE,
                    progress = TRUE,
                    globals = NULL) {

  progressr::with_progress({
    update_progress <- progressr::progressor(steps = n)
    sims <- furrr::future_map(seq.int(n), function(i) {
      sim <- .one_sim(fit, obs_only = obs_only, sim_number = i)
      update_progress()
      sim
    }, .options = furrr::furrr_options(seed = 1, packages = "tinyAM", globals = globals))
  }, enable = progress)
  names(sims) <- seq.int(n)

  pop_nms <- names(sims[[1]]$sim_pop)
  lapply(pop_nms, function(nm) sims)

}

