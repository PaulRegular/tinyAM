
.one_sim <- function(fit, obs_only = FALSE) {
  obj <- fit$obj
  dat <- sim_dat <- fit$dat
  par <- sim_par <- as.list(fit$sdrep, "Estimate")

  # Draw obs | par
  sims <- nll_fun(par, dat, simulate = TRUE)

  # Optionally use simulated random effects
  if (!obs_only) {
    sim_par[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(sim_par, dat, simulate = TRUE)
  }

  # Rebuild report with simulated obs (+ maybe simulated RE)
  sim_dat$log_obs <- sims$log_obs
  make_nll_fun <- function(f, d) function(p) f(p, d)
  sim_obj <- RTMB::MakeADFun(make_nll_fun(nll_fun, sim_dat), sim_par)
  sim_rep <- sim_obj$report()

  # Tidy population tables
  sim_pop <- tidy_rep(list(dat = sim_dat, rep = sim_rep))

  # Tidy observations (catch/index) using the simulated log_obs split
  sim_obs <- fit$dat$obs[c("catch", "index")]
  split_obs <- split(exp(sims$log_obs), fit$dat$obs_map$type)
  sim_obs$catch$obs  <- split_obs$catch
  sim_obs$index$obs  <- split_obs$index

  # Return a flat list of data.frames
  c(sim_obs, sim_pop)
}


#' Simulate from a fitted TAM
#'
#' @description
#' Draws simulated observations (and optionally random effects) from a fitted model
#' by calling the likelihood with `simulate = TRUE`, then recomputes reported quantities
#' under those draws. Results are returned as tidy, stacked tables across `n`
#' simulations (one column `sim = 1..n` is added to each table).
#'
#' @details
#' For each simulation:
#' 1. Extract MLEs: `p <- as.list(fit$sdrep, "Estimate")`.
#' 2. Draw `log_obs` (and, if `obs_only = FALSE`, random effects) using `nll_fun(p, dat, simulate = TRUE)`.
#' 3. Replace `dat$log_obs` with simulated values, rebuild a small RTMB object, and call `report()`.
#' 4. Return:
#'    - observation tables for `catch` and `index` with simulated `obs`, and
#'    - population tables from `report` via [tidy_rep()].
#'
#' All simulation results are then stacked with a `sim` column via [stack_nested()].
#'
#' This function uses [furrr::future_map()] to run the retros in parallel. Remember to plan your
#' session (e.g., `future::plan(multisession, workers = 4)`).
#'
#' @param fit A fitted TAM object returned by [fit_tam()].
#' @param n Integer; number of simulations (default `10`).
#' @param obs_only Logical; if `TRUE` (default), simulate observations only (random effects fixed at MLEs).
#'   If `FALSE`, also simulate random effects and recompute reported quantities under them.
#' @param progress Logical; show a progress bar using [progressr::with_progress()]. Default `TRUE`.
#' @param globals Optional character vector of global objects for parallel workers (forwarded to `furrr`).
#' @param seed Seed for random number generation forwarded to `furrr` (see [furrr::furrr_options()] for details).
#'
#' @return
#' A **named list of data frames** (e.g., `catch`, `index`, `N`, `F`, `M`, `Z`, `ssb`, ...),
#' where each data frame is **row-bound across simulations** and contains a `sim` column (1..n).
#'
#' @examples
#' \dontrun{
#' future::plan(multisession, workers = 4)
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' sims <- sim_tam(fit, n = 5, obs_only = TRUE)
#' head(sims$ssb)     # has column 'sim'
#' head(sims$catch)   # simulated obs for catch with 'sim'
#' }
#'
#' @seealso [fit_tam()], [tidy_rep()], [stack_nested()]
#' @export
sim_tam <- function(
    fit,
    n = 10,
    obs_only = TRUE,
    progress = TRUE,
    globals = NULL,
    seed = TRUE
) {
  sims <- progressr::with_progress({
    update_progress <- progressr::progressor(steps = n)
    furrr::future_map(
      seq_len(n),
      function(i) {
        res <- .one_sim(fit, obs_only = obs_only)
        update_progress()
        res
      },
      .options = furrr::furrr_options(
        seed = seed,
        packages = "tinyAM",
        globals = globals
      )
    )
  }, enable = progress)

  names(sims) <- as.character(seq_len(n))
  stack_nested(sims, id_col = "sim")
}
