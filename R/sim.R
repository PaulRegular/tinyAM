
#' Draw a fixed-parameter vector from sdreport() and map back to par list
#'
#' @keywords internal
#' @importFrom MASS mvrnorm
#' @noRd
.samp_fixed <- function(fit) {
  sdr <- fit$sdrep
  mu <- sdr$par.fixed
  V  <- sdr$cov.fixed
  sim_fixed <- as.numeric(MASS::mvrnorm(1L, mu = mu, Sigma = V))
  all_par <- fit$obj$env$last.par.best
  fit$obj$env$parList(x = sim_fixed, par = all_par)
}

#' Simulate observations (optionally sample fixed and random effects)
#'
#' @keywords internal
#' @noRd
.sim_obs <- function(fit, samp_fixed = TRUE, samp_random = FALSE) {
  obj <- fit$obj
  sdat <- fit$dat
  if (samp_fixed) {
    spar <- .samp_fixed(fit)
  } else {
    spar <- as.list(fit$sdrep, "Estimate")
  }

  # Draw obs | par
  sims <- nll_fun(spar, sdat, simulate = TRUE)

  # Optionally use simulated random effects
  if (samp_random) {
    spar[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(spar, sdat, simulate = TRUE)
  }

  # Rebuild report with simulated obs (+ maybe simulated RE)
  sdat$log_obs <- sims$log_obs
  make_nll_fun <- function(f, d) function(p) f(p, d)
  sobj <- RTMB::MakeADFun(make_nll_fun(nll_fun, sdat), spar)
  srep <- sobj$report()

  # Tidy population tables
  spop <- tidy_rep(list(dat = sdat, rep = srep))

  # Tidy observations (catch/index) using the simulated log_obs split
  sobs <- fit$dat$obs[c("catch", "index")]
  split_obs <- split(exp(sims$log_obs), fit$dat$obs_map$type)
  sobs$catch$obs  <- split_obs$catch
  sobs$index$obs  <- split_obs$index

  # Return a flat list of data.frames
  c(sobs, spop)
}


#' Simulate from a fitted TAM
#'
#' @description
#' Draws simulated observations (and optionally random effects) from a fitted
#' model using [nll_fun()] in simulation mode. Each simulation regenerates
#' observations — and, optionally, fixed and random effects — and recomputes all reported
#' quantities. The results are returned as tidy data frames stacked across `n`
#' simulations, with a column `sim = 1..n`.
#'
#' @details
#' For each simulation:
#' 1. If `samp_fixed = TRUE`, draw fixed-effect parameters from the multivariate
#'    normal distribution defined by `sdreport()` (i.e., parameter uncertainty).
#' 2. Call [nll_fun()] with `simulate = TRUE` to generate new observations.
#' 3. If `samp_random = TRUE`, also simulate new random-effect fields
#'    (e.g., F, N, M deviations) from their process models and recompute
#'    reported quantities under those draws.
#' 4. Replace `dat$log_obs` with the simulated observations, rebuild a minimal
#'    RTMB object, call `report()`, and extract derived quantities via
#'    [tidy_rep()].
#'
#' The resulting observation and population tables are then row-bound with a
#' `sim` column using [stack_nested()].
#'
#' Parallel execution is supported via [furrr::future_map()]. Plan your
#' session first, e.g. `future::plan(multisession, workers = 4)`.
#'
#' @param fit A fitted TAM object returned by [fit_tam()].
#' @param n Integer; number of simulations (default `10`).
#' @param samp_fixed Logical; if `TRUE`, draw fixed-effect parameters from
#'   MVN(`sdrep$par.fixed`, `sdrep$cov.fixed`) before each simulation to
#'   include parameter uncertainty (default `TRUE`).
#' @param samp_random Logical; if `TRUE`, re-simulate random-effect fields from
#'   their process models on each run (default `TRUE`). Set `FALSE` to keep the
#'   fitted random effects and simulate only new observations.
#' @param progress Logical; show a progress bar with
#'   [progressr::with_progress()] (default `TRUE`).
#' @param globals Optional character vector of global objects for parallel
#'   workers (forwarded to [furrr::furrr_options()]).
#' @param seed Seed passed to [furrr::furrr_options()] for reproducibility.
#'
#' @return
#' A **named list of data frames** (e.g. `catch`, `index`, `N`, `F`, `M`,
#' `Z`, `ssb`, …), each containing results stacked across `n` simulations and
#' a `sim` column (1..n).
#'
#' @examples
#' \dontrun{
#' future::plan(multisession, workers = 4)
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#'
#' # Include parameter uncertainty and new process draws (default)
#' sims <- sim_tam(fit, n = 100, samp_fixed = TRUE, samp_random = TRUE)
#'
#' # Conditional simulations (MLE parameters, fixed REs)
#' sims0 <- sim_tam(fit, n = 100, samp_fixed = FALSE, samp_random = FALSE)
#' }
#'
#' @seealso [fit_tam()], [tidy_rep()], [stack_nested()]
#' @export
sim_tam <- function(
    fit,
    n = 10,
    samp_fixed = TRUE,
    samp_random = FALSE,
    progress = TRUE,
    globals = NULL,
    seed = TRUE
) {
  sims <- progressr::with_progress({
    update_progress <- progressr::progressor(steps = n)
    furrr::future_map(
      seq_len(n),
      function(i) {
        res <- .sim_obs(fit, samp_fixed = samp_fixed, samp_random = samp_random)
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
