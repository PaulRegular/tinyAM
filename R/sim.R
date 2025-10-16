
#' Draw MLE parameter list (no uncertainty)
#'
#' @keywords internal
#' @noRd
.draw_none <- function(fit) {
  as.list(fit$sdrep, "Estimate")
}

#' Draw fixed effects from MVN(sdreport) and map back to par list
#'
#' @details Uses the asymptotic MVN with mean `sdrep$par.fixed` and covariance
#' `sdrep$cov.fixed` on the estimation scale; inserts the draw into the full
#' parameter vector via `obj$env$parList()`.
#'
#' @return A parameter list suitable for `MakeADFun`/`nll_fun`.
#' @keywords internal
#' @noRd
#' @importFrom MASS mvrnorm
.draw_fixed <- function(fit) {
  sdr <- fit$sdrep
  mu <- sdr$par.fixed
  V  <- sdr$cov.fixed
  fixed <- as.numeric(MASS::mvrnorm(1L, mu = mu, Sigma = V))
  par <- fit$obj$env$last.par.best
  fit$obj$env$parList(x = fixed, par = par)
}

#' Create a fast joint (random + fixed) drawer using the Laplace joint precision
#'
#' @details Pre-computes a sparse Cholesky of the joint precision. The returned function
#' draws `par_full ~ N(last.par.best, Q^{-1})` via two sparse triangular solves
#' and applies the permutation from the factorization. Much of this code is based on
#' the `sparseMVN::rmvn-sparse()` function.
#'
#' @return A function `function(fit) -> par_list` that generates a joint draw.
#' @keywords internal
#' @noRd
#' @importFrom Matrix Cholesky expand solve
.make_draw_joint <- function(fit) {
  sdr <- fit$sdrep
  if (is.null(sdr$jointPrecision)) {
    sdr <- RTMB::sdreport(fit$obj, getJointPrecision = TRUE)
  }
  Q  <- sdr$jointPrecision
  mu <- fit$obj$env$last.par.best

  CH <- Matrix::Cholesky(Q, LDL = FALSE)
  A  <- Matrix::expand(CH)  # contains $L and permutation matrix $P

  function(fit) {
    z <- rnorm(nrow(Q))
    y <- Matrix::solve(Matrix::t(A$L), z)     # L' y = z
    y <- Matrix::crossprod(A$P, y)            # P' y
    par_full <- as.numeric(mu + y)            # N(mu, Q^{-1})
    fit$obj$env$parList(par = par_full)
  }
}

#' Simulate observations (and optionally re-draw random effects) given a par-drawer
#'
#' @param fit Fitted TAM object.
#' @param par_fun Function taking `fit` and returning a parameter list
#'   (e.g., `.draw_none`, `.draw_fixed`, or the result of `.make_draw_joint(fit)`).
#' @param redraw_random Logical; if `TRUE`, replace the parameter list's random
#'   effects with newly simulated fields from the process before recomputing.
#'
#' @return Flat list of data frames (obs + population summaries).
#' @keywords internal
#' @noRd
.sim_obs <- function(fit, par_fun = NULL, redraw_random = FALSE) {
  obj <- fit$obj
  dat <- fit$dat
  par <- par_fun(fit)

  # Draw obs | par
  sims <- nll_fun(par, dat, simulate = TRUE)

  # Optionally use simulated random effects
  if (redraw_random) {
    par[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(par, dat, simulate = TRUE)
  }

  # Rebuild report with simulated obs (+ maybe simulated RE)
  dat$log_obs <- sims$log_obs
  make_nll_fun <- function(f, d) function(p) f(p, d)
  obj <- RTMB::MakeADFun(make_nll_fun(nll_fun, dat), par)
  rep <- obj$report()

  # Tidy population tables
  pop <- tidy_rep(list(dat = dat, rep = rep))

  # Tidy observations (catch/index) using the simulated log_obs split
  obs <- fit$dat$obs[c("catch", "index")]
  split_obs <- split(exp(sims$log_obs), fit$dat$obs_map$type)
  obs$catch$obs  <- split_obs$catch
  obs$index$obs  <- split_obs$index

  # Return a flat list of data.frames
  c(obs, pop)
}


#' Simulate from a fitted TAM
#'
#' @description
#' Runs the TAM likelihood in simulation mode to generate synthetic observations,
#' and optionally random-effect fields, then recomputes reported quantities under
#' those draws. Results are returned as tidy data frames stacked across `n`
#' simulations with a `sim = 1..n` column.
#'
#' @details
#' The simulation has two orthogonal controls:
#'
#' - **Parameter uncertainty** via `par_uncertainty`:
#'   - `"none"`  — use point estimates `(û, θ̂)`.
#'   - `"fixed"` — sample **fixed effects** `θ ~ MVN(sdrep$par.fixed, sdrep$cov.fixed)`.
#'   - `"joint"` — sample **(random + fixed)** jointly from the Laplace
#'     approximate posterior using the **joint precision** (sparse Cholesky).
#'
#' - **Random-effect handling** via `redraw_random`:
#'   - `FALSE` — keep the sampled/fitted random effects and simulate **observations only**
#'     (posterior-predictive when `par_uncertainty = "joint"`).
#'   - `TRUE`  — generate **new process fields** for the random effects and re-simulate
#'     (projection/HCR style prior-predictive runs).
#'
#' Parallel execution is supported via [furrr::future_map()]. Call
#' `future::plan()` beforehand if you want parallel workers.
#'
#' @param fit A fitted TAM object returned by [fit_tam()].
#' @param n Integer; number of simulations (default `10`).
#' @param par_uncertainty Character; one of `"joint"`, `"fixed"`, `"none"`.
#'   Controls how the parameter list is sampled before each simulation (see Details).
#' @param redraw_random Logical; if `TRUE`, re-draw random-effect fields from their
#'   process models on each run (recommended for projections). If `FALSE`, keep
#'   random effects and simulate observations only.
#' @param progress Logical; show a progress bar using [progressr::with_progress()]
#'   (default `TRUE`).
#' @param globals Optional character vector of global objects for parallel workers
#'   (forwarded to [furrr::furrr_options()]).
#' @param seed Seed passed to [furrr::furrr_options()] for reproducibility.
#'
#' @return
#' A **named list of data frames** (e.g., `catch`, `index`, `N`, `F`, `M`, `Z`, `ssb`, …),
#' each stacked across `n` simulations with a `sim` column.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   # Set-up parallel workers and fit model
#'   future::plan(future::multisession, workers = 4)
#'   fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#'
#'   # Draw fixed-effects uncertainty and redraw random process fields
#'   sims1 <- sim_tam(fit, n = 10, par_uncertainty = "fixed", redraw_random = TRUE)
#'
#'   # Joint draw of random and fixed effects, keep random effects, simulate obs only
#'   sims2 <- sim_tam(fit, n = 10, par_uncertainty = "joint", redraw_random = FALSE)
#'
#'   # Conditional/fast: point estimates, obs only
#'   sims3 <- sim_tam(fit, n = 10, par_uncertainty = "none", redraw_random = FALSE)
#' }
#' }
#'
#' @seealso [fit_tam()], [tidy_rep()], [stack_nested()]
#' @export
sim_tam <- function(
    fit,
    n = 10,
    par_uncertainty = c("joint", "fixed", "none"),
    redraw_random = FALSE,
    progress = TRUE,
    globals = NULL,
    seed = TRUE
) {
  match.arg(par_uncertainty)

  draw_par <- switch(
    par_uncertainty,
    none  = .draw_none,
    fixed = .draw_fixed,
    joint = .make_draw_joint(fit)
  )

  sims <- progressr::with_progress({
    update_progress <- progressr::progressor(steps = n)
    furrr::future_map(
      seq_len(n),
      function(i) {
        res <- .sim_obs(fit, par_fun = draw_par, redraw_random = redraw_random)
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
  stack_nested(sims, label = "sim")
}
