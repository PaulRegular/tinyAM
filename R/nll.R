
#' 2D Gaussian Process Densities and Simulators (age × year)
#'
#' @title Two–dimensional process error: density and simulation
#'
#' @description
#' Helpers for working with simple age × year process-error fields, assuming a
#' separable 2D AR(1) structure.
#'
#' - [dprocess_2d()] returns the log-density contribution for a given matrix `x`.
#' - [rprocess_2d()] generates a matrix draw with the requested dependence structure.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n_y \times n_a}} denote the process on
#' year (\eqn{y}) and age (\eqn{a}) indices.
#' It follows a stationary separable AR(1) in both dimensions with correlations
#' \eqn{\phi_\text{age}} (columns) and \eqn{\phi_\text{year}} (rows), satisfying \eqn{|\phi|<1}.
#'
#' The implied covariance is
#' \deqn{\mathrm{Cov}\{\mathrm{vec}(X)\} =
#'       \frac{\sigma^2}{(1-\phi_\text{age}^2)(1-\phi_\text{year}^2)}
#'       \; \Sigma_\text{age} \otimes \Sigma_\text{year},}
#'
#' where \eqn{\Sigma_\cdot} are AR(1) correlation matrices with entries
#' \eqn{\phi^{|i-j|}}.
#'
#' @param x A numeric matrix (\eqn{n_y \times n_a}) of process residuals
#'   for [dprocess_2d()].
#' @param ny,na Positive integers: numbers of years and ages for
#'   [rprocess_2d()].
#' @param phi Length-2 numeric vector \code{c(phi_age, phi_year)} with
#'   values in \eqn{(0, 1)} for AR(1).
#' @param sd Positive scalar \eqn{\sigma}.
#'
#' @return
#' - [dprocess_2d()]: a single numeric log-density value.
#' - [rprocess_2d()]: a numeric matrix of dimension \eqn{n_y \times n_a}.
#'
#' @examples
#' # Simulate then calculate log-density
#' set.seed(1)
#' X_ar1 <- rprocess_2d(ny = 10, na = 8, sd = 0.3, phi = c(0.5, 0.9))
#' dprocess_2d(X_ar1, sd = 0.3, phi = c(0.5, 0.9))
#'
#' @seealso
#' [RTMB::dseparable()], [RTMB::dautoreg()], [MASS::mvrnorm()]
#'
#' @import RTMB
#' @rdname process_2d
#' @export
dprocess_2d <- function(x, phi = c(0, 0), sd = 1) {

  phi_age  <- phi[1]
  phi_year <- phi[2]
  fa <- function(z) dautoreg(z, phi = phi_age, log = TRUE)
  fy <- function(z) dautoreg(z, phi = phi_year, log = TRUE)
  var <- sd ^ 2 / ((1 - phi_age ^ 2) * (1 - phi_year ^ 2))
  dseparable(fy, fa)(x, scale = sqrt(var))

}

#' @rdname process_2d
#' @export
rprocess_2d <- function(ny, na, phi = c(0, 0), sd = 1) {

  phi_age  <- phi[1]
  phi_year <- phi[2]

  ar1_cor <- function(n, phi) stats::toeplitz(phi ^ (0:(n - 1L)))
  C_age  <- ar1_cor(na, phi_age)    # columns
  C_year <- ar1_cor(ny, phi_year)   # rows

  var <- sd^2 / ((1 - phi_age^2) * (1 - phi_year^2))
  Sigma <- var * kronecker(C_age, C_year)

  z <- MASS::mvrnorm(1L, mu = rep(0, ny * na), Sigma = Sigma)
  return(matrix(z, ny, na))

}


## Simple smooth switch for a ratio that starts at 1 and drops to 0 when r > 1
.logistic_switch <- function(ratio, slope = 10) {
  plogis(-slope * (ratio - 1))
}


#' Solve an F multiplier from a target total catch
#'
#' @description
#' Finds the scalar multiplier `k` such that scaling a terminal fishing-mortality
#' vector `F` by `k` matches a target catch `C_target` (numbers) under the
#' Baranov catch equation with given natural mortality `M` and abundance
#' `N`. Uses a fixed number of Newton updates starting from `k_init`.
#'
#' The total catch as a function of `k` is
#' \deqn{
#'   C(k) \;=\; \sum_a N_a \;
#'   \Big[1 - \exp\!\big(-Z_a(k)\big)\Big] \;
#'   \frac{k F_a}{Z_a(k)}, \qquad Z_a(k)=kF_a+M_a.
#' }
#'
#' The Newton update uses the analytic derivative
#' \deqn{
#'   \frac{dC}{dk}
#'   \;=\;
#'   \sum_a N_a \left[
#'     e^{-Z_a(k)} F_a \frac{kF_a}{Z_a(k)}
#'     \;+\;
#'     \big(1-e^{-Z_a(k)}\big) \frac{F_a M_a}{Z_a(k)^2}
#'   \right].
#' }
#'
#' @param C_target Numeric scalar; target total catch (same units as `N`).
#' @param F Numeric vector of terminal fishing mortality at age (\eqn{\ge} 0).
#' @param M Numeric vector of natural mortality at age (\eqn{\ge} 0).
#' @param N Numeric vector of abundance at age (\eqn{\ge} 0).
#' @param k_init Numeric scalar; initial guess for the multiplier (default `1`).
#' @param n_iter Integer; number of Newton updates (default `7`).
#'
#' @return Numeric scalar multiplier `k`.
#'
#' @details
#' This solves on the natural scale for `k`. If your target is in biomass,
#' pass biomass-at-age in `N` (i.e., `N <- numbers * weight_at_age`)
#' so the identity still holds.
#'
#' @examples
#' F <- rep(0.3, 5); M <- rep(0.2, 5); N <- 1000 * exp(-0.3 * (0:4))
#' C_of_k <- function(k) {
#'   Z <- k*F + M; eZ <- exp(-Z)
#'   sum(N * (1 - eZ) * (k*F / Z))
#' }
#' C_target <- C_of_k(1)   # target matching k = 1
#' solve_F_mult(C_target, F, M, N, k_init = 1)  # ~ 1
#'
#' # Higher target -> k > 1
#' solve_F_mult(1.5 * C_target, F, M, N, k_init = 1)
#'
#' @export
solve_F_mult <- function(C_target, F, M, N,
                         k_init = 1, n_iter = 7) {
  stopifnot(is.numeric(C_target), length(C_target) == 1L,
            is.numeric(k_init), length(k_init) == 1L,
            is.numeric(F), is.numeric(M), is.numeric(N),
            length(F) == length(M),
            length(M)  == length(N))

  # Don't allow more catch than N)
  N_total <- sum(N)
  C_ratio <- C_target / N_total
  C_switch <- .logistic_switch(C_ratio)
  C_max <- C_switch * C_target

  k <- k_init
  for (i in seq_len(n_iter)) {
    Z   <- k * F + M
    eZ  <- exp(-Z)
    Ck  <- sum(N * (1 - eZ) * (k * F / Z))
    dCk <- sum(N * (eZ * F * (k * F / Z) +
                             (1 - eZ) * (F * M / (Z * Z))))
    step <- (C_max - Ck) / dCk
    k <- k + step
  }

  k + 0.001
}





#' Negative log-likelihood (and simulator) for the Tiny Assessment Model
#'
#' @description
#' Core objective function for TAM.
#'
#' - When `simulate = FALSE` (default) it returns the joint negative log-likelihood (JNLL)
#'   of the state–space model given parameters in `par` and data/flags in the captured `dat` list.
#' - When `simulate = TRUE`, it draws the model’s random effects and observations
#'   from the assumed distributions and returns a list of simulated objects
#'   (see **Value**).
#'
#' @details
#' The model follows a standard age–structured state–space formulation:
#'
#' - **Recruitment:** log-recruits \eqn{\log R_y} evolve as a random walk:
#'   \deqn{\Delta \log R_y \sim \mathcal{N}(0,\sigma_R^2).}
#'
#'   If `N_settings$init_N0` is `TRUE`, the initial recruit level
#'   is constrained by:
#'   \deqn{\log R_1 \sim \mathcal{N}(\log R_0,\sigma_R^2).}
#'
#' - **Numbers-at-age:** forward cohort dynamics with plus-group:
#'   \deqn{\log N_{y,a} = \log N_{y-1,a-1} - Z_{y-1,a-1},}
#'
#'   with \eqn{Z_{y,a} = F_{y,a} + M_{y,a}}. The plus-group equation is applied
#'   at the terminal age.
#'   If `N_settings$process != "off"`, residuals
#'   \eqn{\eta^N_{y,a} = \log N_{y,a} - \widehat{\log N}_{y,a}}
#'   are penalized by [dprocess_2d()] according to the chosen process.
#'
#' - **Fishing mortality:**
#'   \deqn{\log F_{y,a} = \mu^F_{y,a} + \eta^F_{y,a},}
#'
#'   where the optional fixed-effects surface \eqn{\mu^F} comes from
#'   \eqn{F_\text{modmat} \cdot \texttt{log\_mu\_f}} (if `F_settings$mu_form` is provided).
#'   Deviations \eqn{\eta^F} are penalized by [dprocess_2d()] using
#'   `F_settings$process` and `logit_phi_f` (AR1) or a RW/IID penalty.
#'
#' - **Natural mortality:**
#'   \deqn{\log M_{y,a} = \log \mu^M_{y,a} + \eta^M_{y,a},}
#'
#'   where \eqn{\log \mu^M = \texttt{log\_mu\_assumed\_m} + M_\text{modmat}\,\texttt{log\_mu\_m}}.
#'   If `M_settings$process != "off"`, deviations \eqn{\eta^M} are penalized
#'   by [dprocess_2d()] and may be grouped by age via `M_settings$age_blocks`.
#'
#' - **Observations:** catch-at-age and index-at-age on the log scale:
#'   \deqn{\log C_{y,a} \sim \mathcal{N}\!\left(
#'       \log\!\left[N_{y,a}\,\frac{F_{y,a}}{Z_{y,a}}\,(1-e^{-Z_{y,a}})\right],
#'       \sigma^2_{\text{obs}}\right),}
#'
#'   \deqn{\log I_{y,a} \sim \mathcal{N}\!\left(
#'       \log q_{a} + \log N_{y,a} - Z_{y,a}\, t_{y,a}, \sigma^2_{\text{obs}}\right).}
#'
#'   Here `sd_obs_modmat %*% log_sd_obs` controls observation SD (by block),
#'   and `q_modmat %*% log_q` controls age- (or block-) specific catchability.
#'
#' **Simulation mode:**
#' When `simulate = TRUE`, the function:
#'
#' 1. Draws `log_r` (RW), optional `log_n` (cohort residual field),
#'    optional `log_m` (M deviations), and `log_f` (F deviations) from
#'    their respective process models via [rprocess_2d()].
#' 2. Regenerates predictions and draws `log_obs` from the observation
#'    model.
#' 3. Returns the simulated objects.
#'
#' Missing observations are preserved (filled and then reset to `NA`).
#'
#' `REPORT()` and `ADREPORT()` calls inside the function make derived
#' quantities (e.g., `N`, `F`, `M`, `Z`, `ssb`, `log_ssb`) available through
#' `obj$report()` / `sdreport()` when used via **RTMB**.
#'
#' @param par Named list of parameters in the format produced by
#'   [make_par()]. This includes scalars (e.g., `log_sd_*`), vectors
#'   (e.g., `log_r`, `log_q`), and matrices (e.g., `log_f`, `log_n`, `log_m`).
#' @param dat Named list of data and setting inputs produced by [make_dat()].
#' @param simulate Logical. If `FALSE`, return the JNLL.
#'   If `TRUE`, simulate random effects and observations and return them (see **Value**).
#'
#' @return
#' - If `simulate = FALSE`: a single numeric JNLL value.
#' - If `simulate = TRUE`: a list with elements:
#'   - `log_f`, `log_r` — always returned;
#'   - `log_n` — if `N_settings$process != "off"`;
#'   - `log_m` — if `M_settings$process != "off"`;
#'   - `log_obs` — simulated observations (NAs restored where input was missing);
#'   - `missing` — the simulated values at missing-observation positions.
#'
#' @section Dependencies and captured data:
#' The function expects a `dat` list in its lexical scope (created by
#' [make_dat()]) containing data matrices/vectors and model matrices
#' (`SW`, `MO`, `obs_map`, `sd_obs_modmat`, `q_modmat`, `F_modmat`,
#' `M_modmat`, settings lists, etc.).
#' It also relies on helper functions [dprocess_2d()] and [rprocess_2d()]
#' for process penalties and simulation.
#'
#' @examples
#' dat <- make_dat(
#'   cod_obs,
#'   years = 1983:2024,
#'   ages  = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' par <- make_par(dat)
#' make_nll_fun <- function(f, d) function(p) f(p, d)
#' obj <- RTMB::MakeADFun(make_nll_fun(nll_fun, dat), par,
#'   random = c("log_n", "log_f","log_r", "missing"), silent = TRUE
#' )
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' rep <- obj$report()
#' sdrep <- RTMB::sdreport(obj)
#'
#' # Simulate from fitted parameters
#' p_hat <- as.list(sdrep, "Estimate")
#' sims  <- nll_fun(p_hat, dat, simulate = TRUE)
#'
#' @seealso
#' [make_dat()], [make_par()], [fit_tam()], [sim_tam()],
#' [dprocess_2d()], [rprocess_2d()]
#' @export
nll_fun <- function(par, dat, simulate = FALSE) {

  "[<-" <- ADoverload("[<-")

  getAll(par, dat)

  log_obs <- OBS(log_obs)
  log_obs[is_na_obs] <- missing

  n_obs <- length(log_obs)
  n_years <- length(years)
  n_ages <- length(ages)
  n_proj <- proj_settings$n_proj

  sd_r <- exp(log_sd_r)
  sd_f <- exp(log_sd_f)

  empty_mat <- matrix(NA, n_years, n_ages,
                      dimnames = list(year = years, age = ages))
  log_F <- log_mu_F <- S <- empty_mat
  N <- log_N <- pred_log_N <- empty_mat
  M <- log_mu_M  <- empty_mat
  Z <- empty_mat

  ## Vital rates ---

  log_R <- log_r
  log_N[, 1] <- log_r

  log_F[!is_proj, ] <- log_f
  log_mu_F[] <- drop(F_modmat %*% log_mu_f)
  mu_F <- exp(log_mu_F)
  F <- exp(log_F)

  log_mu_M[] <- log_mu_assumed_m + drop(M_modmat %*% log_mu_m)
  M <- mu_M <- exp(log_mu_M)
  if (M_settings$process != "off") {
    M[-1, -1] <- exp(log_mu_M[-1, -1] + log_m[, M_settings$age_blocks])
  }
  log_M <- log(M)
  Z <- F + M
  log_Z <- log(Z)


  ## Cohort equation (assumes max age = plus group) ---

  .pred_cohorts <- function(log_N, Z, fill = FALSE) {
    ny <- nrow(log_N); na <- ncol(log_N)
    iy  <- 2:ny; ia <- 2:na
    if (fill) {
      for (a in ia) log_N[iy, a] <- log_N[iy - 1, a - 1] - Z[iy - 1, a - 1]
      log_N[iy, na] <- log(exp(log_N[iy, na]) + exp(log_N[iy - 1, na] - Z[iy - 1, na]))
      return(log_N)
    } else {
      pred <- matrix(NA, nrow = ny, ncol = na, dimnames = dimnames(log_N))
      pred[iy, ia] <- log_N[iy - 1, ia - 1] - Z[iy - 1, ia - 1]
      pred[iy, na] <- log(exp(pred[iy, na]) + exp(log_N[iy - 1, na] - Z[iy - 1, na]))
      return(pred)
    }
  }

  if (N_settings$init_N0) {
    log_N[1, 2:n_ages] <- log_r0 - cumsum(Z[1, 2:n_ages - 1])
  }
  if (N_settings$process == "off") {
    log_N[!is_proj, ] <- .pred_cohorts(log_N[!is_proj, ], Z[!is_proj, ], fill = TRUE)
  } else {
    log_N[, -1] <- log_n
    pred_log_N[!is_proj, ] <- .pred_cohorts(log_N[!is_proj, ], Z[!is_proj, ])
  }
  N <- exp(log_N)


  ## Projections ---

  browser()

  # Replicate terminal selectivity across proj_years
  F_last <- F[sum(!is_proj), ]
  if (proj_settings$n_proj > 0) {
    for (y in as.character(proj_years)) {
      prev_y <- as.character(as.numeric(y) - 1)
      yy <- c(prev_y, y)
      if (N_settings$process == "off") {
        log_N[y, ] <- .pred_cohorts(log_N[yy, ], Z[yy, ], fill = TRUE)[y, ]
      } else {
        pred_log_N[y, ] <- .pred_cohorts(log_N[yy, ], Z[yy, ])[y, ]
      }
      k <- solve_F_mult(C_target = proj_settings$tac[y],
                        F = F_last,
                        M = M[y, ],
                        N = N[y, ])
      F[y, ] <- k * F_last
      Z[y, ] <- F[y, ] + M[y, ]
      log_F[y, ] <- log(F[y, ])
      log_Z[y, ] <- log(Z[y, ])
      N[y, ] <- exp(log_N[y, ])
    }
  }


  ## Transformations etc. ---

  F_full <- apply(F, 1, max)
  S <- sweep(F, 1, F_full, "/")

  ssb_mat <- SW * MO * N * exp(-Z)
  ssb <- rowSums(ssb_mat)
  log_ssb_mat <- log(ssb_mat)
  log_ssb <- log(ssb)


  ## Recruitment deviations (basic random walk) ---

  jnll <- 0

  if (N_settings$init_N0) {
    jnll <- jnll - dnorm(log_N[1, 1], mean = log_r0, sd = sd_r, log = TRUE)
    if (simulate) {
      log_r[1] <- rnorm(1, mean = log_r0, sd = sd_r)
    }
  }
  eta_R <- log_N[2:n_years, 1] - log_N[1:(n_years - 1), 1]
  jnll <- jnll - sum(dnorm(eta_R, 0, sd_r, log = TRUE))
  if (simulate) {
    eta_R <- rnorm(n_years - 1, 0, sd = sd_r)
    log_r[2:n_years] <- log_r[1:(n_years - 1)] + eta_R
  }


  ## Cohort deviations ---

  if (N_settings$process != "off") {
    eta_log_N <- log_N[-1, -1] - pred_log_N[-1, -1]
    sd_n <- exp(log_sd_n)
    phi <- plogis(logit_phi_n)
    jnll <- jnll - dprocess_2d(eta_log_N, sd = sd_n, phi = phi)
    if (simulate) {
      eta_log_N <- rprocess_2d(n_years - 1, n_ages - 1, sd = sd_n, phi = phi)
      log_n[-1, ] <- pred_log_N[-1, -1] + eta_log_N
    }
  }

  ## M deviations ---

  if (M_settings$process != "off") {
    eta_log_m <- log_m - log_mu_M[-1, -1]
    sd_m <- exp(log_sd_m)
    phi <- plogis(logit_phi_m)
    jnll <- jnll - dprocess_2d(eta_log_m, sd = sd_m, phi = phi)
    if (simulate) {
      eta_log_m <- rprocess_2d(nrow(log_m), ncol(log_m), sd = sd_m, phi = phi)
      log_m <- log_mu_M + eta_log_m
    }
  }

  ## F deviations ---

  eta_log_F <- log_F[!is_proj, ] - log_mu_F[!is_proj, ]
  phi <- plogis(logit_phi_f)
  jnll <- jnll - dprocess_2d(eta_log_F, sd = sd_f, phi = phi)
  if (simulate) {
    eta_log_F <- rprocess_2d(nrow(log_f), ncol(log_f), sd = sd_f, phi = phi)
    log_f <- log_mu_F + eta_log_F
  }

  ## Observations ---

  log_pred <- numeric(n_obs)
  iya <- sapply(obs_map[, c("year", "age")], as.character)
  N_obs <- N[iya]
  Z_obs <- Z[iya]
  F_obs <- F[iya]
  sd_obs <- exp(drop(sd_obs_modmat %*% log_sd_obs))
  log_q_obs <- drop(q_modmat %*% log_q)    # length = number of survey index rows
  samp_time <- obs_map$samp_time

  idx <- obs_map$type == "catch"
  log_pred[idx] <- log(N_obs[idx]) - log(Z_obs[idx]) + log(1 - exp(- Z_obs[idx])) + log(F_obs[idx])

  idx <- obs_map$type == "index"
  log_pred[idx] <- log_q_obs + log(N_obs[idx]) - Z_obs[idx] * samp_time[idx]

  jnll <- jnll - sum(dnorm(log_obs, log_pred, sd = sd_obs, log = TRUE))
  if (simulate) {
    log_obs <- rnorm(n_obs, mean = log_pred, sd = sd_obs)
  }

  REPORT(N)
  REPORT(M)
  REPORT(mu_M)
  REPORT(F)
  REPORT(mu_F)
  REPORT(F_full)
  REPORT(S)
  REPORT(Z)
  REPORT(ssb_mat)
  REPORT(ssb)
  REPORT(log_pred)
  REPORT(log_obs)
  REPORT(sd_obs)
  REPORT(log_q_obs)

  ADREPORT(log_ssb)

  if (simulate) {
    sims <- list(log_f = log_f,
                 log_r = log_r,
                 log_obs = log_obs)
    if (length(missing) > 0) {
      sims$missing <- log_obs[is_na_obs]
      sims$log_obs[is_na_obs] <- NA
    }
    if (N_settings$process != "off") {
      sims$log_n <- log_n
    }
    if (M_settings$process != "off") {
      sims$log_m <- log_m
    }
    return(sims)
  }

  jnll

}
