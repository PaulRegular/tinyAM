
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

  ## Vital rates ----

  log_recruitment <- log_r
  log_N[, 1] <- log_r

  log_F[!is_proj, ] <- log_f
  if (n_proj > 0) {
    log_k <- log(proj_settings$F_mult)
    log_f_last <- log_f[rep(nrow(log_f), n_proj), , drop = FALSE]
    proj_log_F <- sweep(log_f_last, 1, log_k, `+`)
    log_F[is_proj, ] <- proj_log_F
  }
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


  ## Cohort equation (assumes max age = plus group) ----

  Y <- 2:n_years
  A <- 2:n_ages
  if (N_settings$init_N0) {
    log_N[1, A] <- log_r0 - cumsum(Z[1, A - 1])
  }
  if (N_settings$process == "off") {
    for (a in A) {
      log_N[Y, a] <- log_N[Y - 1, a - 1] - Z[Y - 1, a - 1]
    }
    log_N[Y, n_ages] <- log(exp(log_N[Y, n_ages]) + exp(log_N[Y - 1, n_ages] - Z[Y - 1, n_ages]))
  } else {
    log_N[, -1] <- log_n
    pred_log_N[Y, A] <- log_N[Y - 1, A - 1] - Z[Y - 1, A - 1]
    pred_log_N[Y, n_ages] <- log(exp(pred_log_N[Y, n_ages]) + exp(log_N[Y - 1, n_ages] - Z[Y - 1, n_ages]))
  }
  N <- exp(log_N)


  ## Recruitment deviations (basic random walk) ----

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


  ## Cohort deviations ----

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

  ## M deviations ----

  if (M_settings$process != "off") {
    eta_log_M <- log_m[, M_settings$age_blocks] - log_mu_M[-1, -1]
    eta_log_m <- eta_log_M[, !duplicated(M_settings$age_blocks)]
    sd_m <- exp(log_sd_m)
    phi <- plogis(logit_phi_m)
    jnll <- jnll - dprocess_2d(eta_log_m, sd = sd_m, phi = phi)
    if (simulate) {
      eta_log_m <- rprocess_2d(nrow(log_m), ncol(log_m), sd = sd_m, phi = phi)
      log_m <- log_mu_M[-1, -1][, !duplicated(M_settings$age_blocks)] + eta_log_m
    }
  }

  ## F deviations ----

  eta_log_F <- log_F[!is_proj, ] - log_mu_F[!is_proj, ]
  phi <- plogis(logit_phi_f)
  jnll <- jnll - dprocess_2d(eta_log_F, sd = sd_f, phi = phi)
  if (simulate) {
    eta_log_F <- rprocess_2d(nrow(log_f), ncol(log_f), sd = sd_f, phi = phi)
    log_f <- log_mu_F + eta_log_F
  }


  ## Observations ----

  log_pred <- numeric(n_obs)
  iya <- sapply(obs_map[, c("year", "age")], as.character)
  N_obs <- N[iya]
  Z_obs <- Z[iya]
  F_obs <- F[iya]
  sd_obs <- exp(drop(sd_obs_modmat %*% log_sd_obs))
  log_q_obs <- drop(q_modmat %*% log_q)    # length = number of survey index rows
  samp_time <- obs_map$samp_time

  ic <- obs_map$type == "catch"
  log_pred[ic] <- log(N_obs[ic]) - log(Z_obs[ic]) + log(1 - exp(- Z_obs[ic])) + log(F_obs[ic])

  ii <- obs_map$type == "index"
  log_pred[ii] <- log_q_obs + log(N_obs[ii]) - Z_obs[ii] * samp_time[ii]

  jnll <- jnll - sum(dnorm(log_obs, log_pred, sd = sd_obs, log = TRUE))
  if (simulate) {
    log_obs <- rnorm(n_obs, mean = log_pred, sd = sd_obs)
  }

  ## Derived quantities ----

  F_full <- apply(F, 1, max)
  S <- sweep(F, 1, F_full, "/")

  ia <- as.character(F_settings$mean_ages)
  F_bar <- rowSums(F[, ia] * N[, ia]) / rowSums(N[, ia])
  log_F_bar <- log(F_bar)
  ia <- as.character(M_settings$mean_ages)
  M_bar <- rowSums(M[, ia] * N[, ia]) / rowSums(N[, ia])
  log_M_bar <- log(M_bar)

  abundance <- rowSums(N)
  log_abundance <- log(abundance)
  biomass_mat <- W * N
  biomass <- rowSums(biomass_mat)
  log_biomass <- log(biomass)
  ssb_mat <- W * P * N
  ssb <- rowSums(ssb_mat)
  log_ssb <- log(ssb)

  C_obs <- C_pred <- empty_mat
  C_obs[] <- exp(log_obs[ic])
  C_pred[] <- exp(log_pred[ic])
  total_catch <- rowSums(C_obs)
  total_catch_pred <- rowSums(C_pred)
  total_yield <- rowSums(C_obs * W)
  total_yield_pred <- rowSums(C_pred * W)


  ## Output ----

  REPORT(N)
  REPORT(abundance)
  REPORT(M)
  REPORT(mu_M)
  REPORT(M_bar)
  REPORT(F)
  REPORT(mu_F)
  REPORT(F_full)
  REPORT(F_bar)
  REPORT(S)
  REPORT(Z)
  REPORT(biomass_mat)
  REPORT(biomass)
  REPORT(ssb_mat)
  REPORT(ssb)

  REPORT(total_catch)
  REPORT(total_catch_pred)
  REPORT(total_yield)
  REPORT(total_yield_pred)

  REPORT(log_pred)
  REPORT(log_obs)
  REPORT(sd_obs)
  REPORT(log_q_obs)

  ADREPORT(log_recruitment)
  ADREPORT(log_abundance)
  ADREPORT(log_biomass)
  ADREPORT(log_ssb)
  ADREPORT(log_F_bar)
  ADREPORT(log_M_bar)

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
