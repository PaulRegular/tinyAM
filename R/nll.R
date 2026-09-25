
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
#' **Latent-state convention:** `log_r`, `log_n0`, `log_n`, `log_f`, and `log_m`
#' represent latent quantities on the log scale, not process deviations.
#' Lowercase names denote compact fitted latent-state parameters; `log_N`,
#' `log_F`, and `log_M` denote full model surfaces after cohort recursion,
#' age-block expansion, and/or projection. Process errors are calculated
#' internally as deviations (`eta_*`) from expected or mean states:
#' `eta_R` uses successive log-recruitment states, `eta_log_N` uses cohort
#' predictions, and `eta_log_f` and `eta_log_m` use log mean surfaces.
#'
#' The model follows a standard age–structured state–space formulation:
#'
#' - **Recruitment:** log-recruits \eqn{\log R_y} evolve as a random walk:
#'   \deqn{\Delta \log R_y \sim \mathcal{N}(0,\sigma_R^2).}
#'
#'   The fixed parameter `log_r0` is the actual first-year state, not a
#'   hypermean. Random states `log_r` contain years 2:Y only. The full path is
#'   `log_recruitment = c(log_r0, log_r)`; the first innovation is
#'   `log_r[1] - log_r0`, followed by successive differences.
#'
#' - **Initial older-age abundance:** independent of `N_settings$process`.
#'   The recursion starts from `log_N[1, 1] = log_r0`. For every older age,
#'   the prediction is `log_N[1, a-1] - Z[1, a-1]`, including the terminal age
#'   without an equilibrium plus-group adjustment. `init = "exp"` uses these
#'   predictions directly. Under `"free"` or `"random"`, `log_n0` contains the
#'   realized initial log abundance at all older ages. The residual is
#'   `eta_log_n0 = log_n0 - (c(log_r0, head(log_n0, -1)) - Z[1, -n_ages])`.
#'   `"free"` estimates these states without a penalty; `"random"` applies an
#'   IID normal density to the residuals with SD `exp(log_sd_n0)`. This SD
#'   describes the initial age margin, separately from temporal recruitment
#'   SD `sd_r` and subsequent cohort-process SD `sd_n`. See [make_dat()].
#'
#' - **Numbers-at-age:** forward cohort dynamics with plus-group:
#'   \deqn{\log N_{y,a} = \log N_{y-1,a-1} - Z_{y-1,a-1},}
#'
#'   with \eqn{Z_{y,a} = F_{y,a} + M_{y,a}}. The plus-group equation is applied
#'   at the terminal age for transitions from year 1 to year 2 onward.
#'   Latent older-age states `log_n` contain years 2:Y only.
#'   If `N_settings$process != "off"`, residuals
#'   \eqn{\eta^N_{y,a} = \log N_{y,a} - \widehat{\log N}_{y,a}}
#'   are penalized by [dprocess_2d()] according to the chosen process.
#'
#' - **Fishing mortality:**
#'   \deqn{\log F_{y,a} = \log \mu^F_{y,a} + \eta^F_{y,a},}
#'
#'   where the log mean surface \eqn{\log \mu^F} comes from
#'   \eqn{F_\text{modmat} \cdot \texttt{log\_mu\_f}} if `F_settings$mu_form`
#'   is provided, and is zero otherwise. The latent state `log_f` represents
#'   realized absolute log fishing mortality in observed years. Its process
#'   deviation is `eta_log_f = log_f - log_mu_F`, with `log_mu_F` restricted
#'   to those years. These deviations are penalized by [dprocess_2d()] using
#'   `F_settings$process` and `logit_phi_f` (AR1) or an approximate RW/IID penalty.
#'
#' - **Natural mortality:**
#'   \deqn{\log M_{y,a} = \log \mu^M_{y,a} + \eta^M_{y,a},}
#'
#'   where \eqn{\log \mu^M = \texttt{log\_mu\_supplied\_m} + M_\text{modmat}\,\texttt{mu\_m}}.
#'   When `M_settings$process != "off"`, the latent state `log_m` represents
#'   realized absolute log natural mortality from `M_settings$first_dev_year`
#'   onward. Its process deviation is `eta_log_m = log_m - log_mu_M`, with
#'   `log_mu_M` restricted to those years and the age-block starts. These
#'   deviations are penalized by [dprocess_2d()] using `M_settings$process`
#'   and `logit_phi_m` (AR1) or an approximate RW/IID penalty.
#'
#' - **Observations:** catch-at-age and index-at-age on the log scale:
#'   \deqn{\log C_{y,a} \sim \mathcal{N}\!\left(
#'       \log\!\left[N_{y,a}\,\frac{F_{y,a}}{Z_{y,a}}\,(1-e^{-Z_{y,a}})\right],
#'       \sigma^2_{\text{catch}}\right),}
#'
#'   \deqn{\log I_{y,a} \sim \mathcal{N}\!\left(
#'       \log q_{a} + \log N_{y,a} - Z_{y,a}\, t_{y,a}, \sigma^2_{\text{index}}\right).}
#'
#'   Here `sd_catch_modmat %*% log_sd_catch` adjusts the supplied observation SDs for catch-at-age,
#'   `sd_index_modmat %*% log_sd_index` does the same for indices-at-age, and `q_modmat %*% log_q`
#'   controls age- (or block-) specific catchability.
#'   With [mono()] terms, `q_mono_modmat %*% dq` is added to this
#'   ordinary component. Non-negative increments enforce non-decreasing q,
#'   with independent steps per `by` group. The survey observation
#'   equation and its SD are otherwise unchanged.
#'
#' **Simulation mode:**
#' When `simulate = TRUE`, the function:
#'
#' 1. Generates latent states: `log_r` using recruitment RW increments from
#'    the fixed first-year anchor `log_r0`,
#'    `log_f` and optional `log_m` by adding process deviations to their log
#'    mean surfaces, and optional `log_n` by adding process deviations to
#'    recursive cohort predictions. Process fields are drawn via [rprocess_2d()];
#'    recruitment increments use [stats::rnorm()]. Initial states without a
#'    specified process distribution retain their supplied values. Random
#'    initial-age residuals are drawn after F/M and Z are constructed, then
#'    `log_n0` is built recursively from `log_r0`. Fixed `log_n0` states under
#'    free initialization are retained.
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
#'   - `log_n0` — if `N_settings$init` is `"free"` or `"random"`;
#'   - `log_n` — if `N_settings$process != "off"`;
#'   - `log_m` — if `M_settings$process != "off"`;
#'   - `log_obs` — simulated observations (NAs restored where input was missing);
#'   - `missing` — the simulated values at missing-observation positions.
#'
#' @section Dependencies and captured data:
#' The function expects a `dat` list in its lexical scope (created by
#' [make_dat()]) containing data matrices/vectors and model matrices
#' (`SW`, `MO`, `obs_map`, `sd_catch_modmat`, `sd_index_modmat`, `q_modmat`,
#' `F_modmat`, `M_modmat`, settings lists, etc.).
#' It also relies on helper functions [dprocess_2d()] and [rprocess_2d()]
#' for process penalties and simulation.
#'
#' @example inst/examples/example_dat_default.R
#' @examples
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
#' @importFrom stats rnorm
#'
#' @seealso
#' [make_dat()], [make_par()], [fit_tam()], [sim_tam()],
#' [dprocess_2d()], [rprocess_2d()]
#' @export
nll_fun <- function(par, dat, simulate = FALSE) {

  "[<-" <- ADoverload("[<-")

  getAll(par, dat)

  observed <- OBS(observed)
  if (any(fill_missing_map)) {
    log_obs[fill_missing_map] <- missing
  }

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

  ## Mean structures and process draws ----

  log_mu_F[] <- drop(F_modmat %*% log_mu_f)
  log_mu_M[] <- log_mu_supplied_m + drop(M_modmat %*% mu_m)
  if (simulate) {
    log_r[] <- log_r0 + cumsum(stats::rnorm(n_years - 1, 0, sd_r))
    log_f[] <- log_mu_F[!is_proj, ] +
      rprocess_2d(nrow(log_f), ncol(log_f), sd = sd_f, phi = plogis(logit_phi_f))
    if (M_settings$process != "off") {
      iy <- rownames(log_m)
      ia <- M_settings$age_block_start
      log_m[] <- log_mu_M[iy, ia, drop = FALSE] +
        rprocess_2d(nrow(log_m), ncol(log_m), sd = exp(log_sd_m), phi = plogis(logit_phi_m))
    }
  }

  ## Vital rates ----

  log_recruitment <- c(log_r0, log_r)
  names(log_recruitment) <- years
  recruitment <- exp(log_recruitment)
  log_N[, 1] <- log_recruitment

  log_F[!is_proj, ] <- log_f
  if (n_proj > 0) {
    log_k <- log(proj_settings$F_mult)
    log_f_last <- log_f[rep(nrow(log_f), n_proj), , drop = FALSE]
    proj_log_F <- sweep(log_f_last, 1, log_k, `+`)
    log_F[is_proj, ] <- proj_log_F
  }
  mu_F <- exp(log_mu_F)
  F <- exp(log_F)

  M <- mu_M <- exp(log_mu_M)
  if (M_settings$process != "off") {
    iy <- rownames(log_m)
    ia <- names(M_settings$age_blocks)
    M[iy, ia] <- exp(log_m[, M_settings$age_blocks, drop = FALSE])
  }
  log_M <- log(M)
  Z <- F + M
  log_Z <- log(Z)


  ## Initial abundance (independent of the subsequent N process) ----

  eta_log_n0 <- numeric(n_ages - 1L)
  if (simulate && N_settings$init == "random") {
    eta_log_n0[] <- stats::rnorm(n_ages - 1L, 0, exp(log_sd_n0))
  }
  for (a in 2:n_ages) {
    pred_log_N[1, a] <- log_N[1, a - 1] - Z[1, a - 1]
    if (N_settings$init == "exp" || (simulate && N_settings$init == "random")) {
      log_N[1, a] <- pred_log_N[1, a] + eta_log_n0[a - 1L]
    } else {
      log_N[1, a] <- log_n0[a - 1L]
    }
  }
  eta_log_n0 <- log_N[1, -1] - pred_log_N[1, -1]
  if (simulate && N_settings$init == "random") {
    log_n0[] <- log_N[1, -1]
  }

  ## Cohort equation (plus group after the initial year) ----

  Y <- 2:n_years
  A <- 2:n_ages
  if (N_settings$process != "off") {
    log_N[-1, -1] <- log_n
  }
  eta_log_N <- matrix(0, n_years - 1, n_ages - 1)
  if (simulate && N_settings$process != "off") {
    eta_log_N <- rprocess_2d(n_years - 1, n_ages - 1,
                            sd = exp(log_sd_n), phi = plogis(logit_phi_n))
  }
  for (y in Y) {
    pred_log_N[y, A] <- log_N[y - 1, A - 1] - Z[y - 1, A - 1]
    pred_log_N[y, n_ages] <- RTMB::logspace_add(pred_log_N[y, n_ages],
                                              log_N[y - 1, n_ages] - Z[y - 1, n_ages])
    if (N_settings$process == "off" || simulate) {
      log_N[y, A] <- pred_log_N[y, A] + eta_log_N[y - 1, ]
    }
  }
  if (simulate && N_settings$process != "off") {
    log_n[] <- log_N[-1, -1, drop = FALSE]
  }
  N <- exp(log_N)


  ## Initial age process ----

  jnll <- 0

  if (N_settings$init == "random") {
    jnll <- jnll - sum(RTMB::dnorm(eta_log_n0, 0, exp(log_sd_n0), log = TRUE))
  }

  ## Recruitment process (basic random walk) ----

  eta_R <- log_N[2:n_years, 1] - log_N[1:(n_years - 1), 1]
  jnll <- jnll - sum(RTMB::dnorm(eta_R, 0, sd_r, log = TRUE))


  ## N process ----

  if (N_settings$process != "off") {
    eta_log_N <- log_N[-1, -1, drop = FALSE] - pred_log_N[-1, -1, drop = FALSE]
    sd_n <- exp(log_sd_n)
    phi <- plogis(logit_phi_n)
    jnll <- jnll - dprocess_2d(eta_log_N, sd = sd_n, phi = phi)
  }

  ## M process ----

  if (M_settings$process != "off") {
    iy <- rownames(log_m)
    ia  <- dat$M_settings$age_block_start
    eta_log_m <- log_m - log_mu_M[iy, ia, drop = FALSE]
    sd_m <- exp(log_sd_m)
    phi  <- plogis(logit_phi_m)
    jnll <- jnll - dprocess_2d(eta_log_m, sd = sd_m, phi = phi)
  }

  ## F process ----

  eta_log_f <- log_F[!is_proj, ] - log_mu_F[!is_proj, ]
  phi <- plogis(logit_phi_f)
  jnll <- jnll - dprocess_2d(eta_log_f, sd = sd_f, phi = phi)


  ## Observations ----

  log_pred <- numeric(n_obs)
  iya <- sapply(obs_map[, c("year", "age")], as.character)
  log_N_obs <- log_N[iya]
  Z_obs <- Z[iya]
  F_obs <- F[iya]
  if (ncol(sd_catch_modmat) > 0) {
    log_sd_catch_eff <- drop(sd_catch_modmat %*% log_sd_catch)
  } else {
    log_sd_catch_eff <- rep(0, nrow(sd_catch_modmat))
  }
  if (ncol(sd_index_modmat) > 0) {
    log_sd_index_eff <- drop(sd_index_modmat %*% log_sd_index)
  } else {
    log_sd_index_eff <- rep(0, nrow(sd_index_modmat))
  }
  sd_catch <- exp(log_sd_catch_supplied + log_sd_catch_eff)
  sd_index <- exp(log_sd_index_supplied + log_sd_index_eff)
  sd_obs <- c(sd_catch, sd_index)
  log_q_obs <- drop(q_modmat %*% log_q) # length = number of survey index rows
  if (!is.null(dat$q_mono_modmat)) {
    log_q_obs <- log_q_obs + drop(dat$q_mono_modmat %*% dq)
  }
  samp_time <- obs_map$samp_time

  ic <- obs_map$type == "catch"
  log_pred[ic] <- log_N_obs[ic] - log(Z_obs[ic]) + log(1 - exp(- Z_obs[ic])) + log(F_obs[ic])

  ii <- obs_map$type == "index"
  log_pred[ii] <- log_q_obs + log_N_obs[ii] - Z_obs[ii] * samp_time[ii]

  jnll <- jnll - sum(RTMB::dnorm(observed, log_pred[is_observed], sd = sd_obs[is_observed], log = TRUE))
  if (any(fill_missing_map)) {
    jnll <- jnll - sum(RTMB::dnorm(missing, log_pred[fill_missing_map], sd = sd_obs[fill_missing_map], log = TRUE))
  }
  if (simulate) {
    log_obs[is_observed] <- stats::rnorm(sum(is_observed), mean = log_pred[is_observed], sd = sd_obs[is_observed])
    if (any(fill_missing_map)) {
      log_obs[fill_missing_map] <- stats::rnorm(sum(fill_missing_map), mean = log_pred[fill_missing_map], sd = sd_obs[fill_missing_map])
    }
  }

  ## Derived quantities ----

  obs <- exp(log_obs)
  obs[is_missing] <- NA

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

  pred <- exp(log_pred)
  C_obs <- C_pred <- empty_mat
  C_obs[] <- obs[ic]
  C_pred[] <- pred[ic]
  total_catch <- rowSums(C_obs, na.rm = TRUE)
  total_catch_pred <- rowSums(C_pred, na.rm = TRUE)
  total_yield <- rowSums(C_obs * W, na.rm = TRUE)
  total_yield_pred <- rowSums(C_pred * W, na.rm = TRUE)


  ## Output ----

  REPORT(recruitment)
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
    if (N_settings$init != "exp") {
      sims$log_n0 <- log_n0
    }
    if (any(fill_missing_map)) {
      sims$missing <- log_obs[fill_missing_map]
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
