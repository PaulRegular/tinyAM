
#' Separable AR1 process: density and simulation
#'
#' @title Separable AR1 process error: density and simulation
#'
#' @description
#' Helpers for a two-dimensional (2D) AR1 process, with separable correlation
#' across years and ages. TAM uses these helpers only for `process = "ar1"`;
#' IID processes use Normal densities and draws directly.
#'
#' - [dprocess_ar1()] returns the log-density contribution for a given matrix `x`.
#' - [rprocess_ar1()] generates a matrix draw with the requested dependence structure.
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
#'   for [dprocess_ar1()].
#' @param ny,na Positive integers: numbers of years and ages for
#'   [rprocess_ar1()].
#' @param phi Length-2 numeric vector \code{c(phi_age, phi_year)} with
#'   values in \eqn{(-1, 1)} for AR(1); zero gives independence.
#' @param sd Positive innovation-scale SD \eqn{\sigma}; the marginal SD is
#'   \eqn{\sigma/\sqrt{(1-\phi_\text{age}^2)(1-\phi_\text{year}^2)}}.
#'
#' @return
#' - [dprocess_ar1()]: a single numeric log-density value.
#' - [rprocess_ar1()]: a numeric matrix of dimension \eqn{n_y \times n_a}.
#'
#' @examples
#' # Simulate then calculate log-density
#' set.seed(1)
#' X_ar1 <- rprocess_ar1(ny = 10, na = 8, sd = 0.3, phi = c(0.5, 0.9))
#' dprocess_ar1(X_ar1, sd = 0.3, phi = c(0.5, 0.9))
#'
#' @seealso
#' [RTMB::dseparable()], [RTMB::dautoreg()], [MASS::mvrnorm()]
#'
#' @import RTMB
#' @rdname process_ar1
#' @export
dprocess_ar1 <- function(x, phi = c(0, 0), sd = 1) {

  phi_age  <- phi[1]
  phi_year <- phi[2]
  fa <- function(z) dautoreg(z, phi = phi_age, log = TRUE)
  fy <- function(z) dautoreg(z, phi = phi_year, log = TRUE)
  var <- sd ^ 2 / ((1 - phi_age ^ 2) * (1 - phi_year ^ 2))
  dseparable(fy, fa)(x, scale = sqrt(var))

}

#' @rdname process_ar1
#' @export
rprocess_ar1 <- function(ny, na, phi = c(0, 0), sd = 1) {

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


#' Temporal random-walk process: density and simulation
#'
#' @description
#' A random walk allows process deviations to persist and change from one year
#' to the next, independently for each age or age block. There is no tendency
#' to return to a mean and no correlation between ages.
#'
#' @details
#' For a year by age matrix `x`, increments
#' `x[y, a] - x[y - 1, a]` are independent Normal variables with mean zero
#' and standard deviation `sd`. The first row has no increment density.
#' [rprocess_rw()] keeps that starting row and simulates subsequent rows
#' conditionally. In TAM these helpers act on process deviations from the
#' existing mean surfaces (F/M) or cohort predictions (N).
#'
#' @param x A numeric matrix with years in rows and ages or age blocks in columns.
#'   For simulation, only the first row supplies starting values.
#' @param sd Positive scalar standard deviation of year-to-year increments.
#' @return [dprocess_rw()] returns the summed log density of increments;
#'   [rprocess_rw()] returns a matrix with the dimensions and names of `x`.
#' @examples
#' x <- matrix(0, 5, 2)
#' set.seed(1)
#' x <- rprocess_rw(x, sd = 0.2)
#' dprocess_rw(x, sd = 0.2)
#' @rdname process_rw
#' @export
dprocess_rw <- function(x, sd = 1) {
  if (nrow(x) < 2L) return(0)
  increments <- x[-1, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  sum(RTMB::dnorm(increments, 0, sd, log = TRUE))
}

#' @rdname process_rw
#' @export
rprocess_rw <- function(x, sd = 1) {
  if (nrow(x) > 1L) {
    for (y in 2:nrow(x)) {
      x[y, ] <- x[y - 1L, ] + stats::rnorm(ncol(x), 0, sd)
    }
  }
  x
}



#' Evaluate or simulate the tinyAM model
#'
#' @description
#' Evaluates the population and observation equations for a parameter list.
#' Most users should use [fit_tam()] to fit a model and [sim_tam()] to generate
#' simulated data. This lower-level function is useful for checking model
#' equations or simulating from chosen parameter values.
#'
#' @details
#' See [tinyAM-model] for the complete state equations, process distributions,
#' likelihood, initial conditions, and derived quantities.
#'
#' **Latent-state convention:** `log_r`, `log_n0`, `log_n`, `log_f`, and `log_m`
#' represent absolute latent quantities on the log scale. Full model surfaces
#' are `log_N`, `log_F`, and `log_M`. Process deviations are separate quantities:
#' recruitment uses successive log states (`eta_R`), N uses cohort predictions
#' (`eta_log_N`), and F/M use their mean log surfaces (`eta_log_f`, `eta_log_m`).
#' In particular, `log_f = log_mu_F + eta_log_f` and
#' `log_m = log_mu_M + eta_log_m` on their represented years and ages/blocks.
#' The mean is not added again when constructing mortality from the latent state.
#'
#' With `simulate = FALSE`, the result is the joint negative log-likelihood.
#' It includes recruitment, active N/F/M processes, random initial abundance,
#' and Gaussian densities for observed or filled log observations. Random
#' effects are integrated out by [fit_tam()], not by this function itself.
#' Normal densities are evaluated on log observations, so predictions on the
#' natural scale are conditional medians rather than arithmetic means.
#' Catchability uses the link chosen by `index_settings$q_link`: the log link
#' exponentiates the formula predictor, while the logit link applies the inverse
#' logit and restricts q to between zero and one. Both estimation and simulation
#' use the same predictor, including any [mono()] increments.
#'
#' With `simulate = TRUE`, recruitment and mortality states are drawn first.
#' N is then constructed through initial-age and cohort recursion, followed by
#' predictions and observation draws with the SD for each matching row.
#' Derived quantities therefore use the same realization as the returned states.
#' If biological inputs extend above the modeled plus age, hidden age abundances
#' are reconstructed after N. Effective terminal W and P preserve biomass and
#' mature biomass; see [tinyAM-model]. No extra parameters or process penalties
#' are introduced. `N_plus`, `W`, and `P` are then included in `report()`.
#' Random N0 is redrawn; free N0 and `log_r0` remain supplied fixed states.
#' RW processes retain their starting state because it has no process density.
#' For N, the first `log_n` row is retained and its starting residual is computed
#' against the newly simulated cohort prediction. Subsequent residuals follow
#' the temporal walk. No new F states are drawn for projection years: projected
#' F is terminal historical F times the specified multiplier.
#'
#' @param par Parameter list with the structure produced by [make_par()].
#' @param dat Data and settings returned by [prepare_tam()].
#' @param simulate Logical; generate process states and observations instead of
#'   returning the likelihood? Defaults to `FALSE`.
#'
#' @return With `simulate = FALSE`, a scalar joint negative log-likelihood.
#'   When used inside [RTMB::MakeADFun()], reported population and observation
#'   quantities are available through the resulting object's `report()` method;
#'   selected log population summaries also receive uncertainty via `sdreport()`.
#'   With `TRUE`, a list containing `log_f`, `log_r`, `log_obs`, and applicable
#'   `log_n0`, `log_n`, `log_m`, and `missing` values. Unfilled missing
#'   observations remain `NA`; filled entries contain simulated log observations.
#'
#' @example inst/examples/example_dat_default.R
#' @examples
#' par <- make_par(dat)
#' set.seed(1)
#' simulated <- nll_fun(par, dat, simulate = TRUE)
#' head(exp(simulated$log_r))
#'
#' @importFrom stats rnorm
#' @seealso [tinyAM-model], [prepare_tam()], [make_par()], [fit_tam()], [sim_tam()]
#' @export
nll_fun <- function(f, d) function(p) f(p, d)
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
#' [prepare_tam()], [make_par()], [fit_tam()], [sim_tam()],
#' [dprocess_ar1()], [rprocess_ar1()]
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

  sd_f <- exp(log_sd_f)

  empty_mat <- matrix(NA, n_years, n_ages,
                      dimnames = list(year = years, age = ages))
  log_F <- log_mu_F <- S <- empty_mat
  N <- log_N <- pred_log_N <- empty_mat
  M <- log_mu_M  <- empty_mat
  Z <- empty_mat

  ## Mean structures and process draws ----

  F_effects <- .mean_effects(par, dat, "F", simulate)
  M_effects <- .mean_effects(par, dat, "M", simulate)
  log_mu_F[] <- drop(F_modmat %*% log_mu_f) + F_effects$contribution
  log_mu_M[] <- log_mu_supplied_m + drop(M_modmat %*% mu_m) + M_effects$contribution
  if (simulate) {
    mu_f <- log_mu_F[!is_proj, , drop = FALSE]
    log_f[] <- mu_f + if (F_settings$process == "rw") {
      rprocess_rw(log_f - mu_f, sd = sd_f)
    } else if (F_settings$process == "iid") {
      matrix(stats::rnorm(length(log_f), 0, sd_f), nrow(log_f), ncol(log_f))
    } else {
      rprocess_ar1(nrow(log_f), ncol(log_f), sd = sd_f, phi = plogis(logit_phi_f))
    }
    if (M_settings$process != "off") {
      iy <- rownames(log_m)
      ia <- M_settings$age_block_start
      mu_m_process <- log_mu_M[iy, ia, drop = FALSE]
      log_m[] <- mu_m_process + if (M_settings$process == "rw") {
        rprocess_rw(log_m - mu_m_process, sd = exp(log_sd_m))
      } else if (M_settings$process == "iid") {
        matrix(stats::rnorm(length(log_m), 0, exp(log_sd_m)), nrow(log_m), ncol(log_m))
      } else {
        rprocess_ar1(nrow(log_m), ncol(log_m), sd = exp(log_sd_m), phi = plogis(logit_phi_m))
      }
    }
  }

  ## Vital rates ----

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


  ## Population states and recruitment ----

  population <- .population_states(par, dat, Z, simulate)
  log_N <- population$log_N
  pred_log_N <- population$pred_log_N
  log_recruitment <- population$log_recruitment
  log_mu_R <- population$log_mu_R
  log_pred_R <- population$log_pred_R
  eta_R <- population$eta_R
  eta_R_state <- population$eta_R_state
  eta_log_n0 <- population$eta_log_n0
  W <- population$W
  P <- population$P
  N <- exp(log_N)
  recruitment <- exp(log_recruitment)
  if (simulate) {
    log_r[] <- log_recruitment[dat$rec$eligible]
    if (N_settings$init == "random") log_n0[] <- log_N[1, -1]
    if (N_settings$process != "off") log_n[] <- log_N[-1, -1, drop = FALSE]
  }
  if (!is.null(dat$plus_ages)) {
    N_plus <- population$N_plus
    REPORT(N_plus)
    REPORT(W)
    REPORT(P)
  }

  ## Initial age process ----

  jnll <- 0

  if (N_settings$init == "random") {
    jnll <- jnll - sum(RTMB::dnorm(eta_log_n0, 0, exp(log_sd_n0), log = TRUE))
  }

  ## Recruitment process ----

  jnll <- jnll + .rec_nll(log_recruitment, log_mu_R, par, dat)


  ## N process ----

  if (N_settings$process != "off") {
    eta_log_N <- log_N[-1, -1, drop = FALSE] - pred_log_N[-1, -1, drop = FALSE]
    sd_n <- exp(log_sd_n)
    jnll <- jnll - if (N_settings$process == "rw") {
      dprocess_rw(eta_log_N, sd = sd_n)
    } else if (N_settings$process == "iid") {
      sum(RTMB::dnorm(eta_log_N, 0, sd_n, log = TRUE))
    } else {
      dprocess_ar1(eta_log_N, sd = sd_n, phi = plogis(logit_phi_n))
    }
  }

  ## M process ----

  if (M_settings$process != "off") {
    iy <- rownames(log_m)
    ia  <- dat$M_settings$age_block_start
    eta_log_m <- log_m - log_mu_M[iy, ia, drop = FALSE]
    sd_m <- exp(log_sd_m)
    jnll <- jnll - if (M_settings$process == "rw") {
      dprocess_rw(eta_log_m, sd = sd_m)
    } else if (M_settings$process == "iid") {
      sum(RTMB::dnorm(eta_log_m, 0, sd_m, log = TRUE))
    } else {
      dprocess_ar1(eta_log_m, sd = sd_m, phi = plogis(logit_phi_m))
    }
  }

  ## F process ----

  eta_log_f <- log_F[!is_proj, ] - log_mu_F[!is_proj, ]
  jnll <- jnll - if (F_settings$process == "rw") {
    dprocess_rw(eta_log_f, sd = sd_f)
  } else if (F_settings$process == "iid") {
    sum(RTMB::dnorm(eta_log_f, 0, sd_f, log = TRUE))
  } else {
    dprocess_ar1(eta_log_f, sd = sd_f, phi = plogis(logit_phi_f))
  }


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
  sd_catch_effects <- .sd_effects(par, dat, "catch", simulate)
  sd_index_effects <- .sd_effects(par, dat, "index", simulate)
  sd_catch <- exp(log_sd_catch_supplied + log_sd_catch_eff + sd_catch_effects$contribution)
  sd_index <- exp(log_sd_index_supplied + log_sd_index_eff + sd_index_effects$contribution)
  sd_obs <- numeric(n_obs)
  sd_obs[obs_map$type == "catch"] <- sd_catch
  sd_obs[obs_map$type == "index"] <- sd_index
  q_coef <- if (identical(index_settings$q_link, "logit")) logit_q else log_q
  q_predictor <- drop(q_modmat %*% q_coef)
  if (!is.null(dat$q_mono_modmat)) {
    q_predictor <- q_predictor + drop(dat$q_mono_modmat %*% dq)
  }
  q_effects <- .q_effects(par, dat, simulate = simulate)
  jnll <- jnll + q_effects$nll + F_effects$nll + M_effects$nll
  jnll <- jnll + sd_catch_effects$nll + sd_index_effects$nll
  q_predictor <- q_predictor + q_effects$contribution
  # Compute log(q) directly to remain stable near the logit boundaries.
  log_q_obs <- if (identical(index_settings$q_link, "logit")) {
    -RTMB::logspace_add(0, -q_predictor)
  } else q_predictor
  log_q_obs <- log_q_obs + q_effects$log_selectivity
  samp_time <- obs_map$samp_time

  ic <- obs_map$type == "catch"
  log_pred[ic] <- log_N_obs[ic] - log(Z_obs[ic]) + log(-expm1(-Z_obs[ic])) + log(F_obs[ic])

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
  F_bar <- rowSums(F[, ia, drop = FALSE] * N[, ia, drop = FALSE]) / rowSums(N[, ia, drop = FALSE])
  log_F_bar <- log(F_bar)
  ia <- as.character(M_settings$mean_ages)
  M_bar <- rowSums(M[, ia, drop = FALSE] * N[, ia, drop = FALSE]) / rowSums(N[, ia, drop = FALSE])
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
  q_link_prediction <- if (identical(index_settings$q_link, "logit")) {
    if (any(vapply(dat$q_terms, function(term) term$type == "logistic", logical(1)))) {
      log_q_obs - log(-expm1(log_q_obs))
    } else q_predictor
  } else log_q_obs
  ADREPORT(q_link_prediction)
  eta_q_increments <- q_effects$rw_increments
  REPORT(eta_q_increments)
  if (length(eta_q_increments)) ADREPORT(eta_q_increments)
  eta_mu_F_increments <- F_effects$rw_increments
  eta_mu_M_increments <- M_effects$rw_increments
  REPORT(eta_mu_F_increments)
  REPORT(eta_mu_M_increments)
  if (length(eta_mu_F_increments)) ADREPORT(eta_mu_F_increments)
  if (length(eta_mu_M_increments)) ADREPORT(eta_mu_M_increments)
  eta_sd_catch_increments <- sd_catch_effects$rw_increments
  eta_sd_index_increments <- sd_index_effects$rw_increments
  REPORT(eta_sd_catch_increments)
  REPORT(eta_sd_index_increments)
  if (length(eta_sd_catch_increments)) ADREPORT(eta_sd_catch_increments)
  if (length(eta_sd_index_increments)) ADREPORT(eta_sd_index_increments)

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
    if (!is.null(par$log_r_init)) sims$log_r_init <- par$log_r_init
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
    return(c(sims, q_effects$parameters, F_effects$parameters, M_effects$parameters,
             sd_catch_effects$parameters, sd_index_effects$parameters))
  }

  jnll

}
