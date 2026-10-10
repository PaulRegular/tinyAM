
#' Initialize parameter list for TAM
#'
#' @description
#' `make_par()` builds a named list of initial values and shapes for all
#' fixed and random-effect parameters used by TAM, based on the structure in
#' a previously constructed `dat` list (see [prepare_tam()]).
#'
#' @details
#' **Latent-state convention:** `log_r`, `log_n0`, `log_n`, `log_f`, and `log_m`
#' represent latent quantities on the log scale. They are compact fitted
#' latent-state parameters; `log_N`, `log_F`, and `log_M` are full model
#' surfaces constructed internally by [nll_fun()]. Process errors are
#' calculated separately as deviations (`eta_*`): recruitment uses successive
#' log states, abundance uses cohort predictions, and F and M use their log
#' mean surfaces. In particular, `log_f` and `log_m` are absolute latent states;
#' their deviations are `eta_log_f = log_f - log_mu_F` and
#' `eta_log_m = log_m - log_mu_M` on the corresponding years and age blocks.
#'
#' The function inspects `dat` to decide which parameters are required and what
#' their dimensions should be. For example, if `dat$F_settings$process == "ar1"`
#' it initializes a 2-vector `logit_phi_f`; if `dat$F_settings$mu_form` is not
#' `NULL` it creates a coefficient vector `log_mu_f` of length
#' `ncol(dat$F_modmat)`, and so on.
#'
#' Numeric parameters are initialized at `0`, except `dq` (initialized at `0.05`)
#' and `log_m`, which starts at its log mean surface so initial M-process
#' residuals are zero. These are starting values, not priors. Matrices are created
#' with appropriate `dimnames` (`year × age` or `year × age_block`).
#'
#' **Created elements (when applicable) include:**
#'
#' - **Recruitment & variability**
#'   - `log_r0` (fixed first-year log recruitment, always present)
#'   - `log_r` (random states for `dat$years[-1]`, length `length(dat$years) - 1`)
#'   - `log_sd_r`
#'
#' - **Initial older-age abundance (independent of the N process)**
#'   - `log_n0`: realized log abundance at initial ages `ages[-1]`, named by age;
#'     absent for `init = "exp"`, fixed for `"free"`, random for `"random"`
#'   - `log_sd_n0`: separate IID initial-age residual SD, only for `init = "random"`
#'   - `eta_log_n0` is calculated internally from adjacent initial log states
#'     and first-year mortality; `log_n0` itself is not a deviation.
#'   - The default `"exp"` initializer is parsimonious survivorship from
#'     `log_r0`, using first-year Z without an equilibrium plus group.
#'     See [prepare_tam()] for the choices and weak-identification warning.
#'
#' - **Abundance states and process variability (N)**
#'   - `log_sd_n` (if `dat$N_settings$process != "off"`)
#'   - `logit_phi_n` length 2 (if `process == "ar1"`)
#'   - `log_n` matrix (`year[-1]` × `age[-1]`) if `process != "off"`
#'
#' - **Fishing mortality (F)**
#'   - `log_sd_f`
#'   - `logit_phi_f` length 2 (if `process == "ar1"`)
#'   - `log_mu_f` coefficients (length `ncol(dat$F_modmat)`) if a mean structure was supplied
#'   - `log_f` matrix (historical `year` × `age`; no projection rows)
#'
#' - **Natural mortality (M)**
#'   - `log_sd_m` (if `dat$M_settings$process != "off"`)
#'   - `logit_phi_m` length 2 (if `process == "ar1"`)
#'   - `mu_m` coefficients (length `ncol(dat$M_modmat)`) if a mean structure was supplied; these act on log-\eqn{M} but are named without the `log_` prefix because they may be positive or negative
#'   - `log_m` matrix (`M_settings$years` × `age_block`) if `process != "off"`, with
#'     `age_block = levels(dat$M_settings$age_blocks)`
#'
#' - **Observation model**
#'   - `log_sd_catch` (length `ncol(dat$sd_catch_modmat)`) adjusting any supplied SDs
#'   - `log_sd_index` (length `ncol(dat$sd_index_modmat)`) adjusting any supplied SDs
#'   - Gaussian observation-SD terms add signed `eta_sd_catch_*` or
#'     `eta_sd_index_*` log-SD effects, initially zero, with process SD parameters
#'     `log_sd_catch_*`/`log_sd_index_*` initially `log(0.1)` unless supplied.
#'     Estimated AR1 correlations start at 0.5. These process SDs measure
#'     variation in log observation SD, not observation SD itself.
#'   - `log_q` for `q_link = "log"`, or `logit_q` for `q_link = "logit"`
#'     (length `ncol(dat$q_modmat)`). Zero coefficients start q at 1 or 0.5,
#'     respectively, before any [mono()] increments.
#'   - `dq` only for [mono()] terms (length `ncol(dat$q_mono_modmat)`):
#'     fixed non-negative increments on the selected q-link scale, initialized to `0.05`
#'     for a nearly flat curve and bounded below by zero in [fit_tam()]. Names
#'     identify transitions and groups; these are not absolute q levels.
#'   - Structured Gaussian q terms have signed `eta_q_*` states, initially zero,
#'     and estimated `log_sd_q_*` parameters initially `log(0.1)` unless SD is
#'     supplied. Estimated `logit_phi_q_*` starts at correlation 0.5.
#'     RW states omit the zero anchor. Forecast states retain their normalized
#'     process densities and integrate out of the historical likelihood.
#'   - Logistic curves have `q_a50_*` midpoints and positive `log_q_slope_*`
#'     coefficients, initialized from each group's observed coordinate range.
#'     Names identify terms, groups and unique states for warm starts.
#'   - `missing` vector of length `sum(dat$fill_missing_map)` (placeholders for
#'     imputed `log_obs`, if any observation type is set to fill missing values)
#'
#' All scalar SD parameters are on the log scale, and AR(1) parameters are on
#' the logit scale (later transformed by `plogis()` in the likelihood).
#' Temporal RW processes have no correlation parameters. Their first process
#' row supplies a conditional starting state during simulation and has no
#' increment penalty. See [fit_tam()] for mean coefficients or SDs held fixed
#' when they have no effect on the RW likelihood.
#'
#' @param dat A data list returned by [prepare_tam()], containing design matrices,
#'   settings, and observation mappings. The shapes and presence/absence of
#'   parameters depend on elements in `dat` (e.g., `F_modmat`, `M_modmat`,
#'   `q_modmat`, `sd_catch_modmat`, `sd_index_modmat`, `N_settings`, `F_settings`,
#'   `M_settings`, `catch_settings`, `index_settings`, `years`, `ages`, and
#'   `M_settings$age_blocks`).
#'
#' @return
#' A named list of initialized parameters suitable to pass to the TAM objective
#' function, with elements as described in **Details**. All numeric entries are
#' initialized to `0` except `dq` (`0.05`) and `log_m` (its log mean). Matrices have
#' informative `dimnames`.
#'
#' @example inst/examples/example_dat_default.R
#' @examples
#' par <- make_par(dat)
#' str(par)
#'
#' @seealso [prepare_tam()], [fit_tam()], [tinyAM-model]
#' @export
make_par <- function(dat) {

  par <- list()
  par$log_r0 <- 0
  if (dat$N_settings$init != "exp") {
    par$log_n0 <- setNames(numeric(length(dat$ages) - 1L), as.character(dat$ages[-1]))
  }
  if (dat$N_settings$init == "random") {
    par$log_sd_n0 <- 0
  }
  if (is.null(dat$rec$sd)) par$log_sd_r <- 0
  if (ncol(dat$rec$matrix)) par$rec_beta <- setNames(numeric(ncol(dat$rec$matrix)), colnames(dat$rec$matrix))
  if (dat$rec$type == "ar1" && is.null(dat$rec$phi)) par$logit_phi_r <- c(year = 0)
  par$log_sd_f <- 0
  if (!is.null(dat$F_settings$mu_form)) {
    par$log_mu_f <- numeric(ncol(dat$F_modmat))
    names(par$log_mu_f) <- colnames(dat$F_modmat)
  }
  if (dat$N_settings$process != "off") {
    par$log_sd_n <- 0
  }
  if (dat$M_settings$process != "off") {
    par$log_sd_m <- 0
  }
  if (!is.null(dat$M_settings$mu_form)) {
    par$mu_m <- numeric(ncol(dat$M_modmat))
    names(par$mu_m) <- colnames(dat$M_modmat)
  }
  if (dat$N_settings$process == "ar1") {
    par$logit_phi_n <- c("age" = 0, "year" = 0)
  }
  if (dat$F_settings$process == "ar1") {
    par$logit_phi_f <- c("age" = 0, "year" = 0)
  }
  if (dat$M_settings$process == "ar1") {
    par$logit_phi_m <- c("age" = 0, "year" = 0)
  }
  par$log_sd_catch <- numeric(ncol(dat$sd_catch_modmat))
  names(par$log_sd_catch) <- colnames(dat$sd_catch_modmat)
  par$log_sd_index <- numeric(ncol(dat$sd_index_modmat))
  names(par$log_sd_index) <- colnames(dat$sd_index_modmat)
  q_par <- if (identical(dat$index_settings$q_link, "logit")) "logit_q" else "log_q"
  par[[q_par]] <- setNames(numeric(ncol(dat$q_modmat)), colnames(dat$q_modmat))
  if (!is.null(dat$q_mono_modmat)) {
    # Small positive link-scale steps start the monotonic curve close to flat.
    par$dq <- setNames(rep(0.05, ncol(dat$q_mono_modmat)),
                       colnames(dat$q_mono_modmat))
  }
  par <- c(par, .q_term_parameters(.formula_terms(dat)))

  if (dat$any_fill_missing) {
    par$missing <- numeric(sum(dat$fill_missing_map))
  }

  par$log_r <- setNames(numeric(length(dat$rec$eligible)), as.character(dat$years[dat$rec$eligible]))
  if (length(dat$rec$boundary) > 1L) {
    par$log_r_init <- setNames(numeric(length(dat$rec$boundary) - 1L), as.character(dat$years[dat$rec$boundary[-1L]]))
  }
  if (!is.null(dat$rec$curve)) {
    par$log_sr_alpha <- par$log_sr_beta <- 0
  }
  if (dat$N_settings$process != "off") {
    par$log_n <- matrix(0, nrow = length(dat$years) - 1L, ncol = length(dat$ages) - 1,
                        dimnames = list(year = dat$years[-1], age = dat$ages[-1]))
  }
  if (dat$M_settings$process != "off") {
    par$log_m <- matrix(0, nrow = length(dat$M_settings$years), ncol = nlevels(dat$M_settings$age_blocks),
                        dimnames = list(year = dat$M_settings$years, age_block = levels(dat$M_settings$age_blocks)))
  }
  par$log_f <- matrix(0, nrow = sum(!dat$is_proj), ncol = length(dat$ages),
                      dimnames = list(year = dat$years[!dat$is_proj], age = dat$ages))

  if (dat$M_settings$process != "off") {
    mu_m_init <- if (is.null(par$mu_m)) dat$mu_m else par$mu_m
    log_mu_M <- matrix(NA, length(dat$years), length(dat$ages),
                       dimnames = list(year = dat$years, age = dat$ages))
    log_mu_M[] <- dat$log_mu_supplied_m + drop(dat$M_modmat %*% mu_m_init)
    par$log_m[] <- log_mu_M[rownames(par$log_m), dat$M_settings$age_block_start, drop = FALSE]
  }

  ## Check for consistent mu M values within age blocks and abort if values are not constant within each block
  if (dat$M_settings$process != "off") {
    getAll(par, dat)
    # Every mean component must be constant within a shared absolute M state.
    mean_parts <- cbind(log_mu_supplied_m, M_modmat)
    for (term in dat$M_terms) {
      mean_parts <- cbind(mean_parts, .q_term_design(term, nrow(dat$obs$weight)))
    }
    for(b in levels(M_settings$age_blocks)) {
      if (sum(M_settings$age_blocks == b) > 1) {
        ia <- names(M_settings$age_blocks)[M_settings$age_blocks == b]
        varies <- vapply(seq_len(ncol(mean_parts)), function(j) {
          surface <- matrix(mean_parts[, j], length(years), length(ages),
                            dimnames = list(year = years, age = ages))
          bmu <- surface[, ia, drop = FALSE]
          any(bmu != bmu[, 1])
        }, logical(1))
        if (any(varies)) {
          cli::cli_abort(c("M mean structure varies within M age_blocks. ",
                           "x" = "mu_form and mu_supplied must be constant within each M process age block. Use finer age_breaks for age-varying M."))
        }
      }
    }
  }

  if (!is.null(dat$rec$curve)) {
    par <- .initialize_rec_curve(par, dat)
  }
  par

}


