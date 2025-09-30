
#' Initialize parameter list for TAM
#'
#' @description
#' `make_par()` builds a named list of initial values and shapes for all
#' fixed and random-effect parameters used by TAM, based on the structure in
#' a previously constructed `dat` list (see [make_dat()]).
#'
#' @details
#' The function inspects `dat` to decide which parameters are required and what
#' their dimensions should be. For example, if `dat$F_settings$process == "ar1"`
#' it initializes a 2-vector `logit_phi_f`; if `dat$F_settings$mu_form` is not
#' `NULL` it creates a coefficient vector `log_mu_f` of length
#' `ncol(dat$F_modmat)`, and so on.
#'
#' All numeric parameters are initialized at `0`, and all matrices are created
#' with appropriate `dimnames` (`year × age` or `year × age_block`).
#'
#' **Created elements (when applicable) include:**
#'
#' - **Recruitment & variability**
#'   - `log_r0` (only if `dat$N_settings$init_N0`)
#'   - `log_r` (length `length(dat$years)`)
#'   - `log_sd_r`
#'
#' - **Abundance deviations (N)**
#'   - `log_sd_n` (if `dat$N_settings$process != "off"`)
#'   - `logit_phi_n` length 2 (if `process == "ar1"`)
#'   - `log_n` matrix (`year` × `age[-1]`) if `process != "off"`
#'
#' - **Fishing mortality (F)**
#'   - `log_sd_f`
#'   - `logit_phi_f` length 2 (if `process == "ar1"`)
#'   - `log_mu_f` coefficients (length `ncol(dat$F_modmat)`) if a mean structure was supplied
#'   - `log_f` matrix (`year` × `age`)
#'
#' - **Natural mortality (M)**
#'   - `log_sd_m` (if `dat$M_settings$process != "off"`)
#'   - `logit_phi_m` length 2 (if `process == "ar1"`)
#'   - `log_mu_m` coefficients (length `ncol(dat$M_modmat)`) if a mean structure was supplied
#'   - `log_m` matrix (`year[-1]` × `age_block`) if `process != "off"`, with
#'     `age_block = levels(dat$M_settings$age_blocks)`
#'
#' - **Observation model**
#'   - `log_q` (length `ncol(dat$q_modmat)`)
#'   - `log_sd_obs` (length `ncol(dat$sd_obs_modmat)`)
#'   - `missing` vector of length `sum(is.na(dat$log_obs))` (placeholders for
#'     imputed `log_obs`)
#'
#' All scalar SD parameters are on the log scale, and AR(1) parameters are on
#' the logit scale (later transformed by `plogis()` in the likelihood).
#'
#' @param dat A data list returned by [make_dat()], containing design matrices,
#'   settings, and observation mappings. The shapes and presence/absence of
#'   parameters depend on elements in `dat` (e.g., `F_modmat`, `M_modmat`,
#'   `q_modmat`, `sd_obs_modmat`, `N_settings`, `F_settings`, `M_settings`,
#'   `years`, `ages`, and `M_settings$age_blocks`).
#'
#' @return
#' A named list of initialized parameters suitable to pass to the TAM objective
#' function, with elements as described in **Details**. All numeric entries are
#' initialized to `0` (or the correct length filled with `0`), and matrices have
#' informative `dimnames`.
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
#' str(par)
#'
#' @seealso [make_dat()]
#' @export
make_par <- function(dat) {

  par <- list()
  if (dat$N_settings$init_N0) {
    par$log_r0 <- 0
  }
  par$log_sd_r <- 0
  par$log_sd_f <- 0
  if (!is.null(dat$F_settings$mu_form)) {
    par$log_mu_f <- numeric(ncol(dat$F_modmat))
  }
  if (dat$N_settings$process != "off") {
    par$log_sd_n <- 0
  }
  if (dat$M_settings$process != "off") {
    par$log_sd_m <- 0
  }
  if (!is.null(dat$M_settings$mu_form)) {
    par$log_mu_m <- numeric(ncol(dat$M_modmat))
  }
  if (dat$N_settings$process == "ar1") {
    par$logit_phi_n <- c(0, 0)
  }
  if (dat$F_settings$process == "ar1") {
    par$logit_phi_f <- c(0, 0)
  }
  if (dat$M_settings$process == "ar1") {
    par$logit_phi_m <- c(0, 0)
  }
  par$log_q <- numeric(ncol(dat$q_modmat))
  par$log_sd_obs <- numeric(ncol(dat$sd_obs_modmat))

  par$missing <- numeric(sum(is.na(dat$log_obs)))

  par$log_r <- rep(13, length(dat$years)) # numeric(length(dat$years))
  if (dat$N_settings$process != "off") {
    par$log_n <- matrix(11, nrow = length(dat$years), ncol = length(dat$ages) - 1,
                        dimnames = list(year = dat$years, age = dat$ages[-1]))
  }
  if (dat$M_settings$process != "off") {
    par$log_m <- matrix(0, nrow = length(dat$years) - 1, ncol = nlevels(dat$M_settings$age_blocks),
                        dimnames = list(year = dat$years[-1], age_block = levels(dat$M_settings$age_blocks)))
  }
  # par$log_f <- matrix(0, nrow = sum(!dat$is_proj), ncol = length(dat$ages),
  #                     dimnames = list(year = dat$years[!dat$is_proj], age = dat$ages))
  par$log_f <- outer(seq_along(dat$years[!dat$is_proj]), dat$ages, function(y, a) {
    log(pmax(1 - ((a - mean(range(dat$ages))) / (diff(range(dat$ages))/2))^2, 1e-3))
    })

  par

}


