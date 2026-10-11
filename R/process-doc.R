#' Mortality changes, age-specific process variation and spawning time
#'
#' @description
#' Neighbouring ages can experience similar annual changes in fishing mortality.
#' `F_settings$process = "cor_rw"` models those changes together, without
#' forcing F to return to an average level. `"rw"` retains independent changes.
#'
#' `sd_form` in N, F and M settings describes how much process variation is
#' allowed at each age. For example, `~ age_group` estimates shared SDs for
#' user-defined age groups stored in `data$weight`; it does not pool abundance
#' or mortality states. Start with few groups and examine their uncertainty.
#' Ordinary numeric and categorical age-based terms are supported. Temporal
#' SDs and random-effect SD terms are not enabled for these latent processes.
#'
#' `ssb_settings = list(spawn_time = 0.25)` evaluates mature biomass after
#' one quarter of annual mortality. This SSB is also used by stock–recruit
#' relationships. The default zero retains start-of-year SSB. This option
#' assumes constant mortality during the year, using the supplied weights and
#' maturity; it does not model seasonal growth or spawning schedules.
#'
#' @details
#' **Correlated increments.** Let \eqn{e_{t,a}=\log F_{t,a}-\log\mu_{F,t,a}}.
#' Independent annual vectors \eqn{e_t-e_{t-1}} have covariance
#' \eqn{\Sigma_{ab}=\sigma_a\sigma_b\rho^{|a-b|}}, where \eqn{-1<\rho<1}.
#' The first residual row is unpenalized. SD is the marginal annual-increment
#' SD, and the independent RW is recovered at \eqn{\rho=0}.
#'
#' **SD formulas.** \eqn{\log\sigma_a=X_a\beta}. Designs must be finite,
#' identifiable, constant across years at each active age, and constant within
#' each fitted M state block. N uses destination ages after recruitment;
#' recruitment and initial-abundance SDs remain separate. IID SDs describe
#' independent residuals and RW SDs describe annual changes. Existing
#' stationary AR1 uses innovation-scale SDs, with marginal SD
#' \eqn{\sigma_a/\sqrt{(1-\phi_a^2)(1-\phi_t^2)}}.
#' Do not compare these interpretations without accounting for correlation.
#'
#' Default `~ 1` retains scalar `log_sd_*` parameters. Other formulas estimate
#' `sd_beta_*` coefficients on the log-SD predictor scale. These coefficients
#' can be negative; they are not individually SDs. `fit$pop$sd_F`, `sd_N` and
#' `sd_M` contain derived positive SD profiles and RTMB confidence intervals
#' for nondefault formulas. For M, the reported ages identify state-block starts.
#'
#' **Spawning biomass.** \eqn{SSB_t=\sum_a N_{t,a}W_{t,a}P_{t,a}e^{-\tau Z_{t,a}}}.
#' The supplied \eqn{\tau\in[0,1]} applies to F and M in every historical and
#' projection year. The chronological population path constructs plus-group
#' biology before calculating SSB and its recruitment contribution. Abundance
#' and total biomass remain start-of-year quantities. Comparison definitions
#' must match age range, biological inputs and spawning time.
#'
#' Numerical convergence does not establish that every age-specific SD is
#' supported. Inspect [check_tam()], intervals and simpler sharing structures,
#' especially when estimating M variability from limited mortality information.
#'
#' @examples
#' data <- cod_obs
#' data$weight$age_group <- factor(ifelse(data$weight$age < 6, "young", "older"))
#' dat <- prepare_tam(data,
#'   F_settings = list(process = "cor_rw", sd_form = ~ age_group),
#'   N_settings = list(process = "iid", sd_form = ~ 1),
#'   ssb_settings = list(spawn_time = 0.25))
#' make_par(dat)$sd_beta_f
#'
#' @name process_variation
#' @seealso [prepare_tam()], [tinyAM-model], [recruitment_formulas], [check_tam()]
NULL
