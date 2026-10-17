source("analysis/comp_assessments/R/bc_canary_obs.R", local = TRUE)
source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1977:2022
ages <- 1:60
catch <- .canary_female_numbers(source$inputs, "catch")
index <- .canary_female_numbers(source$inputs, "index")
weight <- source$inputs[source$inputs$type == "weight" & source$inputs$sex == "female" &
  !is.na(source$inputs$season) & source$inputs$season == "begin_year", ]
maturity <- source$inputs[source$inputs$measure == "maturity_fecundity_at_age" &
  source$inputs$sex == "female", ]
maturity$value <- maturity$value / weight$value[match(maturity$age, weight$age)]
maturity$type <- "maturity"
maturity$measure <- "maturity_at_age"
maturity$basis <- maturity$unit <- "proportion"
maturity$notes <- "Effective female maturity from source reproductive weights divided by January population weight; preserves the executable's age-10 spawning cutoff. Female abundance is modeled directly, so no additional factor 0.5."
inputs <- rbind(catch, index, weight, maturity)
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, assumptions = source$assumptions)
obs$maturity$spawn_time <- 0
obs$pop_covariates$M_prior_group <- factor("female_stock")

## Background and comparisons ----
age_plus_group <- 60
comparison_aggregates <- "ssb"
comparison_outputs <- source$outputs[source$outputs$sex == "female", ]
comparison_definitions <- list(
  ssb = list(status = "approximate", definition = "January 1 mature female biomass, ages 1-60+",
    reason = "Age 0 is immature. Official annual points and 95% credible limits summarize the public 500-draw subset, not the full 2000-draw posterior; tinyAM reports likelihood-based estimates."),
  M = list(status = "matched", definition = "Female constant instantaneous M at ages 1-60+",
    reason = "Official points are full-posterior medians with native 90% limits. The natural-normal prior and posterior integration are approximated in tinyAM; repeated source M values are one shared estimate, not annual independent estimates."))
background <- c(
  "### British Columbia Canary Rockfish: official 2022 Run24", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1935-2022; age-0 recruitment and 60+; sampled ages 1-60+. | 1977-2022, ages 1-60+. | Begin at the first commercial age sample, retaining the full sampled age range. Earlier history is represented by free initial ages. Beginning-2023 source SSB is boundary context, outside historical agreement summaries. |",
  "| Population | Separate female and male states; deterministic cohort survival; unfished equilibrium with early recruitment deviations. | Female abundance with deterministic survival and free initial ages. | Sex-specific numbers are reconstructed jointly before selecting female observations. No assumed adult female fraction or second factor of 0.5 is applied to SSB. |",
  "| Recruitment | Age-0 Beverton-Holt, estimated steepness with beta prior, log SD 0.9, sum-constrained deviations and varying bias adjustment. | Age-1 Beverton-Holt with IID log residuals, SD 0.9, parent SSB lag one year. | Age-1 entry includes juvenile survival. The steepness prior, early deviations, sum constraint, terminal treatment and source bias adjustment are not reproduced. |",
  "| F | Annual intensity in two fleets, Pope catch approximation, static asymptotic double-normal age selectivity; non-trawl shape fixed. | One F surface with a rising logistic age mean, shared annual random walk and IID residuals; Baranov catches. | Approximate the static ascending shape and common annual intensity. Source curve, fleet separation and catch equation differ. |",
  "| M | Separate constant female and male M, natural-normal prior mean 0.06, SD 0.018, bounds 0.02-0.2. | One female log-M effect centered on log(0.06), fixed log SD 0.3. | Lognormal penalty approximates the source prior CV; bounds and Bayesian integration differ. The source-estimated M is not fixed as an input or used to tune observations. |",
  "| Catch | Two biomass streams; joint-sex trawl age samples in 36 years, 1977-2017; no non-trawl ages. | Pooled biomass allocated using each sampled year's trawl joint-sex age profile and both sexes' catch weights, retaining female numbers. | Assumes the unsampled non-trawl catch shares the trawl profile. Non-trawl removals are about 4% of biomass over 1977-2022. No ages are filled for missing years or 2018-2022, and totals alone do not enter the tinyAM likelihood. Independent lognormal age errors replace the Dirichlet-multinomial and ageing error. |",
  "| Indices and q | Six biomass surveys plus commercial CPUE; three surveys have joint-sex ages; profiled q, asymptotic double-normal selectivity and supplied precision. | Female age indices for QCS, WCVI and NMFS, time 5/12; separate rising logistic q curves and log SDs. | CPUE, HS, WCHG and GIG have no corresponding age compositions and are omitted. Composition covariance, ageing error, source precision and selectivity priors are not reproduced. |",
  "| Biology | Fixed sex-specific growth/length variability/allometry; January female mature biomass. | Original female January weights and effective maturity; source fleet weights for biomass conversion. | Matching SS3 evaluated fixed biology only. Preserve the source first-mature-age cutoff at 10, despite nonzero supplied maturity at 4-9. Survey weights retain the source midyear values with timing 5/12. |",
  "", "The official 2022 catch assumption is 780 tonnes, allocated by 2017-2021 fleet shares. It is not final observed 2022 landings. Canonical records retain all original fleets, indices, sexes and supplied compositions. These derived female observations are an approximation to the official likelihood, not independent new samples.",
  "Official annual SSB is a 500-draw public posterior subset; its uncertainty is a 95% credible interval. Published terminal SSB and M retain their native 90% intervals. Numerical official N/F-at-age and recruitment are unavailable; those comparisons remain blank. Default package starts are attempted first. If numerical checks fail, one retry may scale an initial survivorship profile to official first-year SSB; this is an initialization approximation, not a recovered official N estimate."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ bh(ssb), sd = .9),
  F_settings = list(process = "iid", mu_form = ~ logistic(age) + rw(year)),
  M_settings = list(process = "off", mu_form = ~ 0 + iid(M_prior_group, sd = .3),
                    mu_supplied = ~ I(.06)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(q_form = ~ 0 + survey + logistic(age, by = survey),
                        sd_form = ~ 0 + survey),
  proj_settings = list(n_proj = 1),
  silent = silent, ...
)
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    m <- unique(comparison_outputs$value[comparison_outputs$measure == "natural_mortality_at_age"])
    b <- comparison_outputs$value[comparison_outputs$measure == "SSB" &
      comparison_outputs$year == min(years)] * 1000
    survival <- exp(-m * (ages - min(ages)))
    survival[length(ages)] <- tail(survival, 1) / (1 - exp(-m))
    w <- obs$weight$obs[obs$weight$year == min(years)]
    mat <- obs$maturity$obs[obs$maturity$year == min(years)]
    recruitment <- b / sum(survival * w * mat)
    stopifnot(length(m) == 1L, length(b) == 1L, is.finite(recruitment), recruitment > 0)
    list(log_r0 = log(recruitment),
      log_n0 = setNames(log(recruitment * survival[-1]), ages[-1]),
      log_r = setNames(rep(log(recruitment), length(years) - 1L), years[-1]),
      eta_mu_M_iid_M_prior_group = setNames(log(m/.06), "female_stock"))
  }, attempt_file = file.path("analysis/comp_assessments/results/cache",
    source$assessment$assessment_id, "initialization_attempts.rds"))
}
