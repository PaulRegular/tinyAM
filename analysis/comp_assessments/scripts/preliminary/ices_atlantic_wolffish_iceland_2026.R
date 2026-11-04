## Observations ----
years <- 1998:2025
ages <- 4:16
inputs <- source$inputs
relative <- inputs$type == "index" & inputs$measure == "numbers_at_age"
inputs$transformation[relative] <- paste("Uncalibrated age-specific index retained without rescaling; source unit:",
  inputs$unit[relative])
inputs$unit[relative] <- "native survey index"
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, assumptions = source$assumptions)

N_settings <- list(process = "iid", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 10:15)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + survey)
age_plus_group <- 16L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(ssb = list(
  status = "approximate",
  definition = "Common-definition mature biomass: source N with shared stock weights and female maturity ogives over ages 4–16+.",
  reason = "The separate official SSB-weight surface, sex convention and spawning-time definition are not fully recoverable; this is not a direct comparison with native female SSB."))

## Background ----
background <- c(
  "### Icelandic Atlantic wolffish: preliminary comparison with the Official 2026 assessment", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years | 1979–2025; aggregate landings replace unavailable catch-at-age in 1981–1997. | 1998–2025, with no forecast observations. | This is the contiguous recent period with catch-at-age and complete annual biology. Earlier landings are not allocated using estimated N. |",
  "| Ages | Ages 4–16+, with recruitment at age 4. | The same ages and plus group. | No age groups are split or extrapolated. |",
  "| N | SAM; full current process and initialization settings were not recovered. | Random-walk recruitment, IID cohort-survival deviations and free initial abundance. | The covariance and initial-state treatment are approximations. |",
  "| F | Annual mortality-at-age; exact process settings were not recovered. | Independent temporal random walks by age. | No particular SAM covariance or age sharing is asserted. |",
  "| M | Fixed at 0.15 per year at every age. | The same supplied M. | Agreement is by construction. |",
  "| Catch | Reported commercial numbers-at-age in thousand fish. | The same counts with a common lognormal SD. | The source's aggregate/composition likelihood and weights are not reproduced. |",
  "| Index | Spring and autumn age-specific indices on an unresolved native scale. | Retain both surveys, survey-age q, survey-specific SD, and source timing approximations 0.17 and 0.83. | Physical multipliers and exact source likelihood weights remain unknown. The published 2011 autumn row is retained with its source qualification. |",
  "| Weights | Annual stock and catch weights; the separate SSB-weight export is empty. | Use recorded stock weights, converting g/fish to kg/fish. | Missing SSB weights are not invented or replaced in SAIO. |",
  "| Maturity and SSB | Female maturity ogives; native SSB is reported as mature female biomass. | Use the annual ogives and compare a labelled common-definition mature biomass from source N and the same translated biology. | Sex convention, missing SSB weights and unresolved spawning timing prevent claiming native SSB equivalence. Native SSB remains available separately. |",
  "", "The source provides N and F at age without age-specific uncertainty. The 2026 forecast is outside the comparison period. Starts use package defaults; numerical failure permits one aligned source N/F retry."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent)
  fit <- fit_default_first(fit_model, retry_start = function() {
    dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages,
      N_settings = N_settings, rec_settings = rec_settings,
      F_settings = F_settings, M_settings = M_settings,
      catch_settings = catch_settings, index_settings = index_settings)
    source_start_parameters(source, dat)
  }, attempt_file = file.path(.assessment_cache(source$assessment$assessment_id,
    translation_tier, database = source), "initialization_attempts.rds"))
}
