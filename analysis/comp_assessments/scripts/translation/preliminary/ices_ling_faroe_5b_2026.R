## Observations ----
years <- 1996:2025
ages <- 3:12
surveys <- c("Faroe Plateau spring groundfish survey (G1264)",
  "Faroe Plateau summer groundfish survey (G3284)")
inputs <- source$inputs
relative <- inputs$type == "index" & inputs$measure == "relative_abundance_index"
inputs$transformation[relative] <- paste("Age-specific relative index retained without rescaling; source unit:",
  inputs$unit[relative])
inputs$measure[relative] <- "numbers_at_age"
inputs$unit[relative] <- "native survey index"
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, surveys = surveys,
  sampling_times = stats::setNames(c(2 / 12, 7.5 / 12), surveys),
  assumptions = source$assumptions)

N_settings <- list(process = "iid", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 6:10)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + survey)
age_plus_group <- 12L
comparison_aggregates <- c("abundance", "ssb")

## Background ----
background <- c(
  "### Faroes ling: preliminary comparison with the Official 2026 SAM assessment", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1996–2025; ages 3–12+, recruitment at age 3. | Same historical years and ages. | No source years or ages are removed. |",
  "| N | SAM; the native process covariance and initial-state configuration were not recovered. | Random-walk recruitment, IID cohort-survival deviations and free initial abundance. | A documented approximation, not a recovered SAM configuration. |",
  "| F | Annual SAM F states; exact covariance and state sharing remain unresolved. | Independent temporal random walks by age. | No unsupported cross-age correlation or shared state is asserted. |",
  "| M | Fixed at 0.15 per year for all ages. | The same supplied M. | Fixed M agreement is an input check. |",
  "| Catch | Combined catch numbers-at-age. | Retain the reported counts with one lognormal SD. | Exact source likelihood weights are unavailable. |",
  "| Index | Spring ages 4–11 and summer ages 3–11. Gear-change years were downweighted by unknown factors. | Retain both native series, with survey-age q and survey-specific SD. | February/March and August are represented by fractions 2/12 and 7.5/12. No effort division or unavailable downweight is invented. |",
  "| Weights | The benchmark selected catch weights as stock weights; kg/fish is inferred from consistency with reported biomass. | Retain the published annual values on that inferred kg scale. | There is no numerical rescaling beyond the shared unit conversion. |",
  "| Maturity and SSB | Fixed age ogive; report N and biology reproduce aggregate biomass within rounding. | Retain the ogive and compare mature biomass at the start of the year. | Native SSB and arithmetic Fbar remain separate from common-definition comparisons. |",
  "", "Source scale, survey timing and gear-change downweights remain unresolved. Starts come from package defaults; only a numerical failure permits one aligned N/F retry."
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
