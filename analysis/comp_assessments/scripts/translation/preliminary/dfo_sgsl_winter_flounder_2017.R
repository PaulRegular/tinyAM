## Observations ----
years <- 1973:2016
ages <- 2:12
survey <- "September RV survey, offshore strata 415-439"
inputs <- source$inputs
inputs <- inputs[!(inputs$type %in% c("catch", "index") &
  inputs$measure != "numbers_at_age"), ]
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, surveys = survey,
  weight_survey = "September RV survey", assumptions = source$assumptions)
obs$pop_covariates$M_group <- factor(ifelse(obs$pop_covariates$age <= 4, "ages2_4", "ages5_plus"))
obs$pop_covariates$M_period <- factor(1973 + 5 * floor((obs$pop_covariates$year - 1973) / 5))
obs$index$q_period <- factor(ifelse(obs$index$year <= 1984, "1973_1984", "1985_2016"))
obs$catch$catch_log_sd <- ifelse(obs$catch$year == 1984, 2, 0.1)

N_settings <- list(process = "off", init = "exp")
rec_settings <- list(process = "ar1", mu_form = ~1, phi = stats::plogis(-4))
F_settings <- list(process = "off", mu_form = ~logistic(age) + factor(year), mean_ages = 5:12)
M_settings <- list(process = "off", mu_form = ~0 + M_group:M_period, mu_supplied = NULL)
catch_settings <- list(sd_form = ~0, sd_supplied = ~catch_log_sd)
index_settings <- list(q_form = ~logistic(age, by = q_period), sd_form = ~1)
age_plus_group <- 12L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "Model 2 reported SSB compared with beginning-year mature biomass using RV-survey weights.",
  reason = "Native spawning fractions and fishery weight schedules are unavailable; only terminal 2016 source estimates are recovered."))

## Background ----
background <- c(
  "### Southern Gulf winter flounder: preliminary comparison with Official 2017 Model 2", "",
  print_sources(source$assessment), "",
  "The research document's main numerical results concern Model 2, although Model 7 has the lowest AIC. The comparison retains that distinction; it does not assert an unambiguous selected run.", "",
  "| Component | Official Model 2 | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Coverage | 1973–2016, ages 2–12+. | Same years and ages. | No later interim survey indicators are treated as a population-model fit. |",
  "| N and recruitment | Deterministic survival, AR1 recruitment with logit(phi) fixed at −4, equilibrium prehistory. | Deterministic survival and the same recruitment phi; exponential initial abundance. | Initial recruitment history and equilibrium F differ. The supplied phi is about 0.018, irrespective of the report's random-walk label. |",
  "| F | One logistic fishery selectivity with annual fishing intensity. | Static logistic age effect plus annual log-F levels. | Source catch-equation solution and catch-total likelihood are not reproduced. |",
  "| M | Two age groups, 2–4 and 5+, estimated in five-year blocks with priors. | Estimate one M per matching age/year cell. | Source priors are omitted; no arbitrary fixed-M baseline is substituted. |",
  "| Catch | Biomass and number compositions, catch-total log-SD 0.1 except 2 in 1984. | Fit the published count surface with those fixed log-SDs. | Applying total-catch SDs independently to age counts is an explicit approximation; it does not reproduce composition errors or their dependence. Published rounded zero counts remain zero. |",
  "| Survey | RV biomass and composition; logistic selectivity changes after 1984; shared q has a prior. | Fit the original trawlable age counts at fraction 0.71, with separate logistic curves and a shared baseline q. | No duplicate derived composition or unresolved-scale biomass index is added. The q prior is omitted, an important identifiability limitation. |",
  "| Biology and outputs | Annual RV weights and maturity; only terminal Model 2 numerical population outputs recovered. | Retain RV weights and maturity; compare terminal abundance and approximate mature biomass. | Fishery weights and spawning timing are unavailable. A single source year cannot establish trend agreement or historical uncertainty. |",
  "", "Package-default initialization is used. No aligned official age-state estimates are available for a retry. The source model itself reports weak q information without its prior; convergence is not assumed."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(function() tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    silent = silent), retry_start = function() {
      cli::cli_abort("Official N/F age states are unavailable; no aligned source-state retry is possible.")
    }, attempt_file = file.path(.assessment_cache(source$assessment$assessment_id,
      translation_tier), "initialization_attempts.rds"))
}
