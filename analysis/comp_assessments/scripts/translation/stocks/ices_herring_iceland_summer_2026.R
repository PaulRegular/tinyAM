## Observations ----
years <- 1980:2023
ages <- 2:15
survey <- "IS-Her-Aco-4Q/1Q"

inputs <- source$inputs
inputs$year_basis[!is.na(inputs$year)] <- "calendar_year"
inputs <- inputs[!(
  (inputs$type == "catch" & inputs$measure == "numbers_at_age" & inputs$value <= 0) |
    (inputs$type == "index" & inputs$measure == "numbers_at_age" & inputs$value <= 0)
), , drop = FALSE]

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = survey,
  sampling_times = stats::setNames(0, survey),
  assumptions = source$assumptions
)

obs$pop_covariates$F_age_block <- factor(ifelse(
  obs$pop_covariates$age <= 6,
  paste0("age", obs$pop_covariates$age),
  "age7plus"
), levels = c(paste0("age", 2:6), "age7plus"))
obs$index$q_age <- factor(ifelse(
  obs$index$age <= 3,
  "age2to3",
  ifelse(obs$index$age <= 6, paste0("age", obs$index$age), "age7plus")
), levels = c("age2to3", "age4", "age5", "age6", "age7plus"))
obs$maturity$spawn_time <- 0

settings <- list(
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(
    process = "ar1", mu_form = ~0 + F_age_block, mean_ages = 5:10
  ),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_age, sd_form = ~1)
)
prepared <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = years, ages = ages), settings
))

## Model and comparison ----
comparison_outputs <- source$outputs[source$outputs$year %in% years, , drop = FALSE]
age_n <- comparison_outputs[
  comparison_outputs$measure == "numbers_at_age" & comparison_outputs$year %in% years,
  c("year", "age", "value")
]
stock_weight <- source$inputs[
  source$inputs$type == "weight" & source$inputs$measure == "weight_at_age" &
    source$inputs$year %in% years,
  c("year", "age", "value")
]
maturity <- source$inputs[
  source$inputs$type == "maturity" & source$inputs$measure == "maturity_at_age" &
    source$inputs$year %in% years,
  c("year", "age", "value")
]
names(age_n)[[3]] <- "N"
names(stock_weight)[[3]] <- "weight"
names(maturity)[[3]] <- "maturity"
derived <- Reduce(function(x, y) merge(x, y, by = c("year", "age")),
  list(age_n, stock_weight, maturity)
)
derived_output <- function(measure, value, notes) data.frame(
  assessment_id = source$assessment$assessment_id,
  type = "biomass", measure = measure,
  fleet = "", survey = "", sex = "combined", region = "ICES 5.a", season = "",
  year = derived$year, age = derived$age, age_group = "", value = value,
  se = NA_real_, lwr = NA_real_, upr = NA_real_, unit = "tonnes",
  source_type = "derived_official_output",
  source_reference = "Derived from MFRI 2026 tables: stock numbers, stock weights, and maturity-at-age",
  notes = notes, stringsAsFactors = FALSE
)
comparison_outputs <- rbind(
  comparison_outputs,
  derived_output(
    "biomass_at_age", derived$N * derived$weight,
    "Official N-at-age (thousand fish) multiplied by stock weight (kg/fish); no uncertainty reported."
  ),
  derived_output(
    "mature_biomass_at_age", derived$N * derived$weight * derived$maturity,
    "Official N-at-age (thousand fish) multiplied by stock weight and maturity at start of model year; no uncertainty reported."
  )
)
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "approximate",
    definition = "Total abundance derived by summing official N-at-age over ages 2-15.",
    reason = "tinyAM represents the terminal age as a plus group; the 2024 SAM benchmark treats age 15 as a distinct age."
  ),
  biomass = list(
    status = "approximate",
    definition = "Total biomass derived from official N-at-age and the reported stock weights over ages 2-15.",
    reason = "The terminal age treatment differs, and tinyAM uses an annual weight surface on its model-year grid."
  ),
  biomass_at_age = list(
    status = "approximate",
    definition = "Official N-at-age multiplied by the reported stock weight at each age and year.",
    reason = "The terminal age treatment differs; derived rows have no reported uncertainty."
  ),
  ssb_mat = list(
    status = "approximate",
    definition = "Official N-at-age multiplied by reported stock weight and maturity at each age and year.",
    reason = "This start-of-year common definition does not include survival to spawning; the production timing is unknown and derived rows have no uncertainty."
  ),
  ssb = list(
    status = "approximate",
    definition = "Spawning biomass from official N-at-age, reported weights, and maturity at the start of the model year.",
    reason = "The published production-run spawning timing and the current infection-related M surface were not recovered."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Same-year abundance at age 2, in the source-reported units.",
    reason = "The age matches, but tinyAM's recruitment process and initial-state treatment simplify SAM."
  ),
  N = list(
    status = "approximate",
    definition = "Annual abundance-at-age over the common 1980-2023 period.",
    reason = "tinyAM's age-15 state is a plus group, unlike the 2024 benchmark's distinct age-15 class."
  ),
  F = list(
    status = "approximate",
    definition = "Annual fishing mortality at age over the common 1980-2023 period.",
    reason = "The tinyAM process and terminal-age treatment simplify the SAM state structure."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Official Fbar over ages 5-10 versus tinyAM's abundance-weighted Fbar over ages 5-10.",
    reason = "The two Fbar calculations use different weighting conventions."
  ),
  M = list(
    status = "unavailable",
    definition = "The current age-year official M surface is not available for comparison.",
    reason = "The supplied 0.1 per-year value is only baseline M and excludes infection-related mortality."
  ),
  M_bar = list(
    status = "unavailable",
    definition = "The current age-year official M surface is not available for comparison.",
    reason = "The supplied 0.1 per-year value is only baseline M and excludes infection-related mortality."
  )
)
age_plus_group <- NULL

background <- c(
  "### Icelandic summer-spawning herring: Official 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | The machine-readable MFRI assessment tables provide model-year inputs for 1980-2023 and ages 2-15. | Fit that same year-age grid. | Newer raw report catch and acoustic series are not spliced into these tables because their production-model mapping is unresolved. |",
  "| Recruitment and abundance | SAM recruitment is at age 2; its log-N process has a separate age-2 variance and shared ages 3-15 variance. | Use IID age-2 recruitment, IID cohort deviations, and exponential initial abundance. | tinyAM uses a different process and initial-state structure. |",
  "| Fishing mortality | SAM uses an AR(1) age correlation; F means are separate at ages 2-6 and shared at ages 7-15. | Use an AR1 F process with the same mean age groups and Fbar ages 5-10. | tinyAM does not reproduce the complete SAM covariance or variance-sharing structure. |",
  "| Natural mortality | M is 0.1 per year plus infection-related mortality. | Supply the fixed 0.1 baseline only. | The current infection-related age-year M surface was not recovered, so this is not total official M. |",
  "| Catch and acoustic index | Catch-at-age and an acoustic age index are available in the machine-readable tables; the benchmark places the acoustic series at the start of the assessment year. | Keep positive rows and use the reported year labels directly, with survey timing 0. | Zero rows are omitted from tinyAM's lognormal likelihood; cross-age AR(1) observation errors are not represented. |",
  "| q, biology, and derived quantities | q is not estimated at age 2, is age-specific at ages 3-6, and is shared at ages 7-15; weights and maturity are treated as known. | Share q for ages 2-3, estimate q separately at ages 4-6 and share age 7+, and use the reported biology. | tinyAM cannot fix the age-2 q at its official value; its age-15 population state is a plus group. SSB uses start-of-year timing because the current production setting is unavailable. |",
  "| Outputs | Age-specific N and F are available through 2024; aggregate SSB and recruitment are available through 2026. | Compare common years through 2023; derive abundance and biomass from official N-at-age, weights, and maturity. | Age-specific uncertainty and the current infection-M surface are not published in the recovered tables. |",
  "",
  "The model-year labels are mapped directly to tinyAM calendar-year labels without shifting them. The first fit uses package-default starting values. The source tables retain zeros; the analysis omits zero catch and acoustic rows because tinyAM uses lognormal observation likelihoods."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = settings$N_settings,
    rec_settings = settings$rec_settings,
    F_settings = settings$F_settings,
    M_settings = settings$M_settings,
    catch_settings = settings$catch_settings,
    index_settings = settings$index_settings,
    silent = silent
  )
}
