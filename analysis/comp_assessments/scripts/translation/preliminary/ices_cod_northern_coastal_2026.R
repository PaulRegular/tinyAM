## Observations ----
years <- 1995:2025
ages <- 3:10

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = "Coastal swept-area index",
  sampling_times = c("Coastal swept-area index" = 0.8),
  assumptions = source$assumptions
)

obs$index$q_key <- factor(
  ifelse(
    obs$index$age <= 4,
    paste0("age", obs$index$age),
    ifelse(obs$index$age <= 9, "ages5_9", "age10_plus")
  ),
  levels = c("age3", "age4", "ages5_9", "age10_plus")
)
obs$maturity$spawn_time <- 0.8

source_N <- source$outputs[
  source$outputs$measure == "numbers_at_age" &
    source$outputs$year %in% years &
    source$outputs$age %in% ages,
  , drop = FALSE
]
source_weight <- source$inputs[
  source$inputs$type == "weight" &
    source$inputs$measure == "weight_at_age" &
    source$inputs$year %in% years &
    source$inputs$age %in% ages,
  , drop = FALSE
]
state_key <- function(x) paste(x$year, x$age)
weight_index <- match(state_key(source_N), state_key(source_weight))
if (anyNA(weight_index) || anyDuplicated(state_key(source_N)) ||
  anyDuplicated(state_key(source_weight))) {
  stop("Official N and stock weight do not align on the comparison grid.", call. = FALSE)
}
biomass_at_age <- source_N
biomass_at_age$type <- "biomass"
biomass_at_age$measure <- "biomass_at_age"
biomass_at_age$value <- source_N$value * source_weight$value[weight_index]
biomass_at_age$se <- biomass_at_age$lwr <- biomass_at_age$upr <- NA_real_
biomass_at_age$unit <- "tonnes"
biomass_at_age$source_type <- "derived_source_output"
biomass_at_age$source_reference <- "Derived from AFWG 2026 Tables 2.2.13 and 2.2.19"
biomass_at_age$notes <- paste(
  "Point estimate derived as official N-at-age in thousand fish multiplied by",
  "official stock weight in kg per fish; age 10 is the 10+ group.",
  "Age-specific uncertainty is unavailable."
)

## Comparisons and background ----
age_plus_group <- 10L
comparison_outputs <- rbind(source$outputs, biomass_at_age)
comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar", "M_bar")
comparison_definitions <- list(
  abundance = list(
    status = "approximate",
    definition = "Total abundance across ages 3-10, with age 10 representing the source 10+ group.",
    reason = "tinyAM does not expand hidden ages within the source plus group."
  ),
  biomass = list(
    status = "approximate",
    definition = "Biomass across the common modeled ages 3-10, using source stock weights.",
    reason = "Age 2 is excluded to align tinyAM recruitment with the official age-3 recruitment; age 10 is a terminal plus group."
  ),
  biomass_at_age = list(
    status = "matched",
    definition = "Annual biomass at age from official numbers-at-age and stock weight; age 10 is the 10+ group.",
    reason = "Derived point estimate; age-specific N uncertainty is unavailable."
  ),
  ssb = list(
    status = "approximate",
    definition = "Mature biomass over ages 3-10 at survey time 0.8, using source stock weights and maturity.",
    reason = "Age 2 is excluded; the official native SSB over the full age range remains available in the reference."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-3 recruitment in thousand fish.",
    reason = "Both series describe recruitment at age 3; tinyAM uses a random walk and its own initial-state treatment."
  ),
  F = list(
    status = "approximate",
    definition = "Annual fishing mortality at ages 3-10.",
    reason = "The source shares F states at ages 7-9 and uses correlated age innovations; tinyAM models correlated increments without shared states."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Abundance-weighted F over ages 4-8.",
    reason = "This common definition is derived from age-specific F and N; the native SAM output is an arithmetic mean."
  ),
  M = list(
    status = "matched",
    definition = "Fixed annual natural mortality at ages 3-10.",
    reason = "The same values are supplied to tinyAM; this checks input translation, not mortality estimation."
  ),
  M_bar = list(
    status = "matched",
    definition = "Abundance-weighted fixed M over ages 3-10.",
    reason = "The same fixed age-specific M inputs are used; this is an input consistency check."
  )
)

background <- c(
  "### Northern Norwegian coastal cod: official 2026 SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | Catch and population inputs span 1994-2025 at ages 2-10+, with recruitment at age 3. | Fit 1995-2025 at ages 3-10. | Age 2 is excluded so the youngest tinyAM age matches official recruitment age 3. Age 10+ is represented as a terminal age-10 group. |",
  "| Catch and survey | Total catch-at-age and an age-specific swept-area index; both are in thousands of fish. The swept index is sampled at 0.8. | Retain catch and swept-area observations for ages 3-10; use survey timing 0.8. | The aggregate age-2+ acoustic biomass likelihood is not supported by tinyAM and is omitted from the fit; its values remain in the database. |",
  "| Catchability and observation error | Swept q is separate at ages 2, 3, 4 and 10+, shared at ages 5-9. Swept residuals are correlated across ages. | Use the same q age groups and one estimated lognormal SD. | tinyAM does not represent the source cross-age observation correlation. |",
  "| N, recruitment and F | SAM uses log-abundance and log-F random walks; age-2 recruitment variance is separate, ages 3-10 share survival variance, and F states at ages 7-9 are shared. | Use random-walk recruitment and abundance residuals, plus correlated random-walk F increments. | Variance sharing, shared F states, and source initial-state integration are not reproduced. |",
  "| M and biology | Annual M is supplied from the source table; stock weight and maturity are annual. SSB is reported at survey time 0.8. | Retain supplied M, stock weights, maturity, and maturity time 0.8. | M is fixed, not estimated. The dashboard M comparison is an input consistency check. |",
  "| Outputs | The report provides age-specific N and F, age-3 recruitment, SSB, Fbar ages 4-8, total biomass, and q. | Compare available N/F, age-3 recruitment, fixed M, biomass-at-age, SSB, abundance, biomass, and common-definition Fbar/Mbar. | Biomass-at-age is derived for comparison from official N and stock weight; its uncertainty is unavailable. N/F age-specific uncertainty is also unavailable. Native aggregate results are retained separately from common-age comparisons. |",
  "",
  "The current AFWG report names final run NCCN67_AFWG2025_v3. Its cited public TAF repository could not be retrieved, so the exact fitted object and initial-state treatment remain unresolved. The source tables are used as published. The 2000 acoustic tuning value of 35,255 tonnes is retained from the final-run tuning table; the report's acoustic summary table gives 35,243 tonnes. The 2021 acoustic observation is missing from the final tuning series because it was excluded for high uncertainty."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "rw", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(
      process = "cor_rw",
      mu_form = NULL,
      mean_ages = 4:8
    ),
    M_settings = list(
      process = "off",
      mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~1),
    silent = silent
  )
}
