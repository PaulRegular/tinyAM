## Observations ----
years <- 1974:2025
ages <- 0:8

mortality <- source$inputs[
  source$inputs$type == "M" &
    source$inputs$measure == "fixed_mortality_medium",
]
mortality$measure <- "natural_mortality_at_age"
inputs <- source$inputs[
  !(source$inputs$type == "M" & grepl("^fixed_mortality_", source$inputs$measure)),
]
inputs <- rbind(inputs, mortality)
inputs <- inputs[!(inputs$type == "index" & inputs$measure == "numbers_at_age" & inputs$age == 0), ]

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = "BIAS",
  sampling_times = c(BIAS = 0.8),
  assumptions = source$assumptions
)
obs$index$q_key <- factor(obs$index$age)

settings <- list(
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "ar1", mu_form = ~ bh(ssb), sd = 0.5),
  F_settings = list(process = "ar1", mu_form = ~ logistic(age), mean_ages = 3:6),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~1)
)
prepared <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = years, ages = ages), settings
))

## Model and comparison ----
age_plus_group <- 8
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Reported official aggregate SSB compared with tinyAM SSB over ages 0-8.",
    reason = "The official output is an ensemble across three M scenarios; tinyAM uses the medium fixed-M scenario. The source does not provide age-specific population outputs for a shared-biology recalculation."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-0 recruitment in the same calendar year.",
    reason = "The age matches, but the official 2025 value is a stock-recruitment prediction and tinyAM uses a simplified Beverton-Holt AR1 process."
  ),
  M = list(
    status = "non_equivalent",
    definition = "Fixed medium-M scenario copied from the official input surface.",
    reason = "The official assessment reports an ensemble input assumption rather than estimated M outputs; matching values reflect reuse of the same input, not an independent fit comparison."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Official unweighted Fbar over ages 3-6 versus tinyAM abundance-weighted Fbar.",
    reason = "The official report does not provide age-specific F or N needed to construct a matching weighted mean."
  )
)

background <- c(
  "### Central Baltic herring: Official 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1974-2025, ages 0-8+; age 8 is the plus group. | Fit 1974-2025 and ages 0-8+. | The earlier assessment history is not modeled. |",
  "| Recruitment and abundance | SS3 Beverton-Holt recruitment with autoregression; process SD fixed at 0.5. | Use the same curve family and fixed SD with an estimated AR1 correlation; keep age-0 recruitment. | tinyAM has a simplified stock-recruit and initial-population structure. |",
  "| Fishing mortality | One commercial fleet; age selectivity is time-invariant with a logistic older-age pattern. | Use an AR1 F process and a logistic age mean. | tinyAM does not reproduce the SS3 likelihood or all selectivity parameters. |",
  "| Natural mortality | Advice is an ensemble of three fixed-M scenarios. | Use the medium fixed-M surface. | One scenario does not reproduce the official ensemble. |",
  "| Catch and survey | Catch-at-age and BIAS acoustic index; the survey reports October and has age-specific selectivity. | Retain catch and BIAS age values, estimate age-specific q, and use timing 0.8 as an October approximation. | Age-composition sample sizes and annual acoustic CVs were not recovered; the exact survey day is unknown. |",
  "| Biology and outputs | Annual weights, maturity, aggregate SSB and age-0 recruitment are reported. | Use the tabulated source biology and compare SSB and recruitment. | Official age-specific N, F and M outputs are unavailable, so those comparisons cannot be made. |",
  "",
  "The medium-M scenario is selected as a transparent single-model approximation to the official three-model ensemble. The BIAS timing value of 0.8 is an analysis assumption for October, not a value reported by ICES. Age-zero acoustic values are kept in the canonical database and omitted from the lognormal observation table because they are zero. The displayed M comparison is an input-consistency check, not an independent estimate. The fit uses package-default starting values."
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
