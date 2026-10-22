## Observations ----
years <- 1981:2023
ages <- 1:8
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("BTS+COAST", "BSAS")

inputs_for_obs <- source$inputs
source_weight <- inputs_for_obs[
  inputs_for_obs$type == "weight" & inputs_for_obs$measure == "weight_at_age",
  , drop = FALSE
]
fitted_weight <- source$outputs[
  source$outputs$type == "biology" & source$outputs$measure == "weight_at_age",
  , drop = FALSE
]
weight_grid <- expand.grid(year = years, age = ages)
key <- function(x) paste(x$year, x$age, sep = "\034")
missing_weight <- weight_grid[!key(weight_grid) %in% key(source_weight), , drop = FALSE]
fitted_index <- match(key(missing_weight), key(fitted_weight))
if (anyNA(fitted_index)) {
  cli::cli_abort("The SAM-fitted stock-weight surface does not fill every missing source-weight cell.")
}
weight_supplement <- source_weight[rep(1L, nrow(missing_weight)), , drop = FALSE]
weight_supplement$year <- missing_weight$year
weight_supplement$year_basis <- "calendar_year"
weight_supplement$age <- missing_weight$age
weight_supplement$value <- fitted_weight$value[fitted_index]
weight_supplement$source_type <- "derived_source_output"
weight_supplement$source_reference <- paste(
  unique(fitted_weight$source_reference), collapse = "; "
)
weight_supplement$transformation <- paste(
  "For missing year-age cells only, use the official SAM fitted latent stock-weight estimate exp(logSW).",
  "Observed weight inputs are retained unchanged; no canonical input values are replaced."
)
weight_supplement$notes <- "Analysis-only supplement because tinyAM requires a complete fixed weight surface; the fitted SAM surface is an estimated biological quantity."
inputs_for_obs <- rbind(inputs_for_obs, weight_supplement)

obs <- database_to_tam_obs(
  assessment_id,
  inputs_for_obs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("BTS+COAST" = 0.705, BSAS = 0.705),
  assumptions = source$assumptions
)

obs$catch$sd_group <- factor(
  ifelse(obs$catch$age == 1, "age1",
    ifelse(obs$catch$age <= 6, "ages2_6", "ages7_8")
  ),
  levels = c("age1", "ages2_6", "ages7_8")
)

obs$index$q_group <- ifelse(
  obs$index$survey == "BTS+COAST",
  ifelse(obs$index$age <= 5, as.character(obs$index$age), "6+"),
  ifelse(obs$index$age <= 6, as.character(obs$index$age), "7+")
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_group,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_group <- ifelse(
  obs$index$survey == "BTS+COAST",
  ifelse(obs$index$age == 1, "BTS_age1",
    ifelse(obs$index$age <= 5, "BTS_ages2_5", "BTS_ages6_7")
  ),
  ifelse(obs$index$age == 1, "BSAS_age1", "BSAS_ages2_8")
)
obs$index$sd_group <- factor(obs$index$sd_group)
obs$pop_covariates$F_sd_group <- factor(
  ifelse(obs$pop_covariates$age == 1, "age1",
    ifelse(obs$pop_covariates$age == 2, "age2",
      ifelse(obs$pop_covariates$age <= 4, "ages3_4", "ages5_8")
    )
  ),
  levels = c("age1", "age2", "ages3_4", "ages5_8")
)
obs$maturity$spawn_time <- 0

## Comparison and fit ----
age_plus_group <- 8L
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_outputs <- source$outputs
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total abundance summed over ages 1-8+ in the same year.",
    reason = "The source age-specific N and modeled age range are available; uncertainty for the sum is not propagated."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass over ages 1-8+ using the completed fixed weight surface shared by the common-definition calculations.",
    reason = "Common-definition biomass applies official N to the same completed weight surface used by tinyAM; it is not SAM native biomass from its latent weight estimates."
  ),
  ssb = list(
    status = "matched",
    definition = "Start-of-year mature biomass over ages 1-8+ using official N and the shared translated weight and maturity surfaces.",
    reason = "This common-definition value applies official N to the same weight and maturity surfaces used by tinyAM; zero pre-spawning F/M preserves the start-of-year definition, but the value is not native SAM SSB."
  ),
  F_bar = list(
    status = "matched",
    definition = "N-weighted F over ages 2-6, matching tinyAM's population-weighted mean.",
    reason = "The official age-specific N and F surfaces allow calculation under the same age range and weighting."
  ),
  M_bar = list(
    status = "matched",
    definition = "N-weighted fixed M over ages 2-6.",
    reason = "M is fixed to the same source age-specific input in tinyAM; this is an input consistency check, not an estimated fit comparison."
  )
)

background <- c(
  "### North Sea turbot: 2025 ICES SAM benchmark model",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | The final 2025 benchmark run M9 covers 1981-2023 at ages 1-8+, with age 8 as the plus group. | Fit the same historical years and age groups. | The assessment's fitted range is retained. |",
  "| Catch and surveys | One residual catch-at-age fleet; BTS+COAST ages 1-7 and BSAS ages 1-8. Both indices are sampled at 0.705 of the year. | Retain the same catch and survey observations and timing. | The source index scale is relative and its physical multiplier is not stated; tinyAM estimates q. BTS+COAST residuals use an age AR(1) structure that tinyAM does not represent. |",
  "| Recruitment and abundance | SAM uses a random-walk recruitment model and age-specific log-N process variances, with one shared survival variance for ages 2-8. | Use random-walk recruitment, IID cohort-survival deviations, and exponential initial abundance. | The survival process covariance and initial-state treatment are simplified. |",
  "| Fishing mortality | SAM uses AR(1) correlation across age states, four variance groups, and a shared age-7/8 state; Fbar covers ages 2-6. | Use correlated random-walk F increments across ages, four age-based process SD groups, and population-weighted Fbar ages 2-6. | tinyAM does not tie age-7 and age-8 F states and its correlated-increment process is an approximation to the SAM covariance. |",
  "| Natural mortality and biology | M is fixed to age-varying Lorenzen values. SAM estimates latent stock and catch weight surfaces from the supplied annual observations; maturity is fixed. Pre-spawning F and M fractions are zero. | Use the same fixed M and maturity. Retain observed stock weights and fill the 102 missing year-age cells with fitted M9 stock weights; treat the full surface as fixed. | tinyAM does not estimate SAM's latent weight surfaces; only missing stock-weight cells borrow fitted M9 estimates. |",
  "| Observation error and q | Catch SDs are shared across ages 2-6 and 7-8; index SDs follow their configured age groups. BTS+COAST index residuals have AR(1) age correlation. q is shared over ages 6-7 for BTS+COAST and 7-8 for BSAS. | Match the q and SD sharing groups; use independent index observation residuals. | The BTS+COAST residual age correlation is omitted. |",
  "| Outputs | SAM reports N, F, SSB, recruitment, total biomass and Fbar. M is fixed. | Compare age-specific N/F, SSB, recruitment, biomass, abundance and common-definition Fbar; show fixed M as an input. | Derived annual sums have no propagated uncertainty. Native Fbar is retained separately; comparison uses N-weighted F over ages 2-6. |",
  "",
  "The 2026 ICES summary is stored separately and is not used in this 2025 benchmark fit. The M9 fitted stock-weight surface is used only for year-age cells absent from the canonical observed weight series. Observed stock-weight inputs remain unchanged, and tinyAM treats the completed surface as fixed."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(
      process = "cor_rw", mu_form = NULL,
      sd_form = ~0 + F_sd_group, mean_ages = 2:6
    ),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~0 + sd_group),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + sd_group
    ),
    silent = silent
  )
}


