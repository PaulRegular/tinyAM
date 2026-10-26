source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1992:2025
ages <- 6:19
surveys <- c("Ecosystem survey", "Winter bottom trawl")

inputs <- source$inputs
keep <- (
  inputs$type == "catch" &
    inputs$measure == "numbers_at_age" &
    inputs$fleet == "All fisheries (source total)"
) | (
  inputs$type == "index" &
    inputs$measure == "numbers_at_age" &
    inputs$survey %in% surveys
) | (
  inputs$type == "weight" &
    inputs$measure == "weight_at_age" &
    inputs$age %in% ages
) | (
  inputs$type == "maturity" &
    inputs$measure == "maturity_at_age" &
    inputs$age %in% ages
) | (
  inputs$type == "M" &
    inputs$measure == "natural_mortality_at_age" &
    inputs$age %in% ages
)
inputs <- inputs[!is.na(keep) & keep, , drop = FALSE]

sampling_times <- c(
  "Ecosystem survey" = 0.75,
  "Winter bottom trawl" = 0.12
)
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  sampling_times = sampling_times,
  surveys = surveys,
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0

## Comparison and background ----
age_plus_group <- 19
source_outputs <- source$outputs
source_outputs <- source_outputs[
  !(source_outputs$measure == "fishing_mortality_at_age" &
      !is.na(source_outputs$fleet)),
  , drop = FALSE
]
n_age <- source_outputs[
  source_outputs$measure == "numbers_at_age" &
    source_outputs$year %in% years &
    source_outputs$age %in% ages,
  , drop = FALSE
]
stock_weight <- inputs[
  inputs$type == "weight" & inputs$measure == "weight_at_age",
  c("age", "value"),
  drop = FALSE
]
maturity <- inputs[
  inputs$type == "maturity" & inputs$measure == "maturity_at_age",
  c("age", "value"),
  drop = FALSE
]
ssb_at_age <- n_age
ssb_at_age$value <- n_age$value *
  stock_weight$value[match(n_age$age, stock_weight$age)] *
  maturity$value[match(n_age$age, maturity$age)]
ssb_at_age$type <- "biomass"
ssb_at_age$measure <- "mature_biomass_at_age"
ssb_at_age$unit <- "tonnes"
ssb_at_age$se <- ssb_at_age$lwr <- ssb_at_age$upr <- NA_real_
ssb_at_age$source_type <- "derived_source_output"
ssb_at_age$source_reference <- paste(
  source$assessment$assessment_url,
  "Chapter 6, Tables 6.12, 6.14, 6.20a"
)
ssb_at_age$notes <- paste(
  "Common-definition SSB contribution from reported age-6+ N, fixed stock",
  "weights and model-input maturity. The 19+ group uses its published mean",
  "weight and maturity; the source's hidden older-age detail is unavailable."
)
comparison_outputs <- rbind(source_outputs, ssb_at_age)
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  N = list(
    status = "approximate",
    definition = "Beginning-of-year abundance at ages 6-19+",
    reason = "The source models age-2 recruitment and explicit older cohorts; tinyAM starts at age 6 and treats entry to that age range as its recruitment state."
  ),
  F = list(
    status = "approximate",
    definition = "Total source F at ages 6-19+",
    reason = "The source separates demersal and pelagic fishing; tinyAM fits their reported total F with one aggregate fishery."
  ),
  F_bar = list(
    status = "approximate",
    definition = "N-weighted F over ages 12-18",
    reason = "The source reports Fbar over ages 12-18; tinyAM uses the same ages but estimates an aggregate-fishery F surface."
  ),
  ssb = list(
    status = "approximate",
    definition = "Common-definition SSB over ages 6-19+",
    reason = "Ages below 6 are immature in the framework. Published age-6+ biology is used; the 19+ group's hidden age structure causes small differences from native aggregate SSB."
  ),
  recruitment = list(
    status = "non_equivalent",
    definition = "Official age-2 recruitment is outside the fitted age range",
    reason = "The youngest tinyAM state is age 6 and is not labelled as biological age-2 recruitment."
  ),
  abundance = list(
    status = "non_equivalent",
    definition = "Total stock abundance includes ages below 6",
    reason = "The tinyAM model only represents ages 6-19+."
  ),
  biomass = list(
    status = "non_equivalent",
    definition = "Total stock biomass includes ages below 6",
    reason = "The tinyAM model only represents ages 6-19+."
  ),
  M = list(
    status = "matched",
    definition = "Fixed input M = 0.05 per year",
    reason = "This is an input consistency check, not an estimated-mortality comparison."
  )
)

background <- c(
  "### Northeast Arctic beaked redfish: JRN-AFWG 2026 SCAA",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1992-2025, age-2 recruitment and population ages 2-19+ with older cohorts. | 1992-2025, ages 6-19+. | The published fixed weights and maturity begin at age 6; the source framework treats ages below 6 as immature. The age-6 entry is not biological age-2 recruitment. |",
  "| N | AR1 age-2 recruitment, deterministic cohort survival and explicit older-age states. | AR1 age-6 entry and free initial abundance at ages 7-19+. | Age-2 recruitment and survival through ages 2-5 are not estimated by tinyAM; the shifted age-6 entry approximates that boundary. |",
  "| F | Separate demersal and pelagic annual components and selectivity; fixed M = 0.05. | One aggregate-fishery AR1 F process; the same fixed M on ages 6-19+. | The fishery-specific processes are collapsed; only the official total catch-at-age is used once. |",
  "| Catch | Fleet catch-at-age, catch weights and total catch. | Total catch-at-age for ages 6-19+, with one lognormal error scale. | The pelagic series is a subset of total catch and is omitted to avoid double counting. Catch-weight likelihoods are omitted. |",
  "| Indices | Russian autumn cohort survey, ecosystem and winter age indices, plus deep-pelagic age proportions. | Ecosystem and winter number indices for ages 6-15, with survey-specific log-linear age effects on q. | Cohort-based Russian observations and unavailable deep-pelagic proportions are omitted. The framework describes high selectivity for young ages and lower selectivity at older ages; the tinyAM trend allows a declining q with age. |",
  "| Survey q/selectivity | Survey scaling and age-specific selectivity; the ecosystem and winter surveys have high selectivity for young fish that declines at older ages. | Separate survey intercepts and age slopes on log q. | This is a compact approximation to the source's survey-specific age-selectivity curves. |",
  "| Biology and SSB | Fixed M, stock weights and model-input maturity; age-2 to 5 are immature. | Published fixed biology for ages 6-19+, with January-start SSB. | The 19+ group uses published mean weight and maturity; native older-cohort detail is unavailable. |",
  "",
  "Only ages 6-19+ are fitted. The source's age-2 recruitment is not compared to tinyAM's age-6 entry. SSB, age-specific N and F, and Fbar are common-definition approximations; total abundance and total biomass are not comparable.",
  "",
  "Package starting values are tried first. A single retry uses aligned official age-6+ N and total F only if numerical checks fail."
)

## Model ----
index_settings <- list(
  q_form = ~ 0 + survey + survey:age,
  sd_form = ~ 0 + survey
)

fit_model <- function(...) tinyAM::fit_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "ar1", mu_form = ~ 1),
  F_settings = list(process = "ar1", mu_form = ~ 0 + factor(age), mean_ages = 12:18),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ I(0.05)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = index_settings,
  silent = silent,
  ...
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      n <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      total_f <- source$outputs[
        source$outputs$measure == "fishing_mortality_at_age" &
          is.na(source$outputs$fleet),
        , drop = FALSE
      ]
      f <- .translation_start_surface(
        total_f, "mortality", "fishing_mortality_at_age", years, ages
      )
      list(
        log_r0 = log(n[1, 1]),
        log_r = log(n[-1, 1]),
        log_n0 = log(n[1, -1]),
        log_f = log(f)
      )
    },
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache",
      source$assessment$assessment_id,
      "initialization_attempts.rds"
    )
  )
}
