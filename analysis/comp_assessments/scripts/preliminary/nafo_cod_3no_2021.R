## Observations ----

assessment_id <- source$assessment$assessment_id[[1]]
years <- 1960:2019
ages <- 3:12
surveys <- c(
  "Canadian spring RV bottom-trawl survey",
  "Canadian autumn RV bottom-trawl survey",
  "Canadian juvenile survey"
)
sampling_times <- c(
  "Canadian spring RV bottom-trawl survey" = 0.38,
  "Canadian autumn RV bottom-trawl survey" = 0.84,
  "Canadian juvenile survey" = 0.67
)

inputs <- source$inputs
inputs$value <- suppressWarnings(as.numeric(as.character(inputs$value)))
inputs <- inputs[!(inputs$type %in% c("catch", "index") &
  (is.na(inputs$value) | inputs$value <= 0)), , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id, inputs,
  years = years, ages = ages,
  surveys = surveys, sampling_times = sampling_times,
  assumptions = source$assumptions
)
obs$index$q_key <- factor(paste(obs$index$survey, obs$index$age, sep = "_"))

age_plus_group <- NULL
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "approximate",
    definition = "Common-definition total abundance over ages 3–12.",
    reason = "The model omits age 2 and its recruitment age is younger than tinyAM's age-3 recruitment state."
  ),
  biomass = list(
    status = "approximate",
    definition = "Common-definition stock biomass over ages 3–12.",
    reason = "The comparison is restricted to the fitted age range."
  ),
  ssb = list(
    status = "approximate",
    definition = "Common-definition mature biomass over ages 3–12 at start of year.",
    reason = "The fit omits age 2; age-specific dynamics and the source ADAPT reconstruction are simplified."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Beginning-of-year abundance at age 3 in the same calendar year.",
    reason = "The model's youngest age is 3; source recruitment is not defined as age 3."
  ),
  M = list(
    status = "approximate",
    definition = "Fixed M=0.2 per year at each fitted age.",
    reason = "The same fixed input is used in both models; this is not an estimated-M comparison."
  )
)

## Model ----

N_settings <- list(process = "off", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 6:9)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption)
catch_settings <- list(sd_form = ~0, sd_supplied = ~ I(0.01))
index_settings <- list(
  q_form = ~ 0 + q_key,
  sd_form = ~0,
  sd_supplied = ~ I(1)
)

background <- c(
  "### NAFO Divisions 3NO cod: preliminary translation of the 2021 ADAPT assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Catch and survey inputs through 2019; ages 2–12; no plus group. | Fit 1960–2019, ages 3–12. | Stock weights begin in 1960; age 2 is omitted to keep a complete stock-weight surface. tinyAM treats age 12 as its terminal group. |",
  "| Catch | Commercial numbers-at-age, treated as exact. | Direct positive catch-at-age with a fixed log-scale SD of 0.01. | tinyAM requires a finite observation SD; 0.01 is a tight approximation to the source's exact catch constraint, not a source-estimated error. |",
  "| Indices and q | Three separate fish-per-tow age indices; q is separate by survey and age; equal weighting on log residuals. | Retain all three series, with survey-age q and one fixed log-residual SD of 1. | Spring, autumn and juvenile timing are represented by approximate seasonal midpoints (0.38, 0.84 and 0.67); the report's survey dates vary by year. Zero indices are omitted from the lognormal likelihood. |",
  "| N, recruitment and F | ADAPT projects cohorts from catch-at-age and estimates terminal survivors; F at age 12 is constrained to mean F at ages 6–9 in 1959–1993. | Deterministic cohort survival, age-3 recruitment random walk, and an age-specific F random walk. | tinyAM's state/process formulation is not the ADAPT solver; it does not reproduce the historical age-12 F constraint. |",
  "| M and biology | M fixed at 0.2 per year; source stock weights, female maturity and catch weights are reported. | Retain fixed M and source biology over the fitted ages. | Catch weights are retained as source data; only the stock weights enter population biomass. |",
  "| Comparison | Annual N-at-age, F-at-age, total abundance, biomass and SSB are available, with population estimates through 2020. | Compare common ages 3–12 over input years 1960–2019. | Age 2 and the 2020 population-only endpoint are outside the fit; N/F/SSB remain approximate comparisons. |",
  "",
  "The source reports that absolute ADAPT estimates are unreliable and the model is intended to track trends. Treat agreement as descriptive, not validation of absolute population scale."
)

fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent
  )
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      dat <- tinyAM::prepare_tam(
        data = obs, years = years, ages = ages,
        N_settings = N_settings, rec_settings = rec_settings,
        F_settings = F_settings, M_settings = M_settings,
        catch_settings = catch_settings, index_settings = index_settings
      )
      source_start_parameters(source, dat)
    },
    attempt_file = file.path(
      .assessment_cache(assessment_id, translation_tier, database = source),
      "initialization_attempts.rds"
    )
  )
}
