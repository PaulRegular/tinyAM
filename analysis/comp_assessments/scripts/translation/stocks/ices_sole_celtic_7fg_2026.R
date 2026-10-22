source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1971:2025
ages <- 1:10
assessment_id <- source$assessment$assessment_id[[1L]]

inputs <- source$inputs
inputs$value <- suppressWarnings(as.numeric(as.character(inputs$value)))
excluded_catch <- inputs$type == "catch" & inputs$measure == "numbers_at_age" &
  inputs$age %in% c("1", "2") & inputs$year < 2004
excluded_index <- inputs$type == "index" & inputs$measure == "numbers_at_age" &
  !(inputs$age %in% as.character(1:5))
inputs <- inputs[!excluded_catch & !excluded_index, , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = "UK(E&W)-BTS-Q3",
  sampling_times = c("UK(E&W)-BTS-Q3" = 0.625),
  assumptions = source$assumptions
)
obs$index$q_key <- factor(as.character(obs$index$age), levels = as.character(ages))

age_plus_group <- 10L
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Sum of official N-at-age over ages 1-10+; no uncertainty is available for this derived total.",
    reason = "The source N-at-age table covers the full modeled age range; total abundance is derived from those reported values."
  ),
  biomass = list(
    status = "approximate",
    definition = "Total biomass across ages 1-10+ using the reported stock weights.",
    reason = "The report does not make the timing and full total-biomass calculation explicit."
  ),
  ssb = list(
    status = "approximate",
    definition = "Mature biomass across ages 1-10+ using the reported weights and maturity.",
    reason = "The source's exact spawning-time convention is not established; tinyAM uses its start-of-year default."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 abundance in the same calendar year, in thousands of fish.",
    reason = "Age and units match; the process and initial-state treatments differ."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Source Fbar ages 3-8 compared with tinyAM's abundance-weighted F over ages 3-8.",
    reason = "The source report does not establish the same weighting convention as tinyAM."
  ),
  M = list(
    status = "approximate",
    definition = "Fixed natural mortality input of 0.1 per year at all modeled ages.",
    reason = "This checks the supplied assumption; it is not a comparison of independently estimated M."
  )
)

background <- c(
  "### Celtic Sea common sole: 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1971-2025, ages 1-10+, with age 10 as the plus group. | Fit the same years and ages. | The source grids are retained. |",
  "| N and recruitment | SAM uses a random walk on log age-1 abundance; N-at-age is reported. | Exponential initial abundance, an IID cohort process, and a separate random walk for recruitment. | tinyAM simplifies SAM's age-specific process and initial-state integration. |",
  "| F and M | AR(1) F process; Fbar ages 3-8; M is fixed at 0.1 for all ages and years. | AR(1) F with mean ages 3-8; the reported M value is supplied as fixed. | Other SAM process and likelihood details are simplified; M is an input check, not an estimated comparison. |",
  "| Catch | Catch numbers-at-age and catch weights-at-age; the SAM configuration omits ages 1 and 2 before 2004. | Keep source scales and omit those same catch cells only in the analysis observations. | Canonical source values remain unchanged. |",
  "| Indices and q | Five commercial beam-trawl biomass indices and the UK(E&W)-BTS-Q3 age-specific index at ages 1-5. | Fit the age-specific BTS-Q3 series with age-specific q; retain commercial indices in the database only. | tinyAM's current direct index likelihood cannot represent the aggregate biomass series faithfully. |",
  "| Biology and SSB | Annual stock weights, maturity by age, and fixed M. | Use the reported stock weights, maturity and M. | The report's exact SSB timing is not clear; tinyAM uses a start-of-year default. |",
  "| Outputs | N-at-age and F-at-age are reported; recruitment, total biomass, SSB, catch and Fbar include low/value/high values. | Compare N/F and available aggregate series. | The low/high level is not identified in the table; age-specific uncertainty is unavailable. |",
  "",
  "The Q3 timing of 0.625 is a midpoint assumption for tinyAM, not a reported source-model setting. The five commercial biomass series are retained on their native scales because the report does not state physical units. The fit uses package-default starts first; one retry may use aligned source N/F values only if numerical checks fail."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()

  fit_model <- function(start_par = NULL) {
    tinyAM::fit_tam(
      data = obs,
      years = years,
      ages = ages,
      N_settings = list(process = "iid", init = "exp"),
      rec_settings = list(process = "rw", mu_form = NULL),
      F_settings = list(process = "ar1", mu_form = ~factor(age), mean_ages = 3:8),
      M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
      catch_settings = list(sd_form = ~1),
      index_settings = list(q_form = ~0 + q_key, sd_form = ~1),
      silent = silent,
      start_par = start_par
    )
  }

  retry_start <- function() {
    source_N <- .translation_start_surface(
      source$outputs, "population", "numbers_at_age", years, ages
    )
    source_F <- .translation_start_surface(
      source$outputs, "mortality", "fishing_mortality_at_age", years, ages
    )
    dat <- tinyAM::prepare_tam(
      data = obs, years = years, ages = ages,
      N_settings = list(process = "iid", init = "exp"),
      rec_settings = list(process = "rw", mu_form = NULL),
      F_settings = list(process = "ar1", mu_form = ~factor(age), mean_ages = 3:8),
      M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
      catch_settings = list(sd_form = ~1),
      index_settings = list(q_form = ~0 + q_key, sd_form = ~1)
    )
    start_par <- tinyAM::make_par(dat)
    start_par$log_r0 <- log(source_N[1L, 1L])
    start_par$log_r <- log(source_N[-1L, 1L])
    start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
    start_par$log_f <- log(source_F)
    start_par
  }

  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      assessment_id, "initialization_attempts.rds"
    )
  )
}
