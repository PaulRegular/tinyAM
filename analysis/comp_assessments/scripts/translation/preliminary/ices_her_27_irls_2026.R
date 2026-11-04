## Observations ----
years <- 2002:2024
ages <- 1:9
survey <- "Celtic Sea Pelagic Ecosystem Survey (CSHAS)"
inputs <- source$inputs
annual <- function(type, measure) {
  x <- inputs[inputs$type == type & inputs$measure == measure &
    inputs$year %in% years & inputs$age %in% ages, ]
  if (any(x$year_basis != "model_year") || anyDuplicated(paste(x$year, x$age))) {
    cli::cli_abort("Celtic herring requires unique April–March model-year age rows.")
  }
  data.frame(year = x$year, age = x$age, obs = x$value)
}
# Integer labels retain the source April–March origin; no calendar-year shift.
catch <- annual("catch", "numbers_at_age")
catch$unit <- "thousand fish"
index <- annual("index", "numbers_at_age")
index <- index[!index$year %in% c(2004, 2017, 2024), ]
index$survey <- survey
index$samp_time <- 0.5
index$unit <- "thousand fish"
index$q_group <- factor(ifelse(index$age %in% 3:5, "ages3_5", as.character(index$age)))
weight <- annual("weight", "weight_at_age")
weight$unit <- "kg"
grid <- expand.grid(year = years, age = ages)
fixed <- function(type, measure) {
  x <- inputs[inputs$type == type & inputs$measure == measure & is.na(inputs$year), ]
  if (anyDuplicated(x$age) || !setequal(x$age, ages)) {
    cli::cli_abort("Celtic herring requires the complete published static age vector.")
  }
  x$value[match(grid$age, x$age)]
}
maturity <- transform(grid, obs = fixed("maturity", "maturity_at_age"))
pop_covariates <- transform(grid, M_assumption = fixed("M", "natural_mortality_at_age"))
catch$F_age <- factor(pmin(catch$age, 3))
obs <- list(catch = catch, index = index, weight = weight, maturity = maturity,
  pop_covariates = pop_covariates)
attr(obs, "translation") <- list(year_basis = "model_year", year_origin = "1 April",
  notes = "Source numbers are already thousand fish and weights kg/fish. Published values retained, including the unresolved 2019 age-2 index conflict.")

N_settings <- list(process = "off", init = "free")
rec_settings <- list(process = "iid", mu_form = ~1)
F_settings <- list(process = "off", mu_form = ~F_age + factor(year), mean_ages = 2:5)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_group, sd_form = ~1)
age_plus_group <- 9L
comparison_aggregates <- "ssb"
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "Reported native SSB compared with April-origin mature biomass using published spawning weights.",
  reason = "Exact spawning survival fractions are unavailable; one weight surface cannot represent distinct within-year biological schedules."))

## Background ----
background <- c(
  "### Celtic Sea herring: preliminary comparison with the Official 2026 ASAP assessment", "",
  print_sources(source$assessment), "",
  "All annual labels refer to 1 April–31 March, labelled by the first year; they are not calendar-year states. The October acoustic index is at fraction 0.5 of this annual step.", "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Coverage | 1958–2025, winter-ring ages 1–9+. | 2002–2024, the acoustic-era overlap with complete weights. | Missing 2025 age-9+ biology is not filled; earlier catch-only years are not fitted. |",
  "| N and recruitment | ASAP catch-at-age dynamics; no fitted stock–recruitment relationship. | Deterministic survival, IID recruitment around a mean and free initial abundance. | Recruitment penalties and initial-state treatment are approximate. |",
  "| F and q | Static fishery selectivity, fully selected at age 3; survey ages 3–5 share selection. | Annual F effects with ages 3+ shared; one q for ages 3–5 and separate q at ages 2, 6 and 7. | Relative effects are not bounded below the fully selected level. |",
  "| M and biology | Fixed published age-specific M, annual spawning weights and static maturity. | Retain those values. | The model uses one weight surface; beginning-step mature biomass is only an approximate native SSB comparison. |",
  "| Observations | Catch biomass and number compositions; acoustic age inputs with total CV 0.5. | Fit published catch and index counts with estimated lognormal SDs. | Shared totals, composition weights and aggregate CVs are not equivalent to independent age errors. |",
  "| Survey exclusions | 2004, 2017 and 2024 excluded from tuning. | Same exclusions. | No unavailable or excluded zero index is fitted. |",
  "", "The published 2019 age-2 acoustic value exceeds its published total. It is retained without correction and limits interpretation. The conflicting 2025 catch totals are outside this fit. No agreement-based changes or source-state starts are used."
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
      cli::cli_abort("Aligned official N/F states are unavailable; no source-state retry is possible.")
    }, attempt_file = file.path(.assessment_cache(
      source$assessment$assessment_id, translation_tier, database = source), "initialization_attempts.rds"))
}
