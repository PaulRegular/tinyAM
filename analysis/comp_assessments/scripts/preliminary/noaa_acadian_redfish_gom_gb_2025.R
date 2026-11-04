## Observations ----
years <- 1969:2024
ages <- 1:26
inputs <- source$inputs
stock <- inputs$type == "weight" & inputs$measure == "start_weight_at_age"
inputs$measure[stock] <- "weight_at_age"
inputs <- inputs[!(inputs$type == "weight" & inputs$measure == "spawning_weight_at_age"), ]
catch_weight <- inputs$type == "catch_weight"
inputs$measure[catch_weight] <- "weight_at_age"
weights <- .translation_repeat_biology(inputs[catch_weight, ], years, NULL, "Catch weights")
inputs <- rbind(inputs[!catch_weight, ], weights)

# Preserve the unknown total-index multiplier by reconstructing native relative indices.
total <- inputs[inputs$type == "index" & inputs$measure == "total_numbers", ]
composition <- inputs[inputs$type == "index" & inputs$measure == "proportion_at_age", ]
key <- function(x) paste(x$survey, x$year)
if (anyDuplicated(key(total))) cli::cli_abort("Redfish indices need one source total per survey and year.")
matched <- match(key(composition), key(total))
if (anyNA(matched)) cli::cli_abort("Redfish survey compositions need matching source totals.")
index <- composition
index$value <- composition$value * total$value[matched]
index$measure <- "numbers_at_age"
index$unit <- "native survey index"
index$transformation <- "Original number-index total multiplied by its number proportion at age; no physical multiplier inferred."
inputs <- rbind(inputs[inputs$type != "index", ], index)
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, assumptions = source$assumptions)
obs$maturity$spawn_time <- 0.4
obs$catch$F_age <- tinyAM::cut_ages(obs$catch$age, c(1:9, 26))
q_age <- ifelse(grepl("fall", obs$index$survey, ignore.case = TRUE),
  pmin(obs$index$age, 4), pmin(obs$index$age, 7))
obs$index$q_group <- interaction(obs$index$survey, q_age, drop = TRUE)

N_settings <- list(process = "off", init = "exp")
rec_settings <- list(process = "ar1", mu_form = ~1)
F_settings <- list(process = "off", mu_form = ~F_age + factor(year))
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_group, sd_form = ~0 + survey)
age_plus_group <- 26L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "Source N and F with shared January stock weights and maturity, surviving to fraction 0.4.",
  reason = "The official spawning-weight schedule differs from January weights; native SSB remains separate."))

## Background ----
background <- c(
  "### Acadian redfish: preliminary comparison with the Official 2025 WHAM assessment", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1963–2024; ages 1–26+. Catch age compositions begin in 1969 and are absent in many later years. | Fit 1969–2024 at the same ages. | Earlier aggregate catches cannot supply missing age shares. Later unsampled catch-age years remain missing; no age shares are interpolated. |",
  "| N | Deterministic cohort survival, AR1 recruitment about a mean, equilibrium initial abundance. | Deterministic survival and AR1 recruitment about an estimated mean; exponential initial abundance. | The source estimates a separate equilibrium F for initial ages; tinyAM uses first-year mortality. |",
  "| F | A single annual F level and static selectivity; ages 9+ fixed fully selected. | Separable fixed year effects and age groups 1–8 and 9+. | Relative age coefficients are not constrained to be below the terminal level. No additional F process is imposed. |",
  "| M | Fixed at 0.05 per year for all ages. | Same supplied M. | Fixed-M agreement is by construction. |",
  "| Catch | Biomass totals plus logistic-normal AR1 age compositions. | Reconstruct numbers using matched totals, number shares and the original static catch weights in composition years. | The aggregate/composition dependence and likelihood weights are not reproduced; totals from years without compositions do not become invented catch-at-age. |",
  "| Index | Spring and fall number indices, separate number compositions, native multipliers unspecified. | Multiply each total by its own age proportions; retain native scale and timing 3.5/12 and 9.5/12. | Derived ages share a total. Independent age errors with survey SDs approximate the native aggregate/composition likelihood; aggregate CVs are not treated as per-age SDs. |",
  "| q | Fall selectivity fixed at ages 4+; spring at ages 7+. | Retain those q plateaus with separate younger-age groups. | Younger relative q coefficients remain unconstrained. |",
  "| Biology and SSB | Static catch, spawning and January weight schedules; SSB at fraction 0.4. | Use January weights and maturity with spawning fraction 0.4. | Common-definition SSB is approximate because tinyAM uses one weight surface; native source SSB uses separate spawning weights. |",
  "", "Source observations and outputs are never substituted for missing catch compositions. Default starts are attempted first; a numerical failure permits one aligned N/F-state retry, leaving unavailable parameter starts at package defaults."
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
