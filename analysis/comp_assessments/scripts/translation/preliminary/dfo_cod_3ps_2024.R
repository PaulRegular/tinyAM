## Observations ----
years <- 1959:2024
ages <- 2:14
surveys <- c(
  "DFO RV Spring", "ERHAPS", "GEAC",
  "Sentinel Gillnet", "Sentinel Line-trawl"
)
sampling_times <- c(
  "DFO RV Spring" = 0.25,
  "ERHAPS" = 0.17,
  "GEAC" = 0.875,
  "Sentinel Gillnet" = 0.5,
  "Sentinel Line-trawl" = 0.5
)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  sampling_times = sampling_times,
  surveys = surveys,
  assumptions = source$assumptions
)

condition <- source$inputs[
  source$inputs$type == "covariate" &
    source$inputs$measure == "condition_index_for_mortality",
  c("year", "value"), drop = FALSE
]
condition$year <- as.integer(condition$year)
if (!identical(sort(condition$year), years) || anyDuplicated(condition$year)) {
  stop("The annual fish-condition index must contain one value for every modeled year.")
}
obs$pop_covariates$M_condition <- condition$value[
  match(obs$pop_covariates$year, condition$year)
]
obs$pop_covariates$M_condition_age_2_5 <- ifelse(
  obs$pop_covariates$age <= 5,
  obs$pop_covariates$M_condition,
  0
)
obs$pop_covariates$M_condition_age_6_plus <- ifelse(
  obs$pop_covariates$age >= 6,
  obs$pop_covariates$M_condition,
  0
)

comparison_outputs <- source$outputs[
  source$outputs$year %in% 2020:2024, , drop = FALSE
]
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Reported annual SSB, ages 2-14+",
    reason = "The modeled ages and original maturity and stock weights are retained, but population dynamics and spawning calculation details are simplified."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-2 abundance in the same calendar year",
    reason = "Both estimates refer to age 2; recruitment and state-process assumptions differ."
  )
)
age_plus_group <- 14

## Background and model ----
background <- c(
  "### 3Ps Atlantic cod: 2024 DFO HYBRID assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1959-2024, ages 2-14+; catch-at-age ends in 2023. | 1959-2024, ages 2-14+. | Retain the full population and biology span; the final year has no catch-at-age observation. |",
  "| Recruitment and N | Age-2 recruitment follows a random walk; abundance and survival have additional state-space structure. | Age-2 recruitment random walk, IID cohort-survival deviations, exponential initial abundance. | The source survival process is simplified to independent cohort deviations; initial-state integration is omitted. |",
  "| F | Age-correlated random walk for ages 3+, separate age-2 F, with a break at the 1993 moratorium. | Correlated random walk across all modeled ages with one uninterrupted history. | The age-2 separation and moratorium break are not supported by the compact tinyAM F structure used here. |",
  "| M | Time-varying M is informed by the annual fish-condition index, with condition-response coefficients shared across ages 2–5 and 6+. | Estimate age-specific log-M levels and separate log-M responses to condition for ages 2–5 and 6+; no additional M process. | The two age-group slopes approximate the source condition-based M equation, which was not fully recovered. |",
  "| Catch | Censored annual landings with bounds and continuation-ratio catch-at-age composition. | Direct reported catch-at-age observations with lognormal error. | The censored landings and composition likelihoods are not reproduced; source catch-age values are not altered. |",
  "| Indices and q | Five index series; RV q sharing includes offshore/inshore adjustments; zero observations are censored. | Five direct age-specific index series, with separate q by age and survey. Zero values are treated as missing by the log observation likelihood. | Survey units are unresolved, so each survey has its own q. The source q sharing, spatial adjustment and zero-censoring are omitted. |",
  "| Timing and biology | Timing varies among years and surveys; stock weights, catch weights and maturity are reported by year and age. | Original annual stock weights and maturity; RV April, ERHAPS late winter, GEAC late fall, and sentinel midyear timing approximations. | The exact model timing fractions are unavailable; seasonal midpoints are analysis choices, not source inputs. |",
  "| Comparison outputs | Annual aggregate SSB, age-2 recruitment, and mean F/M for ages 5-8 are reported for 2020-2024. | Compare recent SSB and age-2 recruitment where definitions overlap. | Numerical age-specific N/F/M and older annual output trajectories are unavailable; grouped F/M cannot be reconciled to tinyAM's abundance-weighted definitions. |",
  "",
  "The fit uses package defaults first. If numerical checks fail, one retry starts age-2 recruitment and ages 5-8 F from the published 2020-2024 summaries; these values initialize the fit and do not constrain estimates. Other parameters retain package defaults."
)

fit <- NULL
N_settings <- list(process = "iid", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "cor_rw", mu_form = NULL, mean_ages = 5:8)
M_settings <- list(
  process = "off",
  mu_form = ~ 0 + factor(age) + M_condition_age_2_5 + M_condition_age_6_plus,
  mu_supplied = NULL,
  mean_ages = 5:8
)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~ 0 + q_key, sd_form = ~ 0 + survey)

source_start_par <- function(starts, outputs) {
  start_years <- 2020:2024
  recruitment <- outputs[
    outputs$measure == "recruitment" &
      as.integer(outputs$year) %in% start_years &
      as.integer(outputs$age) == 2L,
    , drop = FALSE
  ]
  fbar <- outputs[
    outputs$measure == "Fbar" & as.integer(outputs$year) %in% start_years,
    , drop = FALSE
  ]
  if (!setequal(as.integer(recruitment$year), start_years) ||
      anyDuplicated(recruitment$year) ||
      !all(recruitment$unit == "million fish") ||
      any(!is.finite(recruitment$value) | recruitment$value <= 0)) {
    cli::cli_abort("Published 3Ps age-2 recruitment does not cover 2020-2024 in million fish.")
  }
  if (!setequal(as.integer(fbar$year), start_years) ||
      anyDuplicated(fbar$year) ||
      !all(fbar$unit == "per year") ||
      any(!is.finite(fbar$value) | fbar$value <= 0)) {
    cli::cli_abort("Published 3Ps Fbar does not cover 2020-2024 in per-year units.")
  }

  recruitment <- recruitment[match(start_years, as.integer(recruitment$year)), ]
  fbar <- fbar[match(start_years, as.integer(fbar$year)), ]
  years_in_fit <- as.character(start_years)
  if (!all(years_in_fit %in% names(starts$log_r)) ||
      !all(years_in_fit %in% rownames(starts$log_f)) ||
      !all(as.character(5:8) %in% colnames(starts$log_f))) {
    cli::cli_abort("Published 3Ps starting values do not align with the fitted years and ages.")
  }

  starts$log_r[years_in_fit] <- log(recruitment$value * 1000)
  starts$log_f[years_in_fit, as.character(5:8)] <- matrix(
    log(fbar$value), nrow = length(start_years), ncol = 4L
  )
  starts
}

if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_model <- function(start_par = NULL) {
    tinyAM::fit_tam(
      data = obs,
      years = years,
      ages = ages,
      N_settings = N_settings,
      rec_settings = rec_settings,
      F_settings = F_settings,
      M_settings = M_settings,
      catch_settings = catch_settings,
      index_settings = index_settings,
      silent = silent,
      start_par = start_par
    )
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      model_dat <- tinyAM::prepare_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = N_settings,
        rec_settings = rec_settings,
        F_settings = F_settings,
        M_settings = M_settings,
        catch_settings = catch_settings,
        index_settings = index_settings
      )
      source_start_par(tinyAM::make_par(model_dat), source$outputs)
    },
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      source$assessment$assessment_id[[1]], "initialization_attempts.rds"
    )
  )
}
