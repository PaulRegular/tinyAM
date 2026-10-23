## Observations ----

years <- 1980:2026
ages <- 3:14
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  sampling_times = c("IS-SMB" = 0.20),
  assumptions = source$assumptions
)

## Background and comparisons ----

age_plus_group <- 14L

background <- c(
  "### Icelandic saithe: 2026 assessment", "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Difference |",
  "|---|---|---|---|",
  "| Years | Catch and age-specific F end in 2025; the spring survey and reported population estimates extend through 2026. | Fit 1980–2026, with missing catch in 2026. | The 2026 year has no catch observation. |",
  "| Ages | Ages 3–14+, with age 14 as the plus group. | Fit ages 3–14 and treat age 14 as the plus group. | The published age structure is retained. |",
  "| N | The report describes recruitment as a moving average of earlier years. | Use IID cohort deviations and freely estimated initial abundance. | tinyAM does not reproduce the reported recruitment model or its initial-state integration. |",
  "| F | The report identifies ages 4–9 as the F range; the benchmark describes an age-correlated random-walk process. | Use correlated random-walk F across ages and years; report Fbar over ages 4–9. | The exact 2026 process covariance and parameter sharing are unavailable. |",
  "| M | Fixed at 0.2 per year for all ages. | Supply fixed M = 0.2 for all ages and years. | This biological assumption is retained. |",
  "| Catch | Commercial catch numbers-at-age, in thousands of fish, through 2025. | Retain the combined catch series with one observation SD. | The full SAM observation covariance is unavailable. |",
  "| Survey | The spring SMB index is reported by age for 1985–2026. | Retain the published index scale and fit age-specific catchability. | The physical scale multiplier is not stated. Survey timing 0.20 is a March midpoint approximation; the SAM fraction is unknown. |",
  "| Biology | Stock weight and maturity are reported by age through 2026. | Use the published annual surfaces in grams per fish and proportions. | No biological surfaces are substituted. |",
  "| Outputs | SAG reports SSB, recruitment, age-4+ reference biomass, and FishingPressure; age-specific N and F are published separately. | Compare the reported N, F and aggregate series on their recorded scales. | Reference biomass is not SSB. The 2026 MFRI Table 2 F statistic is not relabeled as Fbar; SAG FishingPressure ends in 2025. |",
  "",
  "The 2025 benchmark documents a SAM configuration but does not establish every setting used in the 2026 production run. The exact 2026 observation covariance, F process settings, recruitment specification, and survey timing remain unresolved. The MFRI report's Table 2 F differs from its HR series; HR agrees with SAG FishingPressure through 2025, so the undefined F series is not used as Fbar. The fit first uses package starting values; because that default attempt failed numerical checks, one retry uses official N and F outputs as initial guesses only. Inputs and model settings remain unchanged."
)

## Model ----

source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

fit_model <- function(start_par = NULL) {
  tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "cor_rw", mean_ages = 4:9),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_block, sd_form = ~1),
    start_par = start_par,
    silent = silent
  )
}

fit_start <- function() {
  model_dat <- tinyAM::prepare_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "cor_rw", mean_ages = 4:9),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_block, sd_form = ~1)
  )
  start <- tinyAM::make_par(model_dat)
  source_N <- .translation_start_surface(
    source$outputs, "population", "numbers_at_age", years, ages
  )
  source_F <- .translation_start_surface(
    source$outputs, "mortality", "fishing_mortality_at_age", 1980:2025, ages
  )
  start$log_r0 <- log(source_N[1L, 1L])
  start$log_n0 <- log(source_N[1L, -1L])
  start$log_r <- log(source_N[names(start$log_r), 1L])
  start$log_n <- log(source_N[rownames(start$log_n), colnames(start$log_n), drop = FALSE])
  positive_F <- is.finite(source_F) & source_F > 0
  start_log_f <- start$log_f[
    rownames(source_F), colnames(source_F), drop = FALSE
  ]
  start_log_f[positive_F] <- log(source_F[positive_F])
  start$log_f[rownames(source_F), colnames(source_F)] <- start_log_f
  start
}

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = fit_start,
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache",
      source$assessment$assessment_id[[1]],
      "initialization_attempts.rds"
    )
  )
}