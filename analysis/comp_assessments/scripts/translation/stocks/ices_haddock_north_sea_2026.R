## Observations ----
years <- 1972:2026
ages <- 0:8
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)
obs$pop_covariates$N_sd_key <- factor(ifelse(obs$pop_covariates$age == 8, "sd_8+", "sd_1-8"))


## Background and comparisons ----

age_plus_group <- 8

comparison_scales <- c(
  ssb = 0.001, recruitment = 0.001, N = 0.001, abundance = 0.001, biomass = 0.001,
  biomass_at_age = 0.001
)


background <- c(
  "### Northern Shelf haddock: official 2026 SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",

  "| Years | Population and survey series span 1972–2026; commercial catch ends in 2025. | Fit 1972–2026. | Retains the final Q1 survey year; commercial catch is unavailable for 2026. |",

  "| Ages | Ages 0–8+, with age 8 as the plus group. | Retain ages 0–8+, including the age-8 plus group. | The official population age structure is preserved. |",

  "| N | Recruitment follows a random walk. Older-cohort survival has separate process-SD groups for ages 1–7 and the 8+ group; the initial states have no ordinary process density. | Random-walk recruitment, IID cohort-survival deviations, integrated older initial abundance with fixed r0, and age-grouped process SDs through N_sd_key. | The process variance grouping is approximated more closely, but initial-state treatment and the full SAM process specification may differ. |",

  "| F | Random-walk F with AR(1) correlation across ages in annual increments, separate F states by age, and one shared innovation variance. | Correlated F random walks (cor_rw), with age-correlated annual log-F increments and one shared process SD. | The principal temporal and cross-age covariance structure is now represented; exact SAM initial-state and covariance conventions may differ. |",

  "| M | Fixed annual age-specific natural mortality from the official assessment. | Retain the supplied fixed M-at-age surface. | The official numerical M inputs are preserved; no mortality process is estimated. |",

  "| Catch | One total-catch fleet, with commercial catch numbers-at-age through 2025. | Retain source catch-at-age with one estimated lognormal observation SD; missing values are eligible for model-based filling. | tinyAM uses a simplified independent lognormal catch likelihood rather than the complete native observation configuration. |",

  "| Indices | Delta-GAM standardized Q1 ages 1–8+ and Q3+Q4 ages 0–8+, sampled at fractions 0.125 and 0.75. Catchability varies by survey and age; selected young ages have estimated density-dependent q powers. Observation errors are independent. | Retain both survey indices, their sampling times, survey-by-age q parameters, survey-specific observation SDs, and supplied relative SD factors. | The official density-dependent q powers are not represented; observation variance sharing may differ, but both models use independent observation errors. |",

  "| Weights and maturity | Annual stock weights and maturity from the official run; spawning fractions for F and M are zero. | Retain the annual biological surfaces and beginning-of-year spawning convention. | The official biological inputs and SSB timing are preserved. |",

  "",
  "The revised translation uses IID cohort-survival deviations with age-grouped process SDs and correlated F random walks. These settings more closely reproduce the official SAM process structure than the previous IID N model with a single survival SD and stationary AR1 F model.",
  "",
  "Recruitment retains a separate random-walk process. The correlated F process represents AR(1) dependence across ages in annual log-F increments, rather than stationary fluctuations around age-specific mean F.",
  "",
  "The main remaining observation-model difference is the official assessment's density-dependent survey catchability at Q1 age 1 and Q3+Q4 ages 0–1. tinyAM uses ordinary survey-by-age catchability without those power parameters.",
  "",
  "Survey observations retain their native numerical scales because calibrated physical units were not specified. Supplied relative precision factors are retained through equivalent log-scale SD multipliers.",
  "",
  "Comparison units: tinyAM numbers, recruitment and biomass are converted to the official assessment's reported units in the dashboard and comparison tables. Official fitted 2026 values remain distinct from advice forecast quantities."
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
    N_settings = list(process = "iid", init = "latent", sd_form = ~N_sd_key),
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + survey, sd_supplied = ~relative_sd),
    silent = silent
  )
}
