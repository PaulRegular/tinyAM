source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
assessment_id <- source$assessment$assessment_id[[1L]]
years <- 1980:2019
ages <- 1:11
surveys <- c(
  "NEFSC spring Albatross",
  "NEFSC spring Bigelow",
  "NEFSC fall Albatross",
  "NEFSC fall Bigelow"
)

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  index_weight_source = "source",
  surveys = surveys,
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0.25
obs$index$q_key <- interaction(obs$index$survey, obs$index$age, drop = TRUE)

## Background and comparisons ----
age_plus_group <- 11
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "ssb")
comparison_age_groups <- list(F = list(`6-9` = 6:9))
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total January 1 abundance summed over ages 1-11+.",
    reason = "Both series use the same model ages; the source total is the sum of fitted marginal N-at-age estimates."
  ),
  recruitment = list(
    status = "matched",
    definition = "January 1 abundance at age 1 in the same calendar year.",
    reason = "The source model represents recruitment at age 1."
  ),
  F_bar = list(
    status = "matched",
    definition = "N-weighted mean F over ages 6-9, calculated from source and tinyAM F-at-age.",
    reason = "The source-reported Fbar is the arithmetic mean over ages 6-9; the dashboard compares a common N-weighted definition instead and retains native Fbar separately."
  ),
  M = list(
    status = "approximate",
    definition = "Fixed M input at age, checked against the source value of 0.3 per year.",
    reason = "The same source value is supplied to tinyAM; this checks input mapping and is not an estimated fit comparison."
  ),
  M_bar = list(
    status = "approximate",
    definition = "N-weighted mean of the fixed M input over ages 1-11+.",
    reason = "The same source value is supplied to tinyAM; this checks input mapping and is not an estimated fit comparison."
  ),
  ssb = list(
    status = "matched",
    definition = "Spawning biomass at 0.25 of the year using the source maturity and stock-weight inputs.",
    reason = "The source and tinyAM use the same annual timing and biological inputs."
  )
)

background <- c(
  "### Gulf of Maine/Georges Bank American plaice: 2022 WHAM Run29F4 candidate",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Run29F4 covers 1980-2019, ages 1-10 and 11+. | Fit 1980-2019, ages 1-10 and 11+. | The candidate run's complete historical grid and age-11 plus group are retained. |",
  "| N and recruitment | WHAM estimates recruitment and all age-year abundance states, with independent process errors and separate SDs for age 1 and ages 2+. | Use IID recruitment and IID cohort-survival deviations with free initial older ages. | This approximates independent state errors but does not reproduce WHAM's joint random-effects structure or initial-state integration. |",
  "| F | One fishery uses logistic selectivity with time-varying selectivity random effects. | Use an age-specific logistic mean F with an AR1 process over years and ages. | The process is a compact approximation; it does not reproduce the source selectivity random-effect structure. |",
  "| M | Fixed at 0.3 per year for all ages and years. | Use the same supplied M schedule with no M process. | Fixed M agreement checks the input translation only. |",
  "| Catch | Aggregate catch numbers are fitted with a lognormal likelihood, alongside number-at-age compositions with a logistic-normal likelihood and effective sample size 100. | Reconstruct catch-at-age by multiplying annual total catch numbers by source number proportions, then fit one lognormal age-specific observation model. | Separate aggregate and composition likelihoods, composition weights, and source catch error are not reproduced. |",
  "| Surveys and q | Four vessel-season series have aggregate biomass indices, number compositions, year-specific log SDs, and effective sample size 50. Bigelow indices are not calibrated to Albatross units. | Reconstruct age-specific abundance indices from each aggregate index, its number composition, and source survey weights; estimate q by survey and age and use source log SDs. | tinyAM does not fit the aggregate and composition likelihoods separately, propagate composition uncertainty, or reproduce selectivity constraints. |",
  "| Survey timing | The input array records 0.25, 0.75, 0.25, and 0.75 for the spring Albatross, spring Bigelow, fall Albatross, and fall Bigelow series. The middle two values conflict with the season labels. | Retain the native timing values. | The mismatch remains unresolved; no timing values are silently corrected. |",
  "| Biology and SSB | Fixed age-specific maturity and annual stock weights; SSB is calculated at 0.25 of the year. | Use the source age-specific maturity and stock weights with spawn time 0.25. | The reported definition and input surfaces are retained. |",
  "| Outputs | The cached candidate fit provides N-at-age, F-at-age, SSB, age-1 recruitment, total abundance, and Fbar ages 6-9 with marginal uncertainty. | Compare N-at-age, F-at-age, recruitment, total abundance, SSB, and N-weighted F over ages 6-9. | Total abundance and recruitment are derived from marginal N-at-age; their uncertainty is unavailable without cross-age covariance. The source arithmetic Fbar is retained separately from the common N-weighted comparison. |",
  "",
  "Run29F4 is one of three 2022 research-track candidates, and the exact candidate selected for management advice is not established. The 2024 Management Track report provides later summary outputs but no age-specific fitted surfaces in the recovered archive; those values remain separate in the database. The 2026 Management Track report is provisional and is not treated as official production advice. This dashboard compares tinyAM only with the 2022 Run29F4 candidate, not with the later summary tables."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = list(process = "iid", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ 1),
  F_settings = list(process = "ar1", mu_form = ~ logistic(age), mean_ages = 6:9),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(
    q_form = ~ 0 + q_key,
    sd_form = ~ 0,
    sd_supplied = ~ relative_sd
  ),
  silent = silent,
  ...
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(
      comparison_outputs, "population", "numbers_at_age", years, ages, 1000
    )
    f <- .translation_start_surface(
      comparison_outputs, "mortality", "fishing_mortality_at_age", years, ages
    )
    list(
      log_r0 = log(n[1L, 1L]),
      log_r = log(n[-1L, 1L]),
      log_n0 = log(n[1L, -1L]),
      log_f = log(f)
    )
  }, attempt_file = file.path(
    "analysis/comp_assessments/results/cache", assessment_id,
    "initialization_attempts.rds"
  ))
}
