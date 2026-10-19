## Observations ----
years <- 1968:2018
ages <- 0:6
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c(
  "NIGFS-WIBTS-Q1 (G7144)",
  "NIGFS-WIBTS-Q4 (G7655)",
  "UK FSP (B7897)"
)
excluded_catch_years <- c(2003:2005, 2020)

inputs <- source$inputs
input_year <- suppressWarnings(as.integer(as.character(inputs$year)))
excluded_catch <- inputs$type == "catch" &
  inputs$measure == "numbers_at_age" &
  input_year %in% excluded_catch_years
inputs <- inputs[!excluded_catch, , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c(
    "NIGFS-WIBTS-Q1 (G7144)" = 0.125,
    "NIGFS-WIBTS-Q4 (G7655)" = 0.875,
    "UK FSP (B7897)" = 0.125
  ),
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0
obs$pop_covariates$recruitment_block <- factor(
  ifelse(obs$pop_covariates$year < 2002, "pre_2002", "2002_plus"),
  levels = c("pre_2002", "2002_plus")
)
obs$pop_covariates$f_age_group <- factor(
  ifelse(obs$pop_covariates$age >= 3, "age3_plus", paste0("age", obs$pop_covariates$age)),
  levels = c("age0", "age1", "age2", "age3_plus")
)

comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Official aggregate SSB compared with tinyAM beginning-of-year mature biomass at ages 0-6+.",
    reason = "The report does not specify the official SSB timing; tinyAM uses beginning-of-year biomass."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Official and tinyAM recruitment at age 0 over 1968-2018.",
    reason = "The source point-estimate convention is not specified, and tinyAM uses a random-walk recruitment process rather than the source Beverton-Holt blocks."
  )
)

## Background and comparisons ----
age_plus_group <- 6

background <- c(
  "### Irish Sea cod: official 2026 SS3 assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | The assessment spans 1968-2025; catch years 2003-2005 and 2020 are excluded. | Fit 1968-2018 and omit the source-excluded catch years within this period. | The weight surface is complete through 2018. The 2019 age-0 weight is missing, and 2020 weights were not included by the source assessment. |",
  "| Ages | Ages 0-6+, with recruitment at age 0. | Fit ages 0-6+, retaining the age-6 plus group. | The official age range is retained for the supported years. |",
  "| N and recruitment | Stock Synthesis uses a Beverton-Holt relationship with recruitment blocks. | Use age-0 recruitment as a random walk and IID cohort-survival deviations. | Recruitment age is aligned, but the source curve, blocks and initial-state treatment are not reproduced. |",
  "| F | SS3 selectivity follows age-based empirical random walks; older ages follow age 3. | Use an AR1 temporal process with separate mean levels for ages 0, 1 and 2 and a shared age-3+ mean. | The source's selectivity process and exact current settings are not reproduced. |",
  "| M | Fixed age-specific natural mortality: 1.7 at age 0, 0.714 at age 1, 0.68 at age 2, and 0.65 at ages 3+. | Use the published fixed M values. | The source assumption is retained. |",
  "| Catch and indices | Catch-at-age is in thousands of fish. Indices include Q1 ages 1-4, Q4 age 0, and UK FSP Q1 ages 2-6+. | Retain catch-at-age and all three age-specific indices. | The source index scales and current likelihood weights are not fully documented. Observation SDs are estimated by survey. |",
  "| Biology | Landings mean weight is used for the stock and catch; maturity is annual. | Retain source weights and maturity for ages 0-6+, with zero spawning-time survival to represent beginning-of-year mature biomass. | Official SSB timing is not stated, so the SSB comparison is approximate. |",
  "| Available comparisons | The report provides aggregate SSB, age-0 recruitment and Fbar ages 2-4, but no age-specific N/F outputs. | Compare SSB and age-0 recruitment over 1968-2018. | Recruitment is approximate because the source point-estimate convention is undocumented. N/F at age are unavailable; tinyAM's Fbar is abundance-weighted and cannot be reconstructed from the reported Fbar alone. |",
  "",
  "This fit preserves the source recruitment age and age-specific observations while avoiding invented values. It does not cover the assessment's 2019-2025 terminal period and does not reproduce the full SS3 likelihood."
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
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(
      process = "iid",
      mu_form = ~ bh(ssb, lag = 1) + recruitment_block
    ),
    F_settings = list(
      process = "ar1",
      mu_form = ~ 0 + f_age_group,
      mean_ages = 2:4
    ),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + survey),
    silent = silent
  )
}
