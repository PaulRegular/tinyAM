## Observations ----
years <- 1968:2025
ages <- 1:6
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("NIGFS-WIBTS-Q1 (G7144)", "UK FSP (B7897)")
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
    "UK FSP (B7897)" = 0.125
  ),
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0
obs$pop_covariates$f_age_group <- factor(
  ifelse(obs$pop_covariates$age >= 3, "age3_plus", paste0("age", obs$pop_covariates$age)),
  levels = c("age1", "age2", "age3_plus")
)

comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Official aggregate SSB compared with tinyAM beginning-of-year mature biomass at ages 1-6+.",
    reason = "The report does not specify the official SSB timing; tinyAM uses beginning-of-year biomass."
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
  "| Years | The model uses 1968-2025 data, with some catch years excluded. | Fit 1968-2025 and omit catch observations in 2003-2005 and 2020. | The report identifies these years as excluded because age sampling was limited. |",
  "| Ages | Ages 0-6+, with recruitment at age 0. | Fit ages 1-6+, retaining the age-6 plus group. | The 2019 age-0 weight is missing. Age-0 catch and survey data are omitted instead of filling that value. Recruitment is therefore not directly comparable. |",
  "| N and recruitment | Stock Synthesis uses a Beverton-Holt relationship with recruitment blocks. | Use age-1 recruitment as a random walk and IID cohort-survival deviations. | The age-0 stock-recruit curve and its blocks are not represented. |",
  "| F | SS3 selectivity follows age-based empirical random walks; older ages follow age 3. | Use an AR1 temporal process with separate mean levels for ages 1 and 2 and a shared age-3+ mean. | The source's selectivity process and exact current settings are not reproduced. |",
  "| M | Fixed age-specific natural mortality: 0.714 at age 1, 0.68 at age 2, and 0.65 at ages 3+. | Use the published fixed M values. | Age 0 is outside the selected tinyAM age range. |",
  "| Catch and indices | Catch-at-age is in thousands of fish. Indices include Q1 ages 1-4, Q4 age 0, and UK FSP Q1 ages 2-6+. | Retain catch-at-age and the Q1 and FSP indices for selected ages; omit the age-0-only Q4 index. | The source index scales and current likelihood weights are not fully documented. Observation SDs are estimated by survey. |",
  "| Biology | Landings mean weight is used for the stock and catch; maturity is annual. | Retain source weights and maturity for ages 1-6+, with zero spawning-time survival to represent beginning-of-year mature biomass. | Official SSB timing is not stated, so the SSB comparison is approximate. |",
  "| Available comparisons | The report provides aggregate SSB, age-0 recruitment and Fbar ages 2-4, but no age-specific N/F outputs. | Compare aggregate SSB approximately. | Recruitment age differs; N/F at age are unavailable; tinyAM's Fbar is abundance-weighted and cannot be reconstructed from the reported Fbar alone. |",
  "",
  "The fit is a limited age-1+ approximation intended to test the recoverable observations. It does not reproduce the official age-0 recruitment process or the full SS3 likelihood."
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
    rec_settings = list(process = "rw"),
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
