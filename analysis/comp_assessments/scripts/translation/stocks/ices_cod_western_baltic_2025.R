## Observations ----
years <- 1985:2025
ages <- 0:7
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("BITS Q4", "BITS Q1", "FEJUCS")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("BITS Q4" = 0.8477816, "BITS Q1" = 0.1653494, "FEJUCS" = 0.85),
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0
obs$pop_covariates$f_age_group <- factor(ifelse(
  obs$pop_covariates$age == 0, "age0",
  ifelse(obs$pop_covariates$age == 1, "age1",
    ifelse(obs$pop_covariates$age == 2, "age2",
      ifelse(obs$pop_covariates$age <= 4, "age3_4", "age5plus")
    )
  )
), levels = c("age0", "age1", "age2", "age3_4", "age5plus"))
obs$catch$sd_block <- factor(ifelse(obs$catch$age == 1, "age1", "age2plus"),
  levels = c("age1", "age2plus"))
obs$index$q_key <- interaction(obs$index$survey, obs$index$age,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- with(obs$index, ifelse(
  survey == "BITS Q4",
  ifelse(age == 0, "q4_age0", ifelse(age == 1, "q4_age1", "q4_age2plus")),
  ifelse(survey == "BITS Q1",
    ifelse(age == 1, "q1_age1", ifelse(age == 2, "q1_age2", "q1_age3plus")),
    "FEJUCS_age0"
  )
))
obs$index$sd_block <- factor(obs$index$sd_block, levels = c(
  "q4_age0", "q4_age1", "q4_age2plus", "q1_age1", "q1_age2", "q1_age3plus", "FEJUCS_age0"
))

## Comparisons ----
age_plus_group <- 7L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total N over ages 0-7+, in thousand fish.",
    reason = "The source total is derived from the same complete SAM N-at-age grid; uncertainty is unavailable for the derived total."
  ),
  ssb = list(
    status = "matched",
    definition = "Start-of-year SSB using N-at-age, stock weights and maturity at ages 0-7+.",
    reason = "The source and tinyAM use the same supplied stock weights, maturity and zero pre-spawning mortality fractions."
  ),
  recruitment = list(
    status = "approximate",
    definition = "tinyAM recruitment at age 0 compared with same-year source N-at-age 0.",
    reason = "The SAM report calls age 1 recruitment; tinyAM's recruitment is the youngest modeled age (0), so the dashboard compares same-year N at age 0 rather than shifting cohorts."
  ),
  F = list(
    status = "approximate",
    definition = "Annual F-at-age over the shared age range.",
    reason = "SAM has no F state at age 0 and shares one state across ages 5-7; tinyAM estimates age 0 and does not share realized states across ages 5-7."
  ),
  M = list(
    status = "matched",
    definition = "Fixed annual natural mortality at ages 0-7+.",
    reason = "The same fixed source M-at-age surface is supplied to tinyAM."
  )
)

background <- c(
  "### Western Baltic cod: official 2025 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM covers 1985-2025 at ages 0-7; catch ends in 2024, while Q1 survey data extend to 2025. | Fit the full 1985-2025 and age 0-7 grid. | The source's unreported 2025 age-0 N cell is excluded from the comparison reference. |",
  "| Catch | Combined landings, discards and recreational catch-at-age in thousand fish, ages 1-7+, with native observation weights of 1 through 2019 and 0.1 in 2020-2024. | Retain observed catch-at-age and estimate two catch SD groups, age 1 and ages 2-7. | tinyAM does not apply SAM's observation weights; their likelihood-scale meaning is not specified in the report. |",
  "| Surveys | BITS Q4 ages 0-4+, BITS Q1 ages 1-4+, and FEJUCS age 0; fitted sample times are 0.8477816, 0.1653494 and 0.85. | Retain all age-specific indices, native scales and timings, with survey-by-age q and SAM-like SD groups. | The BITS Q1 cross-age residual correlation is not represented. |",
  "| N and recruitment | SAM reports age-1 recruitment, with separate age-0 and shared age-1+ N-process variances. | Use age-0 recruitment with its own IID process SD, an exponential initial state and IID older-cohort deviations. | tinyAM's recruitment age and cohort process do not reproduce the SAM state parameterization. |",
  "| F | SAM has no F state at age 0; ages 5-7 share one state and the F process is uncorrelated. The report cautions that F is unreliable for advice. | Use IID F states with shared mean levels for ages 0, 1, 2, 3-4 and 5-7; define Fbar over ages 3-5. | tinyAM estimates F at age 0 and shares means, not realized F states, across age groups. |",
  "| M and biology | M is fixed at age-specific values; the source uses annual stock, catch, landings and discard weights, fixed maturity, and zero F/M before spawning. | Retain the source M, stock weights and maturity, and set spawning time to the start of the year. | These biological inputs and SSB timing are matched. |",
  "| Outputs and advice | The report provides N-at-age, age-1 recruitment, SSB, TSB, Fbar, grouped total mortality and harvest rate; it uses harvest rate rather than F for advice. | Compare N/F at age, SSB, same-year age-0 recruitment, total abundance and fixed M. | F is retained as a fitted output but is not treated as a reliable advice metric. |",
  "",
  "The report and native SAM fit object are cached locally. This translation keeps the source observations unchanged. It does not reproduce SAM's catch observation weights, BITS Q1 correlated residuals, or full N/F process and initial-state structure."
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
    rec_settings = list(process = "iid", mu_form = ~1),
    F_settings = list(
      process = "iid",
      mu_form = ~0 + f_age_group,
      mean_ages = 3:5
    ),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_block),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block),
    silent = silent
  )
}
