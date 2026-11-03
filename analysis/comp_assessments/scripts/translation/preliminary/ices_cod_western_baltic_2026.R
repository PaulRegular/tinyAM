## Observations ----
years <- 1985:2023
ages <- 1:7
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("BITS Q1", "BITS Q4")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("BITS Q1" = 0.1653494, "BITS Q4" = 0.8477816),
  exclude_index_years = list("BITS Q1" = 1996:2000, "BITS Q4" = 1999:2000),
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0
obs$index <- obs$index[obs$index$age <= 3, , drop = FALSE]

obs$index$q_key <- interaction(
  obs$index$survey, obs$index$age, drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- interaction(
  obs$index$survey,
  ifelse(obs$index$age == 1, "age1",
    ifelse(obs$index$age == 2, "age2", "age3plus")),
  drop = TRUE, lex.order = TRUE
)
obs$catch$sd_block <- factor(ifelse(obs$catch$age == 1, "age1", "age2plus"))
obs$pop_covariates$f_age_group <- factor(
  ifelse(obs$pop_covariates$age == 1, "age1",
    ifelse(obs$pop_covariates$age == 2, "age2",
      ifelse(obs$pop_covariates$age <= 4, "age3_4", "age5plus"))),
  levels = c("age1", "age2", "age3_4", "age5plus")
)

## Comparisons ----
age_plus_group <- 7L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Sum of N at ages 1-7+ over 1985-2023, in thousand fish.",
    reason = "Both series use the same reported ages and years; uncertainty for the derived total is unavailable."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 recruitment in the same calendar year.",
    reason = "Age 1 is the youngest modeled age and the source report's recruitment age."
  ),
  ssb = list(
    status = "approximate",
    definition = "SSB over ages 1-7+ using source N-at-age, stock weights and maturity on the translated grid.",
    reason = "Age-0 maturity is not reported, so this common-age calculation is not assumed to equal native SAM SSB."
  ),
  M = list(
    status = "matched",
    definition = "The fixed source M-at-age vector at ages 1-7.",
    reason = "M is supplied from the 2026 report and is not estimated by tinyAM."
  )
)

background <- c(
  "### Western Baltic cod: 2026 WGBFAS assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years | SAM report covers catch through 2025 and estimates through 2026. | Fit 1985-2023. | The report does not provide age-7+ stock weights for 2024-2025; the restricted period retains a complete biological grid. |",
  "| Ages | Reported population is ages 0-7+, with age-1 recruitment. | Fit ages 1-7+, using age 1 as recruitment. | The report does not provide maturity at age 0; age-0 observations and the FEJUCS index are omitted. |",
  "| N | SAM process and initial-state settings were not recovered. | IID N process, exponential initialization and IID recruitment. | Simple approximation; exact 2026 process settings are unknown. |",
  "| F and M | F is considered unreliable for advice; the report supplies an age-specific M vector. | IID F process; fixed age-specific source M. | F-state settings cannot be verified; no F-at-age output is published. |",
  "| Catch | Total catch-at-age, ages 1-7+, in thousand fish. | Use the reported catch series through 2023. | Native observation weights and likelihood encoding are unavailable. |",
  "| Surveys | BITS Q1 ages 1-4+ and Q4 ages 0-4+; the assessment uses years 2001 onward. | Retain Q1 and Q4 at individual ages 1-3 from 2001 onward. | The pooled 4+ index group cannot be distributed across ages 4-7 by tinyAM, and age-0 FEJUCS is outside the translated grid. Sampling times carry forward from the 2025 native SAM run because 2026 timing was not recovered. |",
  "| Biology | Fixed M and maturity; SSB at the start of year with zero pre-spawning mortality. | Use source M, maturity and stock weights over the fit grid; set spawning time to zero. | The restricted SSB comparison excludes age 0 and is a common-definition approximation. |",
  "| Outputs | Reported N-at-age, recruitment, SSB, total biomass and harvest rate; no F-at-age. | Compare common-age N, recruitment and SSB through 2023. | Native 2026 likelihood settings, fitted F, and terminal weight cells are unavailable. |",
  "",
  "This is a restricted approximation to the 2026 report, not a reconstruction of its SAM fit. The latest assessment's report tables remain the source of the observations; no source observations were altered."
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