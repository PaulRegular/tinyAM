## Observations ----
years <- 1981:2025
ages <- 0:7
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("WCIBTS.Q1", "SCO.Q1", "SWC.Q4")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c(
    "WCIBTS.Q1" = 0.125,
    "SCO.Q1" = 0.125,
    "SWC.Q4" = 0.8530996
  ),
  assumptions = source$assumptions
)

obs$catch$sd_block <- factor(
  ifelse(obs$catch$age == 0, "age0",
    ifelse(obs$catch$age == 7, "age7_plus", "ages1_6")
  ),
  levels = c("age0", "ages1_6", "age7_plus")
)

obs$index$q_age <- ifelse(
  obs$index$survey == "SWC.Q4",
  ifelse(obs$index$age <= 5, paste0("age", obs$index$age), "age6_7"),
  ifelse(obs$index$age <= 4, paste0("age", obs$index$age), "age5_6")
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$relative_sd <- 1
has_weight <- is.finite(obs$index$relative_precision_weight)
obs$index$relative_sd[has_weight] <-
  1 / sqrt(obs$index$relative_precision_weight[has_weight])

## Comparison ----
age_plus_group <- 7L
comparison_outputs <- source$outputs

background <- c(
  "### West of Scotland whiting: 2025 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM models 1981–2025 at ages 0–7+, with catch through 2024. | Fit 1981–2025 at ages 0–7+, preserving the age-7 plus group. | The model grid is retained; 2025 has no catch observations. |",
  "| Recruitment and N | Recruitment follows the SAM random-walk specification; older N process variances are shared across ages 1–7. | Use a random-walk recruitment process and an IID N process with exponential initialization. | The N transition covariance and initial-state treatment are simplified. |",
  "| F | SAM uses age-specific F states, random-walk time evolution, and AR(1) correlation across age increments. Fbar is ages 1–3. | Use correlated random walks across ages and calculate Fbar over ages 1–3. | tinyAM uses one shared process SD and does not reproduce all SAM variance settings. |",
  "| M and biology | M is fixed at the supplied age-year values; stock and catch-related weights, maturity, and spawning-time fractions are supplied. | Retain the same M and biological input surfaces. | No M process is estimated in either representation. |",
  "| Catch | Total catch-at-age in thousands of fish, with lognormal errors grouped for age 0, ages 1–6, and age 7+. | Retain the catch series and use the same three age-based SD groups. | The source likelihood's remaining observation structure is simplified. |",
  "| Surveys | Two Q1 fleets and one Q4 fleet remain separate. Timing is 0.125, 0.125, and 0.8530996. q is shared over the source age groups; Q4 residuals have AR(1) age correlation. | Retain the three fleets, timing, q groups, and available relative precision weights. | tinyAM uses independent survey residuals; the Q4 age correlation is omitted. |",
  "",
  "The translated model uses source N, F, and q outputs only for comparison. No source starting values are supplied; tinyAM uses its normal initialization."
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
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mean_ages = 1:3),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~0 + sd_block),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + survey,
      sd_supplied = ~relative_sd
    ),
    silent = silent
  )
}
