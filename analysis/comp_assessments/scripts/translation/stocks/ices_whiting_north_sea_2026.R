## Observations ----
years <- 1978:2026
ages <- 0:8
assessment_id <- source$assessment$assessment_id

inputs <- source$inputs
m_rows <- inputs[inputs$assessment_id == assessment_id &
  inputs$type == "M" &
  inputs$measure == "natural_mortality_at_age", , drop = FALSE]
last_m_year <- max(as.integer(m_rows$year), na.rm = TRUE)
extra_years <- setdiff(years, unique(as.integer(m_rows$year)))
if (length(extra_years)) {
  m_last <- m_rows[as.integer(m_rows$year) == last_m_year, , drop = FALSE]
  extra_m <- m_last[rep(seq_len(nrow(m_last)), each = length(extra_years)), , drop = FALSE]
  extra_m$year <- rep(extra_years, times = nrow(m_last))
  extra_m$source_type <- "translation_assumption"
  extra_m$transformation <- "Held the final numerical M input through the remaining model years as a fixed tinyAM assumption."
  extra_m$notes <- paste(
    "Translation-only baseline extension. The accepted SAM run estimates",
    "later M values with its mortality GMRF."
  )
  inputs <- rbind(inputs, extra_m)
}

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)

obs$index$q_key <- interaction(obs$index$survey, obs$index$q_block,
  drop = TRUE, lex.order = TRUE
)


source_N <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages, 1000)
source_F <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age", years, ages)
source_M <- .translation_start_surface(source$outputs, "mortality", "natural_mortality_at_age", years, ages)
names(dimnames(source_N)) <- names(dimnames(source_F)) <- names(dimnames(source_M)) <- c("year", "age")
q_design <- stats::model.matrix(~ 0 + q_key, data = obs$index)
survey_year <- match(as.character(obs$index$year), rownames(source_N))
survey_age <- match(as.character(obs$index$age), colnames(source_N))
survey_rows <- cbind(survey_year, survey_age)
n_at_survey <- source_N[survey_rows] * exp(
  -(source_F[survey_rows] + source_M[survey_rows]) * obs$index$samp_time
)
q_start <- vapply(seq_len(ncol(q_design)), function(j) {
  selected <- q_design[, j] == 1
  stats::median(log(obs$index$obs[selected] / n_at_survey[selected]))
}, numeric(1))
names(q_start) <- colnames(q_design)

start_par <- list(
  log_r0 = log(source_N[1, "0"]),
  log_r = setNames(log(source_N[-1, "0"]), as.character(years[-1])),
  log_n = log(source_N[-1, as.character(ages[-1]), drop = FALSE]),
  log_f = log(source_F),
  log_q = q_start
)


## Background and comparisons ----

age_plus_group <- 8

comparison_scales <- c(
  N = 0.001, recruitment = 0.001, ssb = 0.001, biomass = 0.001, abundance = 0.001,
  biomass_at_age = 0.001, F_bar = 1
)


background <- c(
  "### North Sea whiting: accepted 2026 SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | Catch covers 1978–2025; IBTS Q1 covers 1983–2026 and Q3 covers 1991–2025. | Fit 1978–2026. | The accepted 2026 survey observation and native terminal states are retained. |",
  "| Ages | Ages 0–8+, with recruitment at age 0. | Ages 0–8+, with age 8 as the population plus group. | The source population age range is retained. Survey plus groups begin at age 6 and are not yet represented exactly in the observation model. |",
  "| N | Random-walk recruitment and cohort-survival processes, with age-specific sharing of survival-process variances. | Random-walk recruitment, IID cohort-survival deviations with one shared SD, and exponential initial abundance. | The recruitment process is retained in simplified form, but older-cohort deviations are independent rather than temporally correlated; variance sharing and initial-state treatment also differ. |",
  "| F | Random-walk F with AR(1) correlation across ages in annual F increments; Fbar covers ages 2–5. | Correlated F random walks (cor_rw), with AR(1) correlation across ages in annual log-F increments and one shared process SD; Fbar covers ages 2–5. | Cross-age correlation is now represented, but SAM's precise F variance grouping and other process details are not necessarily reproduced. |",
  "| M | WGSAM numerical M observations are available through 2022; SAM uses a mortality GMRF to estimate M, including later years. | Use the numerical M surface through 2022 and hold its final-year values fixed through 2026. No M process is estimated. | The fixed extension replaces SAM's estimated later-year M variation; an estimated tinyAM M process did not yield a usable fit. |",
  "| Catch | Total catch-at-age in thousands of fish, with independent lognormal observation errors and source-specific age variance groups. | Convert catch numbers to individual fish and fit independent lognormal errors with one shared SD. | SAM's catch observation-SD age groups are simplified to a single parameter. |",
  "| Indices | IBTS Q1 ages 1–6+ and Q3 ages 0–6+, sampled at fractions 0.125 and 0.625 of the year. Catchability varies by survey and age; observation residuals have AR(1) age correlation. | Retain source observations, survey timing and survey-by-age q groups. Estimate separate baseline lognormal SDs by survey and apply supplied relative SD factors, with independent age residuals. | Cross-age observation correlation is omitted. The native survey 6+ observations are currently represented at age 6 rather than as predictions summed over ages 6–8+. |",
  "| Weights and maturity | Annual stock and catch weights, smoothed maturity, and zero pre-spawning mortality fractions. Biological age-6+ inputs are expanded across older model ages. | Retain supplied biological values and start-year spawning timing. | The published biological input surfaces are retained, including their source age-group conventions. |",
  "",
  "The revised translation uses IID cohort-survival deviations and correlated F random walks, replacing the earlier random-walk N process and independent F random walks. The correlated F process more closely matches SAM's age correlation of annual fishing-mortality increments.",
  "",
  "Accepted SAM N and F surfaces provide starting values only. Initial q values are derived from accepted N, F and M estimates together with observed survey indices; all states and catchability parameters are re-estimated by tinyAM.",
  "",
  "The native number scale is thousands of fish, so catch and population counts are translated to individual fish. Biomass comparisons are scaled back to tonnes.",
  "",
  "The reported 2026 assessment table contains SSB but not the other population summaries for 2026; native fitted state surfaces are nevertheless available through that year.",
  "",
  "The 2022 numerical M surface is extended through 2026 solely as a translation assumption. The canonical accepted-model inputs and outputs remain unchanged."
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
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + survey, sd_supplied = ~relative_sd),
    silent = silent,
    start_par = start_par
  )
}
