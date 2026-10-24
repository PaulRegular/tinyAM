## Observations ----
assessment_id <- source$assessment$assessment_id[[1L]]
years <- 1963:2024
ages <- 1:9
surveys <- c(
  "NEFSC spring bottom trawl survey",
  "NEFSC fall bottom trawl survey",
  "ASMFC shrimp trawl survey",
  "NEFSC spring bottom longline survey",
  "NEFSC fall bottom longline survey"
)

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 0.25
obs$index$q_key <- interaction(obs$index$survey, obs$index$age, drop = TRUE)
q_has_positive_obs <- ave(
  is.finite(obs$index$obs) & obs$index$obs > 0,
  obs$index$q_key,
  FUN = any
)
uninformative_index_rows <- sum(!q_has_positive_obs)
obs$index <- droplevels(obs$index[q_has_positive_obs, , drop = FALSE])

## Comparison and background ----
age_plus_group <- 9
comparison_outputs <- source$outputs
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "matched",
    definition = "Annual SSB at 0.25 of the year.",
    reason = "The source weights, maturity, fixed M, and spawning time are used; population estimates and likelihoods differ."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 abundance in the same calendar year.",
    reason = "Both series report recruitment at age 1; tinyAM uses a different stock-recruit and process formulation."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Source arithmetic Fbar over ages 6-9+; tinyAM reports an abundance-weighted mean.",
    reason = "Official F-at-age and N-at-age were not recovered, so the two averaging definitions cannot be reconciled."
  ),
  M = list(
    status = "matched",
    definition = "Fixed source M input at ages 1-9+.",
    reason = "This checks the fixed input mapping and is not an estimated-M comparison."
  )
)

background <- c(
  "### White Hake: official 2025 NEFSC Management Track Run_005",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | ASAP 3.0, 1963-2024, ages 1-9+; one combined-sex stock. | Fit 1963-2024, ages 1-9, with age 9 as the plus group. | The fitted historical grid is retained. |",
  "| Catch | Total catch biomass with CV is available in every year; age data with ESS 58 are available for 1989-2024. | Use the source catch-at-age values from 1989 onward in tinyAM's lognormal catch likelihood. | The aggregate catch likelihood and age-composition likelihood are combined into one direct age-specific likelihood; total catch biomass before 1989 is not fitted. |",
  "| Surveys and q | Five aggregate indices and age compositions use separate CV and effective-sample-size information; q is survey-specific with age selectivity. | Fit age-specific number indices with q by survey and age and an estimated observation SD shared within survey. | tinyAM does not reproduce the separate aggregate-index and composition likelihoods or use composition ESS; source aggregate CVs are not assigned to age-specific observations. |",
  "| N and recruitment | ASAP fixes steepness at 1, giving constant expected recruitment; annual recruitment deviations have zero penalty weight. Older cohorts survive deterministically. Initial older-age abundance is estimated with a CV 0.2 penalty around equilibrium survival. | Use exponential initial abundance, IID N-process deviations, and IID recruitment around a constant log mean. | tinyAM regularizes recruitment, adds older-cohort process variation, and does not reproduce the source initial-age penalty or its sum-to-zero fixed-effect recruitment vector. |",
  "| F and M | Two age-selectivity blocks; fixed M = 0.20 per year at ages 1-9+. | Use an AR1 F process with age-specific means and the source fixed M surface. | F process and selectivity are simplified; M is an input consistency check. |",
  "| Biology and SSB | Annual stock weights and maturity; SSB is calculated at 0.25 of the year. | Use the source annual weights, maturity, fixed M, and spawning time. | Population states and observation likelihoods still differ. |",
  "| Reported outputs | Annual SSB, age-1 recruitment, and arithmetic Fbar are tabulated; fitted N-at-age and F-at-age are unavailable. | Compare SSB and recruitment; retain fixed M for input checking. | Fbar is not compared as equivalent because its averaging definition cannot be recovered from age-specific source states. No total abundance is derived without N-at-age. |",
  "",
  paste0("The canonical database preserves all five surveys, total catch biomass, CVs, age-composition sample sizes, and source zeros. tinyAM treats zero lognormal observations as missing; ", uninformative_index_rows, " all-zero survey-age rows are omitted from the translated fit because their survey-age q has no positive observation and cannot be estimated. The source report's catch-at-age table labels and the native input scale are inconsistent; the database records the native values as thousands of fish based on the weight cross-check and keeps that caveat visible in the source review."),
  "The package initialization is tried first. If numerical checks fail, one retry uses source recruitment and annual Fbar only as starting values; missing N-at-age and F-at-age starts retain package defaults, and the observations and model settings remain unchanged."
)

## Model ----
source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
model_settings <- list(
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~ 1),
  F_settings = list(
    process = "ar1", mu_form = ~ 0 + factor(age), mean_ages = 6:9
  ),
  M_settings = list(
    process = "off", mu_form = NULL, mu_supplied = ~ M_assumption
  ),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ survey)
)
fit_model <- function(...) do.call(tinyAM::fit_tam, c(
  list(data = obs, years = years, ages = ages, silent = silent),
  model_settings,
  list(...)
))
source_starts <- function() {
  dat <- do.call(tinyAM::prepare_tam, c(
    list(data = obs, years = years, ages = ages), model_settings
  ))
  starts <- tinyAM::make_par(dat)
  recruitment <- source$outputs[
    source$outputs$measure == "recruitment", c("year", "value", "unit"), drop = FALSE
  ]
  fbar <- source$outputs[
    source$outputs$measure == "Fbar", c("year", "value"), drop = FALSE
  ]
  if (any(recruitment$unit != "thousand fish")) {
    cli::cli_abort("White Hake recruitment starts require source values in thousand fish.")
  }
  r <- as.numeric(recruitment$value[match(years, as.integer(recruitment$year))])
  f <- as.numeric(fbar$value[match(years, as.integer(fbar$year))])
  if (any(!is.finite(r) | r <= 0) || any(!is.finite(f) | f <= 0)) {
    cli::cli_abort("Source recruitment or Fbar is incomplete for starting values.")
  }
  starts$log_r0 <- log(r[[1L]])
  starts$log_r[] <- log(r[-1L])
  starts$log_mu_f[] <- log(mean(f))
  starts$log_f[,] <- log(f)
  starts[c("log_r0", "log_r", "log_mu_f", "log_f")]
}
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = source_starts,
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache", assessment_id,
      "initialization_attempts.rds"
    )
  )
}
