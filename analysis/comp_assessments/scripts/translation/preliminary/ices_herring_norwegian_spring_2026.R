## Observations ----
years <- 1988:2025
ages <- 2:12
inputs <- source$inputs
maturity <- inputs$type == "maturity" & inputs$measure == "maturity_at_age"
inputs$year[maturity] <- as.character(
  as.integer(inputs$year[maturity]) + as.integer(inputs$age[maturity])
)
inputs$year_basis[maturity] <- "calendar_year"
natural_mortality <- inputs$type == "M" &
  inputs$measure == "natural_mortality_at_age"
inputs <- inputs[!natural_mortality | as.integer(inputs$age) %in% ages, ]

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = c(
    "NASF_1988_2008", "NASF_2015_onward", "IESNS_Barents",
    "IESNS_Norwegian_Sea", "BESS"
  ),
  assumptions = source$assumptions
)

model_settings <- list(
  N_settings = list(process = "iid", init = "free"),
  F_settings = list(process = "rw", mu_form = NULL, mean_ages = 5:12),
  M_settings = list(
    process = "off", mu_form = NULL, mu_supplied = ~M_assumption
  ),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey)
)

dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = years, ages = ages), model_settings
))

## Background and comparisons ----
age_plus_group <- 12

background <- c(
  "### Norwegian spring-spawning herring: official 2026 SAM assessment",
  "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | Catch-at-age is fitted through 2025; reported N and F estimates extend through 2026. | Fit 1988–2025. | The fit uses the available catch period; comparisons use the common fitted years. |",
  "| Ages | Ages 2–12+, with recruitment at age 2. Source catch and biological tables retain detail through 15+. | Ages 2–12, with age 12 as the plus group. | Catch and index numbers older than 12 are summed; age-12 weight and maturity represent the plus group. |",
  "| N | The official assessment uses SAM's age-structured population process; the exact configuration was not recovered. | IID abundance process with freely estimated initial abundance. | tinyAM does not reproduce the unavailable SAM transition and initial-state settings. |",
  "| F | SAM estimates age-specific fishing mortality; the detailed process configuration was not recovered. | Random-walk F with age-specific states and Fbar ages 5–12. | The F process and age sharing approximate the source model. |",
  "| M | Standard M is 0.9 per year at age 2 and 0.15 per year at ages 3+; the report notes annual stock-annex deviations. | Supply the standard age pattern as fixed M. | Annual stock-annex deviations were not recovered. |",
  "| Catch | One aggregate catch-at-age series in thousand fish, with age-specific sampling errors. | Fit the source catch-at-age in tinyAM's lognormal likelihood with one shared observation SD. | tinyAM estimates its own error rather than applying SAM's external relative errors. |",
  "| Index | Numerical NASF, IESNS Barents, IESNS Norwegian Sea, and BESS abundance-at-age series are available. IESNS Barents is fitted at age 2 only; RFID values are shown only in a figure. | Retain the four numerical series with q by survey and age and observation SD by survey. | RFID is omitted because numeric values could not be recovered; SAM's q sharing and external error weights are not reproduced. |",
  "| Weights and maturity | Annual stock and catch weights; maturity varies by birth cohort. | Use annual stock weights and map maturity cohort to calendar year using cohort = year - age. | Age-12 weight and maturity represent the plus group because older single-age abundance is unavailable. |",
  "| SSB | The report supplies SSB and Fbar summaries; Fbar covers ages 5–12+. | Calculate SSB using tinyAM's fixed M surface and default spawning time. | The source's spawning-time and annual M adjustments are not fully available. |",
  "",
  "Relative standard errors are retained in the database but are not applied in this fit. The fit uses package-default starting values first; one source-based retry is attempted only if numerical checks fail. Source estimates initialize parameters only and do not constrain tinyAM estimates."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_model <- function(...) do.call(tinyAM::fit_tam, c(
    list(data = obs, years = years, ages = ages, silent = silent),
    model_settings,
    list(...)
  ))
  source_starts <- function() {
    starts <- tinyAM::make_par(dat)
    source_N <- .translation_start_surface(
      source$outputs, "population", "numbers_at_age", years, ages, 1000
    )
    source_F <- .translation_start_surface(
      source$outputs, "mortality", "fishing_mortality_at_age", years, ages
    )
    starts$log_r0 <- log(source_N[1L, 1L])
    starts$log_n0 <- log(source_N[1L, -1L])
    starts$log_r <- log(source_N[-1L, 1L])
    starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])
    starts$log_f <- log(source_F)
    index <- obs$index
    n_at_age <- source_N[cbind(match(index$year, years), match(index$age, ages))]
    q_start <- tapply(index$obs / n_at_age, as.character(index$q_key),
      stats::median, na.rm = TRUE
    )
    starts$log_q[paste0("q_key", names(q_start))] <- log(q_start)
    starts
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = source_starts,
    attempt_file = file.path(
      .assessment_cache_root(translation_tier),
      source$assessment$assessment_id[[1L]], "initialization_attempts.rds"
    )
  )
}
