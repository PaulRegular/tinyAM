## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 1987:2025
ages <- 1:4
inputs <- source$inputs

catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" &
  inputs$age %in% ages, ]
annual_catch <- stats::aggregate(value ~ year + age, catch, sum)
annual_rows <- function(rows, values) {
  out <- rows[match(paste(values$year, values$age), paste(rows$year, rows$age)), ]
  out$value <- values$value
  out$season <- "annual"
  out
}
catch_rows <- annual_rows(catch, annual_catch)
catch_rows$transformation <- "Analysis only: sum first- and second-half catch numbers at the same year and age; retain 2025 zeros."

m <- inputs[inputs$type == "M" & inputs$age %in% ages, ]
m_rows <- annual_rows(m, stats::aggregate(value ~ year + age, m, sum))
m_rows$basis <- "per_year"
m_rows$unit <- "per year"
m_rows$transformation <- "Analysis only: add the two integrated half-year M hazards at the same year and age."
biology <- inputs[inputs$type %in% c("weight", "maturity") &
  inputs$season == "1" & inputs$age %in% ages, ]
biology$season <- "annual"
index <- inputs[inputs$type == "index" & inputs$age %in% ages, ]
# Retain the numerical survey scale; these are not absolute fish counts.
index$unit <- "native survey index (unit unresolved)"
translated_inputs <- rbind(catch_rows, m_rows, biology, index)
obs <- database_to_tam_obs(assessment_id, translated_inputs,
  years = years, ages = ages, sampling_times = c(Dredge = .9375, Acoustic = .3),
  surveys = c("Dredge", "Acoustic"), assumptions = source$assumptions)
obs$maturity$spawn_time <- 0
obs$catch$catch_sd_group <- factor(ifelse(obs$catch$age <= 2, "young", "old"))
obs$index$q_group <- factor(ifelse(obs$index$survey == "Dredge", "Dredge1",
  paste0("Acoustic", pmin(obs$index$age, 3))))
obs$index$sd_group <- factor(obs$index$survey)

cw <- inputs[inputs$type == "catch_weight" & inputs$age %in% ages, ]
key <- function(x) paste(x$year, x$age, x$season)
catch$weight <- cw$value[match(key(catch), key(cw))] / 1000
catch$biomass_kg <- catch$value * 1e6 * catch$weight
annual_mass <- stats::aggregate(biomass_kg ~ year + age, catch, sum)
annual_mass$count <- annual_catch$value[match(paste(annual_mass$year, annual_mass$age),
  paste(annual_catch$year, annual_catch$age))] * 1e6
mean_weight <- stats::aggregate(weight ~ year + age, catch, mean)
weight <- mean_weight$weight[match(paste(annual_mass$year, annual_mass$age),
  paste(mean_weight$year, mean_weight$age))]
positive <- annual_mass$count > 0
weight[positive] <- annual_mass$biomass_kg[positive] / annual_mass$count[positive]
stopifnot(all(is.finite(weight)), all(weight > 0),
  all(annual_mass$biomass_kg[!positive] == 0))
# With zero observed catch, use the equal-duration mean only for yield reporting.
catch_reporting <- list(
  weights = data.frame(year = annual_mass$year, age = annual_mass$age, weight = weight),
  totals = stats::setNames(stats::aggregate(biomass_kg ~ year, annual_mass, sum), c("year", "yield")))

## Background and comparisons ----
age_plus_group <- 4
comparison_outputs <- source$outputs[!(source$outputs$measure == "fishing_mortality_at_age" &
  source$outputs$season != "annual"), ]
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", definition = "January abundance at ages 1-4+.",
    reason = "Annual dynamics approximate the seasonal model; July age-0 states are omitted."),
  F = list(status = "approximate", definition = "Annual integrated F at ages 1-4+.",
    reason = "Independent annual log-F changes replace effort-proportional F. Exact zero F in 2025 is not enforced. Zero source F has no percent difference."),
  F_bar = list(status = "approximate", definition = "January-N-weighted annual F over ages 1-2.",
    reason = "Recalculated from annual F and January N, rather than native arithmetic Fbar. Official F in 2025 is zero; terminal percent difference is unavailable."),
  M = list(status = "matched", definition = "Sum of supplied integrated half-year M hazards at each year and age.",
    reason = "Fixed-M input consistency only; time-varying M is supplied, not estimated."),
  M_bar = list(status = "approximate", definition = "January-N-weighted supplied annual M over ages 1-4+.",
    reason = "Differences arise from N weights, not M estimation."),
  recruitment = list(status = "approximate", definition = "January age-1 abundance in the same calendar year.",
    reason = "Official January N1 replaces native July age-0 recruitment for this comparison; no cohort-year shift."),
  ssb = list(status = "approximate", definition = "January SSB over ages 1-4+ under shared first-half biology.",
    reason = "Annual dynamics replace seasonal dynamics, with rounded source inputs."),
  abundance = list(status = "approximate", definition = "January abundance summed over ages 1-4+.",
    reason = "Later July age-0 recruits are excluded from this simultaneous population."),
  biomass = list(status = "approximate", definition = "January biomass over ages 1-4+ under shared first-half weights.",
    reason = "Reported weights are rounded; native aggregates remain separate source context."))
background <- c(
  "### Sandeel Area 3r: official 2026 SMS-effort assessment", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Active window 1987-2025, two half-years; July age 0 and January ages 1-4+. | Calendar years 1987-2025, ages 1-4+. | The loader also requests 1986, but the active start is 1987. July age-0 observations are outside this annual model. |",
  "| N and recruitment | Free initial January abundance; deterministic older-cohort survival; July hockey-stick recruitment, fixed SSB breakpoint 72,713 t, weight 0.05. | Free initial abundance, no older-age N process, IID age-1 recruitment about an estimated mean. | The different entry date and downweighted hockey-stick construction are not reproduced. |",
  "| F and closure | F proportional to seasonal standardized effort; constant age/season coefficients, ages 3-4+ shared. Effort, catch and F are exactly zero in 2025. | Independent random-walk changes in annual log F, without an effort covariate. | Finite log-F states cannot represent structural zero F. No pseudo-effort is added and 2025 is retained. This relaxes common effort variation and shared age-3/4+ states. An earlier correlated-increment approximation failed numerical checks near its correlation boundary. |",
  "| M | Fixed supplied, time-varying multispecies age/season hazards. | Sum both seasonal hazards at each year and age; no M process. | Annual mortality is preserved; within-year exposure is simplified. Matching M is an input check. |",
  "| Catch | Positive seasonal numbers at ages 1-4+, with season/age SD groups; zero-catch flags not recovered. | Sum seasonal counts, convert million to thousand fish; SD groups ages 1-2 and 3-4+. All 2025 zeros are retained. | The lognormal likelihood excludes zero catch, so the closure is not enforced by these observations. Published positive eligible cells number 240 versus 244 in the native diagnostic. |",
  "| Surveys | Dredge ages 0-1 over annual interval 0.875-1; Acoustic ages 1-4+ over 0.25-0.35. | Dredge age 1 at 0.9375 and Acoustic ages 1-4+ at 0.3, original scaled values. Acoustic q shared at ages 3/4+; one Acoustic SD. | Midpoints approximate the integrated predictions. Dredge age 0 is omitted; shared q does not pool acoustic observations. |",
  "| Biology and yield | Seasonal weights/maturity; January SSB. | First-half biology and spawning time zero; seasonal observed-catch-weighted catch weights. | For zero-catch cells, equal-half-year mean weights are used only for predicted-yield reporting. No fitted state modifies source observations. All-age SAG catch remains separate. |",
  "| Comparisons | January older N and July age-0 recruitment; annual F sums seasonal F. | January ages 1-4+, same-year age-1 recruitment, common-definition aggregates. | Native July recruitment/aggregates remain source context. Derived aggregates have no invented uncertainty. Terminal F percent difference is unavailable because official F is zero. |",
  "", "This is a limited annual approximation. The 2025 closure cannot be represented exactly and must be considered when interpreting terminal F and abundance. Exact catch flags, unrounded arrays and the actual production builder revision are unresolved. Defaults are tried first; one aligned source-start retry is permitted only after numerical failure."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "rw", mu_form = NULL, mean_ages = 1:2),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + catch_sd_group),
  index_settings = list(q_form = ~0 + q_group, sd_form = ~0 + sd_group),
  silent = silent, ...)
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages)
    f_rows <- comparison_outputs[comparison_outputs$measure == "fishing_mortality_at_age", ]
    f <- .translation_surface(f_rows, years, ages, "Official annual F")
    # Keep package log-F starts of zero (F = 1) where source F cannot be logged.
    log_f <- matrix(0, length(years), length(ages), dimnames = list(year = years, age = ages))
    positive <- is.finite(f$obs) & f$obs > 0
    log_f[cbind(match(f$year[positive], years), match(f$age[positive], ages))] <- log(f$obs[positive])
    list(log_r0 = log(n[1L, 1L]), log_r = log(n[-1L, 1L]),
      log_n0 = log(n[1L, -1L]), log_f = log_f)
  }, attempt_file = file.path("analysis/comp_assessments/results/cache", assessment_id,
    "initialization_attempts.rds"))
  if (any(!is.finite(fit$pop$ssb$upr))) {
    review_note <- "Exploratory result: SSB has non-finite upper uncertainty limits. The abundance scale is poorly resolved; numerical convergence alone does not make this a usable assessment comparison."
    cli::cli_warn(review_note)
    background <- c(paste0("**", review_note, "**"), "", background)
  }
}
