## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 1993:2024
ages <- 1:4
inputs <- source$inputs
annual_rows <- function(rows, values) {
  out <- rows[match(paste(values$year, values$age), paste(rows$year, rows$age)), ]
  out$value <- values$value
  out$season <- "annual"
  out
}
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" & inputs$age %in% ages, ]
annual_catch <- stats::aggregate(value ~ year + age, catch, sum)
catch_rows <- annual_rows(catch, annual_catch)
catch_rows$transformation <- "Analysis only: sum the two source half-year age-count cells; retain zero-coded counts without equating them to total landings."
m <- inputs[inputs$type == "M" & inputs$age %in% ages, ]
m_rows <- annual_rows(m, stats::aggregate(value ~ year + age, m, sum))
m_rows$basis <- "per_year"; m_rows$unit <- "per year"
m_rows$transformation <- "Analysis only: sum the two integrated same-age half-year M hazards."
biology <- inputs[inputs$type %in% c("weight", "maturity") & inputs$season == "1" & inputs$age %in% ages, ]
biology$season <- "annual"

# With constant supplied M0 and F0=0, the juvenile index is proportional to next January N1.
m0 <- inputs[inputs$type == "M" & inputs$age == 0 & inputs$season == "2", ]
f0 <- source$outputs[source$outputs$measure == "fishing_mortality_at_age" &
  source$outputs$age == 0 & source$outputs$season == "annual", ]
stopifnot(setequal(m0$year, years), length(unique(m0$value)) == 1L, all(m0$value > 0),
  setequal(f0$year, years), all(f0$value == 0))
index <- inputs[inputs$type == "index", ]
index_mapping <- data.frame(source_survey = index$survey, source_year = index$year,
  source_age = index$age, value = index$value)
index$survey <- paste0(index$survey, " (source age ", index$age, ")")
juvenile <- index$age == 0
index$year[juvenile] <- index$year[juvenile] + 1L
index$age[juvenile] <- 1L
index$season <- "annual"
index$transformation[juvenile] <- "Analysis only: link the original age-0 index to next January age-1 abundance. Constant supplied M0 and zero age-0 F imply a constant q multiplier; observed index values are unchanged."
index_mapping$survey <- index$survey
index_mapping$model_year <- index$year
index_mapping$model_age <- index$age
index_mapping$included <- index$year %in% years & !is.na(index$value)
index <- index[index$year %in% years, ]
index$unit <- "native survey index (unit unresolved)"
sampling_times <- c("Dredge (source age 0)" = 0, "Dredge (source age 1)" = .9375,
  "Dredge2 (source age 0)" = 0, "Dredge2 (source age 1)" = .9375)
translated_inputs <- rbind(catch_rows, m_rows, biology, index)
obs <- database_to_tam_obs(assessment_id, translated_inputs, years = years, ages = ages,
  surveys = names(sampling_times), sampling_times = sampling_times, assumptions = source$assumptions)
obs$maturity$spawn_time <- 0
obs$catch$catch_sd_group <- factor(ifelse(obs$catch$age <= 2, "young", "old"))
obs$index$q_group <- obs$index$sd_group <- factor(obs$index$survey)
mapping <- index_mapping[index_mapping$included, ]
key_index <- function(z, year) paste(z$survey, year)
obs$index$source_year <- mapping$source_year[match(key_index(obs$index, obs$index$year),
  key_index(mapping, mapping$model_year))]
obs$index$source_age <- mapping$source_age[match(key_index(obs$index, obs$index$year),
  key_index(mapping, mapping$model_year))]

cw <- inputs[inputs$type == "catch_weight" & inputs$age %in% ages, ]
key <- function(x) paste(x$year, x$age, x$season)
catch$weight <- cw$value[match(key(catch), key(cw))] / 1000
catch$biomass_kg <- catch$value * 1e6 * catch$weight
mass <- stats::aggregate(biomass_kg ~ year + age, catch, sum)
mass$count <- annual_catch$value[match(paste(mass$year, mass$age), paste(annual_catch$year, annual_catch$age))] * 1e6
mean_weight <- stats::aggregate(weight ~ year + age, catch, mean)
weight <- mean_weight$weight[match(paste(mass$year, mass$age), paste(mean_weight$year, mean_weight$age))]
positive <- mass$count > 0
weight[positive] <- mass$biomass_kg[positive] / mass$count[positive]
totals <- stats::setNames(stats::aggregate(biomass_kg ~ year, mass, sum), c("year", "yield"))
landings <- inputs[inputs$type == "catch" & inputs$measure == "total_biomass", ]
# A zero-coded age array with positive landings does not establish a zero biomass observation.
totals <- totals[totals$yield > 0 | landings$value[match(totals$year, landings$year)] == 0, ]
catch_reporting <- list(weights = data.frame(year = mass$year, age = mass$age, weight = weight), totals = totals)

## Background and comparisons ----
age_plus_group <- 4
comparison_outputs <- source$outputs[!(source$outputs$measure == "fishing_mortality_at_age" & source$outputs$season != "annual"), ]
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", definition = "January N at ages 1-4+.",
    reason = "Annual dynamics replace seasonal dynamics; source M and biology are rounded."),
  F = list(status = "approximate", definition = "Annual integrated F at ages 1-4+.",
    reason = "Independent random-walk F increments replace common effort-proportional states. Exact zero F in 2009/2024 is not enforced; zero denominators have no percent difference."),
  F_bar = list(status = "approximate", definition = "January-N-weighted annual F over ages 1-2.",
    reason = "Recalculated from official annual F/N, rather than native arithmetic Fbar. Terminal source F is zero."),
  M = list(status = "matched", definition = "Sum of the supplied integrated seasonal M hazards.",
    reason = "Fixed-M input consistency only; M is not estimated."),
  M_bar = list(status = "approximate", definition = "January-N-weighted fixed annual M at ages 1-4+.",
    reason = "Any difference arises from N weights, not M estimation."),
  recruitment = list(status = "approximate", definition = "January age-1 abundance in the same calendar year.",
    reason = "Reference uses official January N1; native July age-0 recruitment remains separate. Reference years are not shifted."),
  ssb = list(status = "approximate", definition = "January SSB at ages 1-4+ under shared first-half biology.",
    reason = "Annual dynamics approximate the two-season assessment."),
  abundance = list(status = "approximate", definition = "January abundance summed over ages 1-4+.",
    reason = "July age-0 recruits are not added to January abundance."),
  biomass = list(status = "approximate", definition = "January biomass at ages 1-4+ under shared first-half weights.",
    reason = "Rounded source weights; native aggregates remain separate source context."))
background <- c(
  "### Sandeel Area 4: official 2025 SMS-effort assessment", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1993-2024, two half-years, July recruitment at age 0 and January ages 1-4+. | Calendar-year model, January ages 1-4+. | Newer 2026 category-5 advice and revised survey scaling are kept separate. |",
  "| N and recruitment | Free initial January abundance, deterministic older survival; weak July hockey-stick recruitment, fixed 44,716 t breakpoint, weight 0.05. | Free initial abundance, no older-age N process, IID January age-1 recruitment about an estimated mean. | Recruitment date, mean relationship, likelihood weight and integration differ. |",
  "| Juvenile survey | Dredge age 0 averages abundance over annual interval 0.875-1. Its next January age-1 survivors have zero age-0 F and constant supplied M. | Keep original index values and a separate q/SD; predict them from next January N1 at time zero. Original survey year/age are retained in the data. | A constant known survival/averaging factor is absorbed by the free q. This is an observation-mean re-expression, not a shift of reference recruitment. The final 2024 juvenile index is outside the historical state window and omitted. |",
  "| Older survey | Separate Dredge eras (1999-2003, 2008-2024), age 1, averaged over the same autumn interval. | Keep original age-1 values at midpoint 0.9375, with separate q/SD by era. | A midpoint approximates the integrated prediction. Dredge age 1 in 2003 remains missing. |",
  "| F and effort | Common standardized-effort driver, constant age/season coefficients, ages 3/4+ shared; exact zero F in 2009 and 2024. | Independent random-walk changes in annual log F; no effort floor or covariate. | This relaxes source state sharing/proportionality. Finite log F cannot enforce structural zeros. Source effort includes a verified second-half header correction, not a fitted-data adjustment. |",
  "| M and biology | Fixed multispecies averages by age/season, supplied seasonal weights/maturity; January SSB. | Sum the seasonal M hazards; January biology and spawning time zero. | Annual mortality is retained, while exposure within the year is simplified. Fixed M is an input check. |",
  "| Observation errors | Separate q/SD by survey era and source age; survey SDs estimated with a 0.3 lower-bound penalty. Catch log-SDs profiled by age group and season. | Four q and index-SD groups; catch SDs estimated for ages 1-2 and 3-4+. | The survey floor and profiled seasonal catch SD construction are not reproduced. Boundary estimates at 0.3 are not treated as known errors. |",
  "| Catch and yield | Seasonal positive age-count likelihood, age 0 excluded. Zero-coded counts can coexist with positive total landings. Native Yield sums predictions only at positive-observation cells. | Sum seasonal age counts and convert million to thousand fish. Catch reporting uses observed-catch-weighted seasonal weights, or duration means for zero cells. | Do not equate zero-coded age arrays to zero landings. Biomass observations are left missing in those years unless SAG also reports zero. Unrestricted tinyAM yield is not the native conditional Yield. |",
  "| Comparisons | January adults, July recruitment, annual F and native aggregates. | Same-year January N1 recruitment and January ages 1-4+ common-definition aggregates. | No cohort-year shift is applied to outputs. Derived sums have no invented uncertainty; terminal F percent differences are unavailable for official zero F. |",
  "", "This is a limited annual approximation. Native data files, exact zero-catch flags and the fitted object are missing. The source fit has 50 parameters; tinyAM does not reproduce its likelihood. Defaults are tried first, with one aligned source-start retry only after numerical failure."
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
    f <- .translation_surface(comparison_outputs[comparison_outputs$measure == "fishing_mortality_at_age", ], years, ages, "Official annual F")
    log_f <- matrix(0, length(years), length(ages), dimnames = list(year = years, age = ages))
    positive <- is.finite(f$obs) & f$obs > 0
    log_f[cbind(match(f$year[positive], years), match(f$age[positive], ages))] <- log(f$obs[positive])
    list(log_r0 = log(n[1L, 1L]), log_r = log(n[-1L, 1L]), log_n0 = log(n[1L, -1L]), log_f = log_f)
  }, attempt_file = file.path("analysis/comp_assessments/results/cache", assessment_id, "initialization_attempts.rds"))
  if (any(!is.finite(fit$pop$ssb$upr))) {
    note <- "Exploratory result: SSB has non-finite upper uncertainty limits. Numerical convergence alone does not establish a usable abundance scale."
    cli::cli_warn(note)
    background <- c(paste0("**", note, "**"), "", background)
  }
}
