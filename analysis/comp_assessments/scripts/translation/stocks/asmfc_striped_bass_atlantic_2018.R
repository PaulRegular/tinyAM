## Observations and derived inputs ----

assessment_id <- source$assessment$assessment_id[[1L]]
years <- 1982:2017
ages <- 1:15
inputs <- source$inputs
inputs$year <- suppressWarnings(as.numeric(inputs$year))
inputs$age <- suppressWarnings(as.numeric(inputs$age))
inputs$value <- suppressWarnings(as.numeric(inputs$value))

catch_totals <- inputs[inputs$type == "catch" & inputs$measure == "total_numbers", ]
catch_props <- inputs[inputs$type == "catch" & inputs$measure == "proportion_at_age", ]
catch_key <- paste(catch_totals$year, catch_totals$fleet)
catch_props$total <- catch_totals$value[match(paste(catch_props$year, catch_props$fleet), catch_key)]
if (anyNA(catch_props$total)) stop("Catch proportions do not all have a matching fleet total.")
catch_props$value <- catch_props$value * catch_props$total
catch_age <- aggregate(value ~ year + age, catch_props, sum)
catch_rows <- inputs[rep(1L, nrow(catch_age)), , drop = FALSE]
catch_rows$type <- "catch"
catch_rows$measure <- "numbers_at_age"
catch_rows$basis <- "numbers"
catch_rows$fleet <- catch_rows$survey <- catch_rows$sex <- catch_rows$region <- catch_rows$season <- ""
catch_rows$year <- catch_age$year
catch_rows$year_basis <- "calendar_year"
catch_rows$age <- catch_age$age
catch_rows$value <- catch_age$value
catch_rows$unit <- "fish"
catch_rows$source_type <- "analysis_derived_input"
catch_rows$source_reference <- paste(unique(c(catch_totals$source_reference, catch_props$source_reference)), collapse = "; ")
catch_rows$transformation <- "Multiply each fleet's reported total removals by its reported age proportions, then sum the two fleets by year and age."
catch_rows$notes <- "Analysis-only pooled catch-at-age. The official model estimates fleet-specific F and selectivity; tinyAM fits one pooled catch series."

index_totals <- inputs[inputs$type == "index" & inputs$measure == "total_numbers", ]
index_props <- inputs[inputs$type == "index" & inputs$measure == "proportion_at_age", ]
index_map <- c(
  "MRIP" = "MRIP CPUE", "CT Trawl" = "CT LISTS", "NYOHS" = "NY OHS",
  "NJ Trawl" = "NJ OT", "DE 30 Trawl" = "DE 30'", "MD SSN" = "MD SSN",
  "DE SSN" = "DE SSN", "ChesMMAP" = "ChesMMAP"
)
index_props$source_survey <- index_props$survey
index_props$survey <- unname(index_map[index_props$source_survey])
if (anyNA(index_props$survey)) stop("An age-composition index could not be matched to its aggregate survey.")
index_key <- paste(index_totals$year, index_totals$survey)
index_props$total <- index_totals$value[match(paste(index_props$year, index_props$survey), index_key)]
index_props <- index_props[is.finite(index_props$total) & index_props$year %in% years, ]
index_props$composition_sum <- ave(index_props$value, index_props$year, index_props$survey, FUN = sum)
index_props <- index_props[is.finite(index_props$composition_sum) & index_props$composition_sum > 0, ]
index_props$age_index <- index_props$total * index_props$value / index_props$composition_sum
age_index_rows <- index_props[is.finite(index_props$age_index) & index_props$age_index > 0, ]

recruitment_index <- index_totals[
  !is.na(index_totals$age) & index_totals$year %in% years &
    index_totals$survey %in% c("NY YOY", "NY Age-1", "NJ YOY", "MD YOY", "MD Age-1", "MDVA YOY"),
]
recruitment_index$age_index <- recruitment_index$value
index_rows <- rbind(
  age_index_rows[c("year", "age", "survey", "age_index", "sampling_time", "source_reference")],
  recruitment_index[c("year", "age", "survey", "age_index", "sampling_time", "source_reference")]
)
index_rows <- index_rows[is.finite(index_rows$age_index) & index_rows$age_index > 0, ]
index_data <- inputs[rep(1L, nrow(index_rows)), , drop = FALSE]
index_data$type <- "index"
index_data$measure <- "numbers_at_age"
index_data$basis <- "index_scale"
index_data$year <- index_rows$year
index_data$year_basis <- "calendar_year"
index_data$age <- index_rows$age
index_data$value <- index_rows$age_index
index_data$fleet <- index_data$sex <- index_data$region <- index_data$season <- ""
index_data$survey <- index_rows$survey
index_data$unit <- "source index scale per effort unit"
index_data$sampling_time <- index_rows$sampling_time
index_data$source_type <- "analysis_derived_input"
index_data$source_reference <- index_rows$source_reference
index_data$transformation <- "Multiply the aggregate relative index by the normalized published age composition; recruitment indices are linked to their Table B10 age."
index_data$notes <- "Analysis-only age-specific index. Source zero proportions are omitted because tinyAM uses a lognormal likelihood; source composition rows are normalized over available ages."

maturity <- inputs[inputs$type == "maturity" & inputs$measure == "maturity_at_age", ]
female <- inputs[inputs$type == "biology" & inputs$measure == "proportion_female_at_age", ]
maturity$value <- maturity$value * female$value[match(maturity$age, female$age)]
maturity$source_type <- "analysis_derived_input"
maturity$source_reference <- paste(unique(c(maturity$source_reference, female$source_reference)), collapse = "; ")
maturity$transformation <- "Female maturity multiplied by the fixed proportion-female-at-age input."
maturity$notes <- "Effective mature-female fraction for tinyAM's single maturity surface; original maturity and female proportions remain separate in the canonical database."

fit_inputs <- inputs[
  !((inputs$type == "catch" & inputs$measure %in% c("total_numbers", "relative_standard_error", "proportion_at_age")) |
    (inputs$type == "index" & inputs$measure %in% c("total_numbers", "abundance_index_cv", "proportion_at_age")) |
    (inputs$type == "maturity" & inputs$measure == "maturity_at_age")),
]
fit_inputs <- rbind(fit_inputs, catch_rows, index_data, maturity)
surveys <- unique(index_data$survey)
obs <- database_to_tam_obs(
  assessment_id, fit_inputs, years = years, ages = ages,
  surveys = surveys, assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0
obs$pop_covariates$age_factor <- factor(obs$pop_covariates$age, levels = ages)
obs$index$q_key <- interaction(obs$index$survey, factor(obs$index$age, levels = ages), drop = TRUE, lex.order = TRUE)
zero_index_rows <- sum(obs$index$obs <= 0, na.rm = TRUE)
obs$index <- droplevels(obs$index[is.finite(obs$index$obs) & obs$index$obs > 0, , drop = FALSE])

## Comparisons and model ----

age_plus_group <- 15L
comparison_outputs <- source$outputs
n_source <- source$outputs[source$outputs$type == "population" & source$outputs$measure == "numbers_at_age", ]
weight <- inputs[inputs$type == "weight" & inputs$measure == "weight_at_age", ]
weight_key <- paste(weight$year, weight$age)
weight_at_age <- weight$value[match(paste(n_source$year, n_source$age), weight_key)]
if (anyNA(weight_at_age)) stop("A source January-1 weight is missing for a reported abundance age.")
biomass_age <- source$outputs[rep(1L, nrow(n_source)), , drop = FALSE]
biomass_age$type <- "biomass"
biomass_age$measure <- "biomass_at_age"
biomass_age$year <- n_source$year
biomass_age$age <- n_source$age
biomass_age$age_group <- n_source$age_group
biomass_age$value <- n_source$value * weight_at_age / 1000
biomass_age$se <- biomass_age$lwr <- biomass_age$upr <- NA_real_
biomass_age$unit <- "metric tons"
biomass_age$source_type <- "derived_source_output"
biomass_age$source_reference <- paste(unique(c(n_source$source_reference, weight$source_reference)), collapse = "; ")
biomass_age$notes <- "Derived January-1 biomass-at-age: official N-at-age times official January-1 Rivard weight-at-age; uncertainty is unavailable."
biomass_total <- aggregate(value ~ year, biomass_age, sum)
biomass_total_rows <- biomass_age[rep(1L, nrow(biomass_total)), , drop = FALSE]
biomass_total_rows$measure <- "total_biomass"
biomass_total_rows$year <- biomass_total$year
biomass_total_rows$age <- biomass_total_rows$age_group <- ""
biomass_total_rows$value <- biomass_total$value
biomass_total_rows$notes <- "Derived January-1 total biomass: sum of derived source biomass-at-age; uncertainty is unavailable."
comparison_outputs <- rbind(comparison_outputs, biomass_age, biomass_total_rows)
comparison_aggregates <- "biomass"
comparison_definitions <- list(
  N = list(status = "matched", definition = "Official January-1 abundance-at-age", reason = "The source and tinyAM report abundance at the beginning of the year; the process and likelihood formulations differ."),
  F = list(status = "matched", definition = "Official total instantaneous fishing mortality-at-age", reason = "The translated model pools two fleets and simplifies age selectivity and temporal variation."),
  M = list(status = "matched", definition = "Fixed age-specific source M input", reason = "This is an input consistency check, not an estimated-M comparison."),
  ssb = list(status = "approximate", definition = "Female spawning biomass using the translated effective mature-female fraction", reason = "tinyAM uses January-1 weights and the printed timing convention; the source reports distinct spawning weights and does not provide a complete within-year calculation for direct equivalence."),
  biomass = list(status = "matched", definition = "January-1 biomass calculated from the same source N-at-age and source January-1 weights", reason = "The official age-specific N and weight inputs are used; uncertainty is unavailable."),
  biomass_at_age = list(status = "matched", definition = "January-1 biomass-at-age calculated from the same source N-at-age and weights", reason = "The official age-specific N and weight inputs are used; uncertainty is unavailable."),
  recruitment = list(status = "matched", definition = "Same-calendar-year January-1 abundance at age 1", reason = "The source and tinyAM recruitment definitions align by age and year; the process formulations differ."),
  F_bar = list(status = "approximate", definition = "Population-weighted mean F over ages 3-8", reason = "The source's reported comparison is an unweighted mean over ages 3-8; tinyAM reports an abundance-weighted mean.")
)
background <- c(
  "### Atlantic striped bass: 2018 SAW-66 recommended single-stock SCA", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1982-2017 catch; indices are described from 1970, ages 1-15+. | Fit 1982-2017, ages 1-15, with age 15 as the plus group. | The report's numeric index tables begin in 1982; the earlier Maryland YOY and age-1 series are not numerically available in the published tables, and weights and catch begin in 1982. |",
  "| Catch | Chesapeake Bay and Coast fleet totals with age proportions. | Reconstruct each fleet's catch-at-age from its own total and proportions, then pool the fleets. | tinyAM uses one catch series and one observation SD; the source fits fleet-specific F and selectivity and a separate catch-composition likelihood. |",
  "| Indices | Eight age-composition surveys plus young-fish indices, with source CVs and composition effective sample sizes. | Convert aggregate indices and available age compositions to age-specific index values; link YOY to age 1 and age-1 indices to age 2. | tinyAM fits one lognormal likelihood by age and survey. Source CV weighting and multinomial composition likelihoods are not reproduced; the Table B10 “VA Age 1” label conflicts with Table B5.6's “NY Age-1” label. |",
  "| N and recruitment | Deterministic cohort survival with annual age-1 recruitment deviations. | Free initial abundance at older ages and IID recruitment deviations. | The source sum constraint, recruitment bias correction and early index data are not represented. |",
  "| F and M | Two fleets, annual F and time-blocked age selectivity; fixed age-specific M. | One total F-at-age AR1 process; use the fixed M schedule. | Fleet separation, selectivity functions and source observation weights are simplified. Fixed M is checked as an input, not as an estimated match. |",
  "| Biology and SSB | January-1 weights, fixed female maturity, female proportions and separate female SSB weights. | Use January-1 weights and female maturity multiplied by the female fraction in tinyAM's single maturity surface. | SSB is approximate because tinyAM has one weight surface; source spawning weights and the full spawning-time calculation are not represented. |",
  "| Outputs | N-at-age, total F-at-age, recruitment, total abundance and female SSB with natural-scale SE are reported. | Compare available N, F, recruitment, fixed M, approximate SSB and derived January-1 biomass. | Derived biomass has no uncertainty; source point estimates and available SSB uncertainty remain in the reference. |",
  "",
  "The 2018 historical benchmark is stored separately from the official 2024 update. The 2018 report recommends this single-stock model for management; its alternative two-stock migration model is not used. Catch-at-age and age-specific indices are analysis-only transformations of the original canonical totals and proportions. Zero survey-age observations are omitted from tinyAM's lognormal likelihood. The fit uses default starting values first; only a failed numerical check permits one retry initialized from official N and F outputs, without changing observations or model settings."
)

source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
model_settings <- list(
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "ar1", mu_form = ~0 + age_factor, mean_ages = 3:8),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~survey)
)
fit_model <- function(start_par = NULL) {
  args <- c(list(data = obs, years = years, ages = ages, silent = silent), model_settings)
  if (!is.null(start_par)) args$start_par <- start_par
  do.call(tinyAM::fit_tam, args)
}
source_starts <- function() {
  model_data <- do.call(tinyAM::prepare_tam, c(
    list(data = obs, years = years, ages = ages), model_settings
  ))
  starts <- tinyAM::make_par(model_data)
  source_n <- .translation_start_surface(source$outputs, "population", "numbers_at_age",
    years, ages, multiplier = 1 / 1000)
  f_rows <- source$outputs[source$outputs$type == "mortality" & source$outputs$measure == "fishing_mortality_at_age", ]
  f_rows$value <- as.numeric(f_rows$value)
  positive <- f_rows$value[is.finite(f_rows$value) & f_rows$value > 0]
  f_rows$value[!is.finite(f_rows$value) | f_rows$value <= 0] <- min(positive) / 2
  source_f <- .translation_start_surface(f_rows, "mortality", "fishing_mortality_at_age", years, ages)
  starts$log_r0 <- log(source_n[1L, 1L])
  if (length(starts$log_r)) starts$log_r[] <- log(source_n[-1L, 1L])
  if (length(starts$log_n0)) starts$log_n0[] <- log(source_n[1L, -1L])
  starts$log_mu_f[] <- colMeans(log(source_f))
  start_years <- intersect(rownames(starts$log_f), rownames(source_f))
  start_ages <- intersect(colnames(starts$log_f), colnames(source_f))
  starts$log_f[start_years, start_ages] <- log(source_f[start_years, start_ages, drop = FALSE])
  starts
}

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = source_starts,
    attempt_file = file.path("analysis/comp_assessments/results/cache", assessment_id,
      "initialization_attempts.rds")
  )
}
