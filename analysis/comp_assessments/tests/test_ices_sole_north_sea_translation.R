source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_sole_north_sea_2025"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
source_inputs <- source_data$inputs
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat

tinyAM::check_obs(obs)

source_catch <- source_inputs[
  source_inputs$type == "catch" &
    source_inputs$measure == "numbers_at_age",
  , drop = FALSE
]
expected_catch <- stats::aggregate(value ~ year + age, source_catch, sum)
expected_catch_index <- match(
  paste(obs$catch$year, obs$catch$age),
  paste(expected_catch$year, expected_catch$age)
)

source_index <- source_inputs[
  source_inputs$type == "index" &
    source_inputs$measure == "numbers_at_age" &
    source_inputs$survey %in% c("BTS Q3", "COAST (SNS + DYFS)") &
    source_inputs$age %in% 1:10,
  , drop = FALSE
]
expected_index <- match(
  paste(obs$index$survey, obs$index$year, obs$index$age),
  paste(source_index$survey, source_index$year, source_index$age)
)

expected_m <- c(.754, .465, .368, .321, .295, .278, .267, .260, .254, .248)
actual_m <- vapply(1:10, function(age) {
  unique(obs$pop_covariates$M_assumption[obs$pop_covariates$age == age])
}, numeric(1))
readiness <- read.csv(
  file.path(root, "results", "observation_readiness.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
sole_readiness <- readiness[
  readiness$assessment_id == assessment_id, , drop = FALSE
]

stopifnot(
  identical(translated$years, 1957:2024),
  identical(translated$ages, 1:10),
  identical(translated$age_plus_group, 10L),
  nrow(obs$catch) == 68L * 10L,
  nrow(obs$index) == 35L * 15L,
  nrow(obs$weight) == 68L * 10L,
  nrow(obs$maturity) == 68L * 10L,
  all(!is.na(expected_catch_index)),
  isTRUE(all.equal(obs$catch$obs, expected_catch$value[expected_catch_index])),
  all(!is.na(expected_index)),
  isTRUE(all.equal(obs$index$obs, source_index$value[expected_index])),
  setequal(unique(obs$index$survey), c("BTS Q3", "COAST (SNS + DYFS)")),
  all(obs$index$samp_time == 0.625),
  all(obs$index$unit == "native survey index (unit unresolved)"),
  length(unique(obs$index$q_key)) == 10L,
  all(as.character(obs$index$q_block[obs$index$age >= 5]) == "5+"),
  all(obs$maturity$spawn_time == 0),
  isTRUE(all.equal(actual_m, expected_m)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$rec$type, "rw"),
  is.null(dat$rec$curve),
  identical(dat$F_settings$process, "rw"),
  identical(dat$F_settings$mean_ages, 2:6),
  identical(dat$M_settings$process, "off"),
  is.null(translated$start_par),
  !"relative_biomass_index" %in% obs$index$measure,
  identical(source_data$inputs, source_inputs),
  nrow(sole_readiness) == 1L,
  isTRUE(sole_readiness$weight_full_year_age_grid),
  isTRUE(sole_readiness$catch_weight_full_year_age_grid),
  isTRUE(sole_readiness$maturity_full_year_age_grid),
  isTRUE(sole_readiness$M_numerical_values_available),
  !isTRUE(sole_readiness$M_estimated_in_source)
)

cat("North Sea sole translation structure passed.\n")
