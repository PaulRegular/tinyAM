root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "ices_herring_west_scotland_ireland_2020"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
translated <- run_assessment(
  assessment_id, database = database, fit = FALSE, silent = TRUE
)
obs <- translated$obs

source_catch <- source_data$inputs[
  source_data$inputs$type == "catch" &
    source_data$inputs$measure == "numbers_at_age",
  , drop = FALSE
]
source_catch_key <- paste(source_catch$year, source_catch$age, sep = "/")
source_catch_total <- tapply(
  as.numeric(source_catch$value), source_catch_key, sum
)
translated_key <- paste(obs$catch$year, obs$catch$age, sep = "/")
expected_catch <- as.numeric(source_catch_total[translated_key])

source_F <- source_data$outputs[
  source_data$outputs$type == "mortality" &
    source_data$outputs$measure == "fishing_mortality_at_age",
  , drop = FALSE
]
source_F_key <- paste(source_F$year, source_F$age, sep = "/")
source_F_total <- tapply(as.numeric(source_F$value), source_F_key, sum)
comparison_F <- translated$comparison$comparison_outputs[
  translated$comparison$comparison_outputs$source_type == "derived_source_output" &
    translated$comparison$comparison_outputs$measure == "fishing_mortality_at_age",
  , drop = FALSE
]
comparison_F_key <- paste(comparison_F$year, comparison_F$age, sep = "/")

stopifnot(
  identical(translated$diagnostics$status, "not_fitted"),
  is.null(translated$fit),
  inherits(obs, "tam_data"),
  tinyAM::check_obs(obs),
  sum(is.finite(obs$catch$obs)) == 567L,
  isTRUE(all.equal(obs$catch$obs, expected_catch, tolerance = 1e-12)),
  identical(unique(source_catch$fleet), c("Catch 6aN", "Catch 6aS/7bc")),
  setequal(unique(obs$index$survey), c("MS_HERAS", "IBTS_Q1", "IBTS_Q4")),
  sum(is.finite(obs$index$obs)) == 640L,
  length(unique(obs$index$q_key[obs$index$survey == "MS_HERAS"])) == 8L,
  length(unique(obs$index$q_key[obs$index$survey == "IBTS_Q1"])) == 7L,
  length(unique(obs$index$q_key[obs$index$survey == "IBTS_Q4"])) == 7L,
  length(unique(obs$index$sd_block[obs$index$survey == "MS_HERAS"])) == 4L,
  length(unique(obs$index$sd_block[obs$index$survey == "IBTS_Q1"])) == 4L,
  length(unique(obs$index$sd_block[obs$index$survey == "IBTS_Q4"])) == 3L,
  all(obs$index$unit == "native survey index"),
  all(obs$maturity$spawn_time == 0.67),
  all(is.finite(obs$pop_covariates$M_assumption)),
  nrow(comparison_F) == 567L,
  isTRUE(all.equal(
    comparison_F$value,
    as.numeric(source_F_total[comparison_F_key]),
    tolerance = 1e-12
  )),
  identical(translated$comparison$age_plus_group, 9L),
  any(grepl("independent survey-age residuals", translated$background, fixed = TRUE))
)

cat("West of Scotland/Ireland herring tinyAM translation preparation passed.\n")
