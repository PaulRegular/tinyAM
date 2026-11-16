source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")

# This preliminary stock is absent from the SAIO release. Synthetic values
# on its documented age and year grids exercise the recipe independently.
assessment_id <- "argentina_merluccius_hubbsi_south_41s_2025"
age_grid <- expand.grid(year = 1986:2024, age = 1:7)
survey_grid <- expand.grid(year = 2011:2024, age = 1:7)
age_vector <- data.frame(year = NA_integer_, age = 1:7)
fixture_rows <- function(type, measure, basis, grid, value, unit,
                         survey = NA_character_, season = NA_character_) {
  data.frame(
    assessment_id, type, measure, basis, survey,
    fleet = NA_character_, sex = NA_character_, region = NA_character_, season,
    year = grid$year, year_basis = ifelse(is.na(grid$year), NA_character_, "calendar_year"),
    age = grid$age, value, unit, sampling_time = NA_real_,
    source_type = "synthetic", source_reference = "test fixture",
    transformation = NA_character_, notes = NA_character_
  )
}
record <- list(
  assessment = data.frame(assessment_id),
  inputs = rbind(
    fixture_rows("catch", "numbers_at_age", "numbers", age_grid, 100, "fish"),
    fixture_rows("index", "numbers_at_age", "numbers", survey_grid, 200, "fish",
      survey = "INIDEP global winter survey"),
    fixture_rows("weight", "weight_at_age", "kg_per_fish", age_vector, 1, "kg/fish",
      season = "1 January"),
    fixture_rows("weight", "weight_at_age", "kg_per_fish", age_vector, 2, "kg/fish",
      season = "1 July"),
    fixture_rows("maturity", "maturity_at_age", "proportion", age_vector, 0.5, "proportion"),
    fixture_rows("M", "natural_mortality_at_age", "rate", age_vector, 0.3, "per year")
  ),
  assumptions = NULL,
  outputs = data.frame(type = "biomass", measure = "spawning_biomass")
)
original <- record
stock <- .test_stock(record)

stopifnot(
  identical(record, original),
  identical(stock$years, 1986:2024),
  identical(stock$ages, 1:7),
  tinyAM::check_tam_data(stock$obs),
  nrow(stock$obs$catch) == 273L,
  nrow(stock$obs$index) == 98L,
  identical(unique(stock$obs$index$survey), "INIDEP global winter survey"),
  all(stock$obs$index$samp_time == 0.5),
  all(stock$obs$weight$obs == 2),
  all(stock$obs$pop_covariates$M_assumption == 0.3),
  identical(stock$comparison_outputs$measure, "SSB"),
  identical(stock$settings$N_settings$init, "exp"),
  identical(stock$settings$rec_settings$process, "iid"),
  identical(stock$settings$F_settings$process, "rw"),
  identical(stock$settings$M_settings$process, "off"),
  identical(deparse(stock$settings$M_settings$mu_supplied), "~M_assumption"),
  identical(deparse(stock$settings$index_settings$q_form), "~1"),
  identical(deparse(stock$settings$index_settings$sd_form), "~1"),
  is.null(stock$fit)
)

cat("Southern Argentine hake translation preparation and source-M mapping passed.\n")