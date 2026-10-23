source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))

database <- read_database()
source_data <- read_assessment("nefsc_winter_flounder_gb_2022_vpa", database)
original_inputs <- source_data$inputs
translated <- .test_stock(source_data)
obs <- translated$obs
index <- obs$index

tinyAM::check_obs(obs)

stopifnot(
  identical(translated$years, 1982:2021),
  identical(translated$ages, 1:7),
  translated$age_plus_group == 7L,
  translated$zero_index_cells == 99L,
  nrow(index) == 741L,
  all(index$obs > 0),
  identical(as.integer(sort(unique(index$year))), 1982:2021),
  identical(as.integer(sort(unique(index$age))), 1:7),
  setequal(as.character(unique(index$survey)), c("DFO spring", "NEFSC fall", "NEFSC spring")),
  all(index$samp_time == 0),
  all(as.character(index$q_key) == as.character(interaction(index$survey, index$age, drop = TRUE))),
  identical(translated$settings$N_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$F_settings$process, "ar1"),
  identical(translated$settings$M_settings$process, "off"),
  identical(translated$comparison_definitions$ssb$status, "approximate"),
  identical(translated$comparison_aggregates, "ssb"),
  identical(source_data$inputs, original_inputs)
)

cat("Georges Bank winter flounder translation checks passed.\n")
