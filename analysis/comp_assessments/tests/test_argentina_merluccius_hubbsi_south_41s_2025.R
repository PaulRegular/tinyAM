source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/read_database.R")

assessment_id <- "argentina_merluccius_hubbsi_south_41s_2025"
database <- read_database()
record <- read_assessment(assessment_id, database)
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
  all(stock$obs$pop_covariates$M_assumption == 0.3),
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