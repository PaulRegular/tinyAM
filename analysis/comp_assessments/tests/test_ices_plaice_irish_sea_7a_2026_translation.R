source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

database <- read_database()
source_record <- read_assessment("ices_plaice_irish_sea_7a_2026", database)
translated <- .test_stock(source_record)

stopifnot(
  identical(translated$years, 1981:2025),
  identical(translated$ages, 1:8),
  identical(translated$age_plus_group, 8L),
  nrow(translated$obs$catch) == 360L,
  sum(is.na(translated$obs$catch$obs)) == 2L,
  all(translated$obs$catch$unit == "thousand fish"),
  nrow(translated$obs$index) == 224L,
  identical(unique(translated$obs$index$survey), "UK BTS-Q3"),
  identical(unique(translated$obs$index$samp_time), 0.8),
  length(unique(translated$obs$index$q_key)) == 7L,
  identical(translated$settings$N_settings$process, "rw"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$rec_settings$process, "rw"),
  identical(translated$settings$F_settings$process, "cor_rw"),
  identical(translated$settings$M_settings$process, "off"),
  identical(paste(deparse(translated$settings$index_settings$q_form), collapse = ""), "~0 + q_key"),
  translated$comparison_definitions$F_bar$status == "non_equivalent",
  translated$comparison_definitions$ssb$status == "matched",
  !is.null(translated$fit_call),
  is.null(translated$warm_call)
)

cat("Irish Sea plaice translation prepares with the missing catch cells and the documented model settings.\n")
