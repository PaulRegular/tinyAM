source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

database <- read_database()
source <- read_assessment("ices_plaice_english_channel_7d_2026", database)
model <- .test_stock(source)

stopifnot(
  nrow(model$obs$catch) == 322L,
  nrow(model$obs$index) == 431L,
  setequal(unique(model$obs$index$survey), c("UK BTS", "FR GFS")),
  setequal(unique(model$obs$index$samp_time[model$obs$index$survey == "UK BTS"]), 0.625),
  setequal(unique(model$obs$index$samp_time[model$obs$index$survey == "FR GFS"]), 0.875),
  identical(model$settings$N_settings$process, "rw"),
  identical(model$settings$N_settings$init, "exp"),
  identical(model$settings$rec_settings$process, "rw"),
  identical(model$settings$F_settings$process, "cor_rw"),
  identical(model$settings$F_settings$mean_ages, 3:6),
  identical(model$settings$M_settings$process, "off"),
  length(unique(model$obs$index$q_key[model$obs$index$survey == "UK BTS"])) == 5L,
  length(unique(model$obs$index$q_key[model$obs$index$survey == "FR GFS"])) == 5L,
  !is.null(model$fit_call),
  is.null(model$warm_call),
  model$comparison_aggregates == "ssb",
  model$comparison_definitions$F_bar$status == "non_equivalent",
  model$comparison_definitions$ssb$status == "matched"
)

cat("Eastern English Channel plaice translation prepares without a fit or source-based starts.\n")
