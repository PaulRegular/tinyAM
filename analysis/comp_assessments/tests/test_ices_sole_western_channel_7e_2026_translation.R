source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
assessment_id <- "ices_sole_western_channel_7e_2026"
assessment <- read_assessment(assessment_id, database)
trans <- .test_stock(assessment)
translated_index <- trans$obs$index
stopifnot(
  identical(trans$years, 1969:2025),
  identical(trans$ages, 2:12),
  identical(trans$age_plus_group, 12L),
  all(grepl("fish", trans$obs$catch$unit)),
  all(translated_index$obs > 0 & is.finite(translated_index$obs)),
  setequal(unique(translated_index$survey), c("UK-CBT-late", "UK-COT", "Q1SWBeam-nonoffset", "FSP-UK")),
  max(translated_index$age[translated_index$survey == "UK-CBT-late"]) == 11L,
  min(translated_index$age[translated_index$survey == "Q1SWBeam-nonoffset"]) == 2L,
  max(translated_index$age[translated_index$survey == "FSP-UK"]) == 11L,
  nrow(translated_index) == 868L,
  !any(translated_index$year[translated_index$survey == "Q1SWBeam-nonoffset"] == 2022L),
  min(translated_index$year[translated_index$survey == "FSP-UK"]) == 2004L,
  max(translated_index$year[translated_index$survey == "UK-COT"]) == 2015L,
  all(translated_index$samp_time[translated_index$survey %in% c("UK-CBT-late", "UK-COT")] == 0.5),
  all(translated_index$samp_time[translated_index$survey == "Q1SWBeam-nonoffset"] == 0.175),
  all(translated_index$samp_time[translated_index$survey == "FSP-UK"] == 0.725),
  all(translated_index$q_key == interaction(translated_index$survey, pmin(translated_index$age, 7L), drop = TRUE, sep = ".")),
  all(trans$obs$pop_covariates$M_assumption == 0.1),
  identical(trans$settings$N_settings$process, "iid"),
  identical(trans$settings$N_settings$init, "exp"),
  identical(trans$settings$rec_settings$process, "rw"),
  identical(trans$settings$F_settings$process, "ar1"),
  identical(trans$settings$F_settings$mean_ages, 3:9),
  identical(trans$settings$M_settings$process, "off"),
  is.null(trans$warm_call),
  !is.null(trans$fit_call)
)
stopifnot(
  identical(trans$comparison_definitions$M$status, "approximate"),
  identical(trans$comparison_definitions$M_bar$status, "approximate"),
  grepl("input check", trans$comparison_definitions$M$reason, ignore.case = TRUE)
)

cat("ICES Western Channel sole tinyAM preparation passed.\n")