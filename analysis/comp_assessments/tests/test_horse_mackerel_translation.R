pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/run_assessment.R")

source_record <- read_assessment("ices_horse_mackerel_9a_2024")
translation <- .test_stock(
  source_record,
  script = "analysis/comp_assessments/scripts/translation/stocks/ices_horse_mackerel_9a_2024.R"
)

stopifnot(
  identical(as.integer(sort(unique(translation$obs$catch$year))), 1992:2023),
  identical(as.integer(sort(unique(translation$obs$catch$age))), 0:11),
  nrow(translation$obs$catch) == 32L * 12L,
  all(translation$obs$catch$year %in% 1992:2023),
  all(translation$obs$index$age %in% 1:11),
  all(translation$obs$index$obs > 0),
  all(translation$obs$index$relative_sd == 0.30),
  all(translation$obs$index$samp_time == 0.8),
  identical(as.integer(sort(unique(translation$obs$index$year))), setdiff(1992:2023, c(2012, 2019, 2020))),
  identical(translation$settings$N_settings$init, "exp"),
  identical(translation$settings$rec_settings$process, "iid"),
  identical(translation$settings$F_settings$process, "rw"),
  is.null(translation$fit_call$start_par),
  identical(translation$comparison_definitions$ssb$status, "non_equivalent"),
  identical(translation$comparison_definitions$F_bar$status, "non_equivalent")
)

cat("Southern horse mackerel translation preparation and source limitations passed.\n")
