source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
record <- read_assessment("afsc_pollock_ai_2024")
before <- record
stock <- .test_stock(record)
stopifnot(identical(record, before), is.null(stock$fit),
  identical(stock$years, 1978:2018), identical(stock$ages, 1:15),
  all(stock$obs$maturity$spawn_time == 2 / 12),
  all(stock$obs$index$samp_time == .5),
  all(stock$obs$index$q_block[stock$obs$index$age >= 8] == "8"),
  !"M_assumption" %in% names(stock$obs$pop_covariates),
  length(tinyAM::make_par(stock$dat)$mu_m) == 1,
  identical(stock$dat$rec$curve$type, "bh"), stock$dat$rec$curve$lag == 1L,
  !"start_par" %in% names(stock$fit_call),
  stock$obs$maturity$obs[stock$obs$maturity$year == 1978 & stock$obs$maturity$age == 4] == .269 * .5,
  stock$obs$weight$obs[stock$obs$weight$year == 1978 & stock$obs$weight$age == 1] == .080)
cat("AI pollock preparation, source preservation, biology and initialization passed.\n")
