source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
database <- read_database(tier = "preliminary")

testthat::test_that("the first expansion recipes preserve source values and never fit during preparation", {
  ids <- c("ices_ling_faroe_5b_2026", "ices_atlantic_wolffish_iceland_2026")
  testthat::local_mocked_bindings(fit_tam = function(...) {
    stop("Preparation must not fit or construct a source-start retry")
  }, .package = "tinyAM")
  for (id in ids) {
    record <- read_assessment(id, database)
    original <- record
    stock <- .test_stock(record)
    testthat::expect_null(stock$fit)
    testthat::expect_identical(record, original)
    testthat::expect_identical(stock$settings$N_settings$init, "free")
    testthat::expect_identical(stock$settings$F_settings$process, "rw")
    testthat::expect_true(all(abs(stock$obs$pop_covariates$M_assumption - 0.15) < 1e-12))
    testthat::expect_true(all(stock$obs$index$unit == "native survey index"))
    observed <- stock$obs$index
    source_index <- record$inputs[record$inputs$type == "index" &
      record$inputs$measure %in% c("numbers_at_age", "relative_abundance_index") &
      record$inputs$year %in% stock$years, ]
    key <- function(x) paste(x$survey, x$year, x$age)
    testthat::expect_equal(observed$obs, source_index$value[match(key(observed), key(source_index))])
    if (id == "ices_ling_faroe_5b_2026") {
      testthat::expect_equal(stock$years, 1996:2025)
      testthat::expect_equal(stock$ages, 3:12)
      spring <- grepl("spring", observed$survey)
      testthat::expect_equal(observed$samp_time, ifelse(spring, 2 / 12, 7.5 / 12))
    } else {
      testthat::expect_equal(stock$years, 1998:2025)
      testthat::expect_equal(stock$ages, 4:16)
      testthat::expect_false(any(stock$obs$catch$year %in% 1981:1997))
      raw_weight <- record$inputs[record$inputs$type == "weight" &
        record$inputs$year %in% stock$years, ]
      index <- match(paste(stock$obs$weight$year, stock$obs$weight$age),
        paste(raw_weight$year, raw_weight$age))
      testthat::expect_equal(stock$obs$weight$obs, raw_weight$value[index] / 1000)
    }
  }
})
