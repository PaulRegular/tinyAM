source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

testthat::test_that("translations use thousands, kg and matching count-rate scales", {
  db <- read_database()
  for (id in c("ices_whiting_north_sea_2026", "afsc_pollock_goa_2024", "dfo_cod_2j3kl_2025")) {
    source <- read_assessment(id, db)
    original <- source
    stock <- .test_stock(source)
    testthat::expect_identical(source, original)
    testthat::expect_equal(unique(stock$obs$catch$unit), "thousand fish")
    testthat::expect_equal(unique(stock$obs$weight$unit), "kg")
    testthat::expect_equal(tinyAM:::.data_units(stock$obs)$mass, "tonnes")
    for (u in unique(stock$obs$index$unit)) {
      count <- tinyAM:::.count_unit(u)
      if (is.finite(count$scale)) testthat::expect_equal(count$scale, 1000)
    }
  }
})

testthat::test_that("report conversions preserve log SEs and natural intervals", {
  x <- data.frame(est = 2000, lwr = 1000, upr = 3000, se = .2,
    se_scale = "log", unit = "fish")
  y <- tinyAM:::.convert_report_units(x, "thousand fish", "numbers")
  testthat::expect_equal(y$est, 2)
  testthat::expect_equal(y$lwr, 1)
  testthat::expect_equal(y$upr, 3)
  testthat::expect_equal(y$se, .2)
  x$se <- 100; x$se_scale <- "reported"
  testthat::expect_equal(tinyAM:::.convert_report_units(x, "thousand fish", "numbers")$se, .1)
})
