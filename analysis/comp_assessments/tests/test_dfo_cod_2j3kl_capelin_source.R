testthat::test_that("Northern cod Capelin covariate preserves the digitized 3L index", {
  root <- normalizePath(getwd())
  while (!file.exists(file.path(root, "analysis", "comp_assessments", "database", "inputs.csv"))) {
    parent <- dirname(root)
    if (identical(parent, root)) stop("Could not locate the assessment database.")
    root <- parent
  }
  inputs <- read.csv(
    file.path(root, "analysis", "comp_assessments", "database", "inputs.csv"),
    stringsAsFactors = FALSE, na.strings = "", check.names = FALSE
  )
  capelin <- inputs[
    inputs$assessment_id == "dfo_cod_2j3kl_2025" &
      inputs$measure == "capelin_biomass_index", , drop = FALSE
  ]
  years <- c(1985:1992, 1996, 1999:2005, 2007:2015, 2017:2019, 2022:2024)
  log_biomass_kt <- c(
    8.13, 8.20, 7.83, 8.34, 8.20, 8.65, 4.92, 4.92, 3.85,
    5.22, 4.45, 4.80, 4.33, 4.66, 4.71, 4.54,
    5.50, 5.52, 5.50, 3.08, 5.29, 5.15, 6.48, 6.64, 5.90,
    4.87, 5.64, 5.62, 5.55, 5.81, 6.46
  )

  testthat::expect_equal(nrow(capelin), length(years))
  testthat::expect_equal(as.integer(capelin$year), years)
  testthat::expect_true(all(is.finite(capelin$value) & capelin$value > 0))
  testthat::expect_equal(log(capelin$value), log_biomass_kt, tolerance = 1e-8)
  testthat::expect_true(all(capelin$type == "covariate"))
  testthat::expect_true(all(capelin$basis == "native_covariate"))
  testthat::expect_true(all(capelin$survey == "3L spring acoustic survey"))
  testthat::expect_true(all(capelin$region == "NAFO Division 3L"))
  testthat::expect_true(all(capelin$season == "spring"))
  testthat::expect_true(all(capelin$unit == "kt"))
  testthat::expect_true(all(capelin$source_type == "digitized"))
  testthat::expect_true(all(grepl("Figure 32", capelin$source_reference)))
  testthat::expect_true(all(grepl("round the plotted log value to 0.01", capelin$transformation)))
  testthat::expect_true(all(grepl("not the original", capelin$notes)))
})
