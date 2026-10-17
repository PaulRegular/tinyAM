source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

testthat::test_that("translation units describe converted values and preserve source records", {
  database <- read_database()
  for (id in c("afsc_pollock_ebs_2024", "ices_cod_northeast_arctic_2026")) {
    source <- read_assessment(id, database)
    original <- source
    translated <- .test_stock(source)
    testthat::expect_identical(source, original)
    testthat::expect_true(all(translated$obs$catch$unit == "thousand fish"))
    testthat::expect_true(all(translated$obs$weight$unit == "kg"))
    testthat::expect_true(all(vapply(split(translated$obs$index$unit,
      translated$obs$index$survey), function(u) length(unique(u)) == 1L, logical(1))))
    testthat::expect_true(tinyAM::check_tam_data(translated$obs))
  }
  source <- read_assessment("dfo_cod_2j3kl_2025", database)
  stock <- .test_stock(source)
  native <- source$inputs[source$inputs$type == "index" &
    source$inputs$measure == "numbers_at_age" & source$inputs$unit == "mean numbers per tow", ]
  for (survey in unique(native$survey)) {
    rows <- stock$obs$index$survey == survey
    if (any(rows)) testthat::expect_true(all(stock$obs$index$unit[rows] == "thousand fish per tow"))
  }
})

testthat::test_that("retained source-checked recipes use latent initialization", {
  database <- read_database()
  for (id in c("ices_bluewhiting_northeast_atlantic_2026",
    "ices_haddock_north_sea_2026", "ices_whiting_north_sea_2026")) {
    stock <- .test_stock(read_assessment(id, database))
    testthat::expect_identical(stock$dat$N_settings$init, "latent")
    testthat::expect_length(tinyAM::make_par(stock$dat)$log_n0, length(stock$dat$ages) - 1L)
  }
  stock <- .test_stock(read_assessment("ices_cod_northeast_arctic_2026", database))
  testthat::expect_identical(stock$dat$N_settings$init, "free")
})
