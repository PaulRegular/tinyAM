pkgload::load_all(quiet = TRUE, helpers = FALSE)
source("analysis/comp_assessments/R/source_start_parameters.R")

testthat::test_that("source starts align year-age names and leave missing values at defaults", {
  dat <- list(years = 2000:2002, ages = 2:4, is_proj = rep(FALSE, 3))
  defaults <- list(log_r0 = 0, log_n0 = c(`3` = 0, `4` = 0),
    log_r = c(`2001` = 0, `2002` = 0),
    log_n = matrix(0, 2, 2, dimnames = list(c("2001", "2002"), c("3", "4"))),
    log_f = matrix(0, 3, 3, dimnames = list(as.character(2000:2002), as.character(2:4))),
    log_sd_n = 0, log_q = c(survey = 0))
  testthat::local_mocked_bindings(make_par = function(dat) defaults, .package = "tinyAM")
  n <- expand.grid(year = 2000:2002, age = 2:4)
  n$type <- "population"
  n$measure <- "numbers_at_age"
  n$value <- n$age
  n$unit <- "million fish"
  f <- n
  f$type <- "mortality"
  f$measure <- "fishing_mortality_at_age"
  f$value <- 0.2
  f$unit <- "per year"
  f$value[f$year == 2000 & f$age == 3] <- 0
  n <- n[!(n$year == 2002 & n$age == 4), ]
  source <- list(outputs = rbind(n[nrow(n):1, ], f[nrow(f):1, ]))
  original <- source
  starts <- source_start_parameters(source, dat)
  testthat::expect_equal(starts$log_r0, log(2000))
  testthat::expect_equal(unname(starts$log_n0), log(c(3000, 4000)))
  testthat::expect_equal(unname(starts$log_r), rep(log(2000), 2))
  testthat::expect_equal(starts$log_n["2001", "4"], log(4000))
  testthat::expect_identical(starts$log_n["2002", "4"], 0)
  testthat::expect_identical(starts$log_f["2000", "3"], 0)
  testthat::expect_equal(starts$log_f["2002", "4"], log(0.2))
  testthat::expect_identical(starts$log_q, defaults$log_q)
  testthat::expect_identical(starts$log_sd_n, defaults$log_sd_n)
  testthat::expect_identical(source, original)
  source$outputs <- rbind(source$outputs, source$outputs[1, ])
  testthat::expect_error(source_start_parameters(source, dat), "unique year-age")
})
