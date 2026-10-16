source("analysis/comp_assessments/R/fit_default_first.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

testthat::test_that("default starts run first and successful fits do not construct starts", {
  calls <- list()
  model <- function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    structure(list(passed = TRUE), class = "tam_fit")
  }
  testthat::local_mocked_bindings(check_tam = function(fit) {
    list(is_converged = fit$passed, max_gradient = 0)
  }, .package = "tinyAM")
  fit <- fit_default_first(model, function() stop("Starts must remain lazy"))
  testthat::expect_length(calls, 1L)
  testthat::expect_length(calls[[1L]], 0L)
  testthat::expect_identical(attr(fit, "initialization_attempts")$initialization, "default")
})

testthat::test_that("failed numerical checks allow exactly one initialization-only retry", {
  calls <- list()
  starts_built <- 0L
  model <- function(...) {
    arguments <- list(...)
    calls[[length(calls) + 1L]] <<- arguments
    structure(list(passed = length(arguments) > 0L), class = "tam_fit")
  }
  testthat::local_mocked_bindings(check_tam = function(fit) {
    list(is_converged = fit$passed, max_gradient = if (fit$passed) 0 else 1)
  }, .package = "tinyAM")
  starts <- list(log_r0 = 10)
  fit <- fit_default_first(model, function() {
    starts_built <<- starts_built + 1L
    starts
  })
  testthat::expect_identical(starts_built, 1L)
  testthat::expect_identical(calls, list(list(), list(start_par = starts)))
  attempts <- attr(fit, "initialization_attempts")
  testthat::expect_identical(attempts$is_converged, c(FALSE, TRUE))
  testthat::expect_equal(attempts$max_abs_gradient, c(1, 0))

  calls <- list()
  testthat::local_mocked_bindings(check_tam = function(fit) {
    list(is_converged = FALSE, max_gradient = 1)
  }, .package = "tinyAM")
  fit <- fit_default_first(model, function() starts)
  testthat::expect_length(calls, 2L)
  testthat::expect_true(all(!attr(fit, "initialization_attempts")$is_converged))
})

testthat::test_that("BSAI POP preparation leaves the source untouched and runs no fit", {
  record <- read_assessment("afsc_pop_bsai_2024")
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) {
    stop("Fitting-disabled preparation must not fit")
  }, .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_null(stock$fit)
  testthat::expect_equal(stock$ages, 3:40)
  testthat::expect_equal(stock$years, 1960:2024)
  testthat::expect_true(all(stock$obs$maturity$spawn_time == 0.25))
  testthat::expect_equal(stock$obs$weight$obs[stock$obs$weight$age == 3], rep(0.094, 65))
  testthat::expect_true(all(stock$obs$index$samp_time == 0.5))
  testthat::expect_true(all(stock$obs$maturity$obs <= 0.5))
  testthat::expect_true(all(stock$settings$M_settings$process == "off"))
  testthat::expect_equal(stock$settings$M_settings$mu_form, ~ 1, ignore_attr = TRUE)
  par <- tinyAM::make_par(stock$dat)
  testthat::expect_length(par$mu_m, 1L)
  ref <- database_to_tam_ref("afsc_pop_bsai_2024", record$outputs,
    obs = stock$obs, years = stock$years, ages = stock$ages,
    age_plus_group = 40, assumptions = record$assumptions)
  testthat::expect_equal(sum(ref$pop$M$age == 40), 65L)
  testthat::expect_true(all(ref$pop$M$est[ref$pop$M$age == 40] == 0.0505762))
})
