pkgload::load_all(quiet = TRUE)
source("analysis/formula_expansion/sd_validation.R")

testthat::test_that("replicated SD recovery exercises the integrated scale likelihood", {
  for (type in c("iid", "rw", "ar1")) {
    result <- sd_simple_recovery(type, 30L, 240101L)
    testthat::expect_true(is.finite(result$diagnostics$objective))
    testthat::expect_true(all(is.finite(result$curve$estimate) & result$curve$estimate > 0))
    testthat::expect_lt(result$log_sd_rmse, .3)
    testthat::expect_equal(result$parameters$truth[1], .35)
  }
})

testthat::test_that("the tail experiment can distinguish a common SD from an age curve", {
  common <- sd_tail_recovery("common", 380001L)
  quadratic <- sd_tail_recovery("quadratic", 380001L)
  testthat::expect_lt(quadratic$log_sd_rmse, common$log_sd_rmse / 2)
  testthat::expect_equal(common$curve$truth, quadratic$curve$truth)
  set.seed(380001L)
  truth <- rep(quadratic$curve$truth, 30L)
  y <- rnorm(length(truth), 0, truth)
  predicted <- rep(quadratic$curve$estimate, 30L)
  # The intercept score is zero at the MLE, not at sdreport's last perturbation.
  testthat::expect_lt(abs(sum(1 - (y / predicted)^2)), .01)
})
