pkgload::load_all(quiet = TRUE)
source("analysis/formula_expansion/recruitment_cross_validation.R")

testthat::test_that("all five recruitment truths have explicit self- and cross-fits", {
  for (truth in rec_cross_cases) {
    result <- rec_cross_dataset("isolated", truth, 1L)
    attempts <- do.call(rbind, lapply(result$results, `[[`, "attempt"))
    testthat::expect_setequal(attempts$fitted, rec_cross_cases)
    testthat::expect_true(all(is.finite(attempts$objective)))
    testthat::expect_true(result$results[[truth]]$attempt$success)
  }
})

testthat::test_that("a RW truth has no curve parameters to recover", {
  simulated <- rec_cross_simulation("rw", "isolated", 12L)
  testthat::expect_length(rec_cross_parameters(simulated, "bh_iid"), 0L)
  testthat::expect_length(rec_cross_parameters(simulated, "ricker_ar1"), 0L)
  testthat::expect_equal(rec_cross_parameters(simulated, "rw")$log_sd_r, log(.35))
  d <- do.call(prepare_tam, c(list(data = simulated$obs), rec_settings(~ rw(year))))
  expected <- -sum(dnorm(diff(simulated$log_R), 0, .35, log = TRUE))
  testthat::expect_equal(tinyAM:::.rec_nll(simulated$log_R, rep(0, length(d$years)), simulated$truth, d), expected)
})

testthat::test_that("noise and parent contrast are separate generating assumptions", {
  small <- rec_cross_simulation("bh_iid", "isolated", 72L, sigma = .15, contrast = "wider")
  large <- rec_cross_simulation("bh_iid", "isolated", 72L, sigma = .6, contrast = "wider")
  testthat::expect_equal(small$S, large$S)
  testthat::expect_equal(small$truth$log_sr_alpha, large$truth$log_sr_alpha)
  testthat::expect_equal(small$truth$log_sr_beta, large$truth$log_sr_beta)
  testthat::expect_equal(small$truth$log_sd_r, log(.15))
  testthat::expect_equal(large$truth$log_sd_r, log(.6))
  full <- rec_cross_simulation("bh_iid", "full", 72L, contrast = "wider")
  testthat::expect_true("f_pressure" %in% names(full$obs$catch))
  testthat::expect_true(all(is.finite(full$log_R)))
})
