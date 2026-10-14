test_that("missing and zero observations add no nuisance states or likelihood terms", {
  data <- cod_obs
  rows <- which(data$catch$year == 2002 & data$catch$age %in% 2:3)
  data$catch$obs[rows] <- c(NA, 0)
  dat <- make_test_dat(data = data, years = 2000:2005, ages = 2:6)
  par <- make_par(dat)
  expect_false("missing" %in% names(par))
  expect_false(any(c("fill_missing_map", "any_fill_missing") %in% names(dat)))
  expect_true(any(dat$is_missing))
  no_obs <- dat
  no_obs$is_observed[] <- FALSE
  no_obs$observed <- numeric()
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expected <- -sum(stats::dnorm(dat$observed, report$log_pred[dat$is_observed],
                                report$sd_obs[dat$is_observed], log = TRUE))
  expect_equal(nll_fun(par, dat) - nll_fun(par, no_obs), expected)
  expect_true(all(is.finite(report$log_pred)))
  expect_true(all(is.finite(report$N)))
  set.seed(182)
  simulated <- nll_fun(par, dat, simulate = TRUE)
  expect_true(all(is.na(simulated$log_obs[dat$is_missing])))
  expect_true(all(is.finite(simulated$log_obs[dat$is_observed])))
  expect_identical(data$catch$obs[rows], c(NA_real_, 0))
})

test_that("hindcast RMSE uses real observations, not missing or simulated forecast observations", {
  d <- data.frame(year = c(2021, 2021), age = 2, fold = c(2020, 2021),
                  is_proj = c(TRUE, FALSE), obs = c(NA, 100), pred = c(110, 99))
  expect_equal(compute_hindcast_rmse(d), log(1.1))
  expect_equal(compute_hindcast_rmse(d, log = FALSE), 10)
  d$obs[1] <- 999
  expect_equal(compute_hindcast_rmse(d), log(1.1))
  d$obs[] <- NA
  expect_true(is.nan(compute_hindcast_rmse(d)))
})
