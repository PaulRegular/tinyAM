
testthat::skip_if_not_installed("RTMB")
testthat::skip_if_not(exists("northern_cod_data"), "northern_cod_data not available")

set.seed(1)

YEARS <- 1983:2024
AGES  <- 2:14

test_fit <- function(inputs = northern_cod_data,
                     years = YEARS,
                     ages  = AGES,
                     N_settings = list(process = "iid", init_N0 = FALSE),
                     F_settings = list(process = "approx_rw",  mu_form = NULL),
                     M_settings = list(process = "off", mu_form = NULL, assumption = ~I(0.3)),
                     obs_settings = list(q_form = ~ q_block, sd_form = ~ sd_obs_block),
                     silent = TRUE) {
  args <- mget(ls())
  do.call(fit_tam, args)
}

default_fit <- test_fit()

test_that("fit_tam runs on a cod dataset and returns expected structure", {
  fit <- default_fit

  # Structure
  expect_type(fit, "list")
  expect_named(
    fit,
    c("call","dat","obj","opt","rep","sdrep"),
    ignore.order = TRUE
  )

  # Optimizer status
  expect_true(is.finite(fit$opt$objective))
  expect_equal(round(fit$opt$objective, 4), 989.6083)
  expect_true(is.list(fit$rep))
  expect_s3_class(fit$sdrep, "sdreport")

  # Reported matrices have expected dims
  expect_equal(dim(fit$rep$F), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$M), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$Z), c(length(YEARS), length(AGES)))
  expect_equal(length(fit$rep$ssb), length(YEARS))
})

test_that("fit_tam emits warning if random effects far exceed observations", {
  expect_warning(
    test_fit(
      N_settings = list(process = "iid", init_N0 = TRUE),
      F_settings = list(process = "iid"),
      M_settings = list(process = "iid", mu_form = NULL, assumption = ~I(0.3))
    ),
    regexp = "Number of random effects .* exceed 1.5 times the number of observations",
    fixed  = FALSE
  )
})

test_that("fit_retro produces a list of fits with decreasing end years", {
  retro_fit <- fit_retro(default_fit, folds = 2)
  expect_true("retro" %in% names(retro_fit))
  expect_type(retro_fit$retro, "list")
  expect_equal(length(retro_fit$retro), 3) # max year, and two peels
})

test_that("sim_tam returns simulated observations and (optionally) re-computed reports", {
  # Simulate obs only (quick)
  rep_obs <- sim_tam(default_fit, obs_only = TRUE)
  expect_true(is.list(rep_obs))
  expect_true("log_obs" %in% names(rep_obs))
  expect_true(all(is.finite(rep_obs$log_obs[!is.na(rep_obs$log_obs)])))

  # Simulate and also regenerate report using simulated random effects
  rep_full <- sim_tam(default_fit, obs_only = FALSE)
  expect_true(all(c("F", "M", "Z", "ssb", "log_obs") %in% names(rep_full)))
  expect_equal(length(rep_full$ssb), length(YEARS))
})

