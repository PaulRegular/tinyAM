
## sim_tam ----

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

