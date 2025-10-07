
fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14, silent = TRUE)

## sim_tam ----

test_that("sim_tam: structure and labeling for obs_only = TRUE", {
  sims <- sim_tam(fit, n = 3, obs_only = TRUE, progress = FALSE, seed = 1)

  # Expected top-level tables exist
  expect_true(all(c("catch", "index") %in% names(sims)))
  # A couple of core pop tables should exist too (from tidy_rep)
  expect_true(all(c("N", "F", "M", "Z", "ssb") %in% names(sims)))

  # Label column present and correct
  expect_true("sim" %in% names(sims$catch))
  expect_true("sim" %in% names(sims$index))
  expect_equal(sort(unique(sims$catch$sim)), 1:3)
  expect_equal(sort(unique(sims$ssb$sim)),   1:3)

  # Simulated obs are finite where not missing
  expect_true(all(is.finite(sims$catch$obs[!is.na(sims$catch$obs)])))
  expect_true(all(is.finite(sims$index$obs[!is.na(sims$index$obs)])))

  # is_proj flag carried through
  expect_true("is_proj" %in% names(sims$catch))
  expect_true("is_proj" %in% names(sims$index))
  expect_type(sims$catch$is_proj, "logical")
})

test_that("sim_tam: schema is same for obs_only = FALSE, and values look sensible", {
  sims <- sim_tam(fit, n = 2, obs_only = FALSE, progress = FALSE, seed = 1)

  # Same key tables exist
  expect_true(all(c("catch", "index", "N", "F", "M", "Z", "ssb") %in% names(sims)))

  # Label column present and correct
  expect_true("sim" %in% names(sims$ssb))
  expect_equal(sort(unique(sims$ssb$sim)), 1:2)

  # Population estimates are finite (not all NA)
  expect_true(any(is.finite(sims$ssb$est)))
  expect_true(any(is.finite(sims$F$est)))
})

test_that("sim_tam: n controls the number of stacked simulations", {
  sims5 <- sim_tam(fit, n = 5, obs_only = TRUE, progress = FALSE, seed = 1)
  expect_equal(sort(unique(sims5$catch$sim)), 1:5)
  expect_equal(sort(unique(sims5$ssb$sim)),   1:5)
})

test_that("sim_tam works when projecting", {
  proj_fit <- update(fit, proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1))
  proj_sim <- sim_tam(proj_fit, n = 5, obs_only = TRUE, progress = FALSE, seed = 1)
  expect_false(all(!proj_sim$index$is_proj))
})
