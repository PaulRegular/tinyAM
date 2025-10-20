
source(system.file("examples/example_fits.R", package = "tinyAM"))
fit <- N_dev

## sim_tam ----

test_that("sim_tam: structure and labeling for redraw_random = FALSE (obs-only)", {
  sims <- sim_tam(
    fit, n = 3,
    par_uncertainty = "fixed",
    redraw_random   = FALSE,
    progress = FALSE, seed = 1
  )

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

test_that("sim_tam: schema is same for redraw_random = TRUE, and values look sensible", {
  sims <- sim_tam(
    fit, n = 2,
    par_uncertainty = "fixed",
    redraw_random   = TRUE,
    progress = FALSE, seed = 1
  )

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
  sims5 <- sim_tam(
    fit, n = 5,
    par_uncertainty = "fixed",
    redraw_random   = FALSE,
    progress = FALSE, seed = 1
  )
  expect_equal(sort(unique(sims5$catch$sim)), 1:5)
  expect_equal(sort(unique(sims5$ssb$sim)),   1:5)
})

test_that("sim_tam works when projecting (is_proj propagated)", {

  ## N deviations
  proj_fit <- update(N_dev, proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1))
  proj_sim <- sim_tam(
    proj_fit, n = 5,
    par_uncertainty = "fixed",
    redraw_random   = FALSE,
    progress = FALSE, seed = 1
  )
  expect_true(any(proj_sim$index$is_proj))

  ## Also check that it works for M deviations
  proj_fit <- update(M_dev, proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1))
  proj_sim <- sim_tam(
    proj_fit, n = 5,
    par_uncertainty = "fixed",
    redraw_random   = FALSE,
    progress = FALSE, seed = 1
  )
  expect_true(any(proj_sim$index$is_proj))
})

test_that("sim_tam: seed makes results deterministic for a given mode", {
  # obs-only, point estimates (no parameter uncertainty)
  sims_a1 <- sim_tam(
    fit, n = 3,
    par_uncertainty = "none",
    redraw_random   = FALSE,
    progress = FALSE, seed = 42
  )
  sims_a2 <- sim_tam(
    fit, n = 3,
    par_uncertainty = "none",
    redraw_random   = FALSE,
    progress = FALSE, seed = 42
  )

  # identical structure and values
  expect_identical(sims_a1$catch$obs, sims_a2$catch$obs)
  expect_identical(sims_a1$index$obs, sims_a2$index$obs)
  expect_identical(sims_a1$ssb$est,   sims_a2$ssb$est)

  # projections with parameter uncertainty + new REs, same seed => same results
  sims_b1 <- sim_tam(
    fit, n = 3,
    par_uncertainty = "fixed",
    redraw_random   = TRUE,
    progress = FALSE, seed = 99
  )
  sims_b2 <- sim_tam(
    fit, n = 3,
    par_uncertainty = "fixed",
    redraw_random   = TRUE,
    progress = FALSE, seed = 99
  )

  expect_identical(sims_b1$catch$obs, sims_b2$catch$obs)
  expect_identical(sims_b1$index$obs, sims_b2$index$obs)
  expect_identical(sims_b1$ssb$est,   sims_b2$ssb$est)
})

test_that("sim_tam: joint + keep behaves (posterior predictive)", {
  sims <- sim_tam(
    fit, n = 2,
    par_uncertainty = "joint",
    redraw_random   = FALSE,   # keep sampled (u,θ), simulate obs only
    progress = FALSE, seed = 7
  )
  expect_true(all(c("catch", "index", "ssb") %in% names(sims)))
  expect_true("sim" %in% names(sims$ssb))
  expect_equal(sort(unique(sims$ssb$sim)), 1:2)
  # sanity: obs should be finite
  expect_true(all(is.finite(sims$catch$obs[!is.na(sims$catch$obs)])))
  expect_true(all(is.finite(sims$index$obs[!is.na(sims$index$obs)])))
})

test_that("sim_tam: par_uncertainty must be one of 'none','fixed','joint'", {
  expect_error(
    sim_tam(fit, n = 1, par_uncertainty = "nope", redraw_random = FALSE,
            progress = FALSE, seed = 1),
    regexp = "arg"
  )
})
