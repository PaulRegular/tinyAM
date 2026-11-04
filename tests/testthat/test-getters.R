test_that("quantity getters select the maintained annual and age-specific tables", {
  getters <- list(ssb = get_ssb, recruitment = get_recruitment,
    F_bar = get_fbar, M_bar = get_mbar, F = get_f, M = get_m)
  original <- default_fit
  expected <- .population_tables(list(base = default_fit), .95)[[1L]]
  for (name in names(getters)) {
    value <- getters[[name]](default_fit)
    reference <- expected[[name]]
    rownames(reference) <- NULL
    expect_equal(value, reference)
    expect_false("model" %in% names(value))
    expect_true(all(c("year", "est", "is_proj") %in% names(value)))
  }
  expect_true("age" %in% names(get_f(default_fit)))
  expect_false("age" %in% names(get_fbar(default_fit)))
  expect_identical(default_fit, original)
})

test_that("quantity getters retain labels, intervals and projections without diagnostics", {
  testthat::local_mocked_bindings(
    tidy_obs_pred = function(...) stop("Unexpected observation diagnostics"),
    tidy_par = function(...) stop("Unexpected parameter tidying"),
    tidy_recruitment = function(...) stop("Unexpected recruitment diagnostics"),
    .tidy_process_residuals = function(...) stop("Unexpected process diagnostics"))
  original <- N_dev
  value <- get_ssb(Base = default_fit, Forecast = N_dev, interval = .8)
  expect_identical(names(value)[1L], "model")
  expect_identical(unique(value$model), c("Base", "Forecast"))
  expect_true(any(value$is_proj))
  expected <- tidy_pop(N_dev, interval = .8)$ssb
  expect_equal(value$upr[value$model == "Forecast"], expected$upr)
  expect_identical(get_ssb(list(Base = default_fit))$model,
    rep("Base", nrow(default_fit$pop$ssb)))
  expect_identical(get_ssb(model_list = list(`2020` = default_fit))$model,
    rep(2020L, nrow(default_fit$pop$ssb)))
  expect_identical(N_dev, original)
})

test_that("population selectors share explicit reference conversions and SE scales", {
  fit <- default_fit
  fit$pop$ssb$unit <- "tonnes"
  reference <- structure(list(pop = list(ssb = fit$pop$ssb),
    comparison_scales = c(ssb = 1e-3)), class = c("tam_ref", "list"))
  reference$pop$ssb$unit <- "thousand tonnes"
  original <- list(fit = fit, reference = reference)
  value <- get_ssb(model_list = list(Base = fit, Accepted = reference))
  tiny <- value[value$model == "Base", ]
  expect_equal(tiny$est, fit$pop$ssb$est * 1e-3)
  expect_equal(tiny$se, fit$pop$ssb$se)
  expect_true(all(value$unit == "thousand tonnes"))
  fit$pop$ssb$se_scale <- "reported"
  expect_equal(get_ssb(fit, reference)$se[seq_len(nrow(fit$pop$ssb))],
    fit$pop$ssb$se * 1e-3)
  expect_identical(reference, original$reference)
  bad <- reference
  bad$comparison_scales <- NULL
  expect_error(get_ssb(default_fit, bad), "incompatible units")
  bad <- reference
  bad$comparison_scales <- c(ssb = .01)
  expect_error(get_ssb(default_fit, bad), "require")
  missing <- structure(list(pop = list()), class = "tam_ref")
  expect_equal(get_f(default_fit, missing)$est, get_f(default_fit)$est)
  expect_error(get_f(missing), "unavailable")
})

test_that("quantity getters preserve configured means and missing uncertainty", {
  fit <- N_dev
  expect_equal(get_fbar(fit)$est, fit$pop$F_bar$est)
  expect_equal(get_mbar(fit)$est, fit$pop$M_bar$est)
  fit$pop <- NULL
  fit$sdrep <- NULL
  value <- get_ssb(fit)
  expect_true(all(is.na(value$lwr)))
  expect_true(all(is.na(value$upr)))
  expect_equal(value$est, unname(fit$rep$ssb))
  expect_error(get_ssb(fit, interval = 1), "strictly between")
  expect_error(get_ssb(fit, model_list = list(Base = fit)), "either")
  expect_error(get_ssb(model_list = list(fit)), "named")
})
