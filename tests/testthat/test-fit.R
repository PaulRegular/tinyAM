
## fit_tam ----

testthat::skip_if_not_installed("RTMB")
testthat::skip_if_not(exists("cod_obs"), "cod_obs not available")
testthat::skip_if(is.null(default_fit), "default_fit unavailable")

set.seed(1)

test_that("fit_tam runs on a cod dataset and returns expected structure", {
  fit <- default_fit

  # Structure
  expect_type(fit, "list")
  expect_s3_class(fit, "tam_fit")
  expect_named(
    fit,
    c("call", "dat", "obj", "opt", "rep", "sdrep", "obs_pred", "pop", "is_converged",
      "fixed_par", "random_par", "refit_args", "grad_tol", "diagnostics",
      "sdreport_error", "gradient", "parameter_values", "parameter_map", "bounds"),
    ignore.order = TRUE
  )

  # Optimizer status
  expect_true(is.finite(fit$opt$objective))
  expect_true(fit$opt$objective > 0)
  expect_true(is.list(fit$rep))
  expect_s3_class(fit$sdrep, "sdreport")

  # Reported matrices have expected dims
  expect_equal(dim(fit$rep$F), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$M), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$Z), c(length(YEARS), length(AGES)))
  expect_equal(length(fit$rep$ssb), length(YEARS))
  expect_equal(fit$grad_tol, 0.1)
})


test_that("check_tam detects a fit with a failed Hessian check", {
  # A particular process combination need not fail under every initializer.
  bad_fit <- default_fit
  bad_fit$sdrep$pdHess <- FALSE
  expect_false(check_tam(bad_fit)$is_converged)
})

test_that("fit_tam works when an survey does not provide an index for all ages", {
  obs <- cod_obs
  sub_ages <- 2:10
  obs$index <- obs$index[obs$index$age %in% sub_ages, ]
  fit <- update(
    default_fit,
    data = obs,
    F_settings = list(process = "ar1", mu_form = ~ 1),
    silent = TRUE
  )
  expect_equal(range(fit$obs_pred$index$age), range(sub_ages))
})


test_that("fit_tam objective is unaffected by projections", {
  fit <- update(
    default_fit,
    proj_settings = list(n_proj = 20, n_mean = 20, F_mult = 1),
    start_par = as.list(default_fit$sdrep, "Estimate"),
    silent = TRUE
  )
  expect_equal(fit$opt$objective, default_fit$opt$objective, tolerance = 1e-6)

  # "missing" random effects in projections = predictions
  is_proj <- fit$dat$obs_map$is_proj
  expect_equal(fit$rep$log_obs[is_proj], fit$rep$log_pred[is_proj])
})

test_that("fit_tam does not estimate missing values when fill_missing = FALSE", {
  fit <- update(
    default_fit,
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(sd_form = ~1, q_form = ~ q_block, fill_missing = FALSE),
    silent = TRUE
  )
  expect_false("missing" %in% fit$obj$env$.random)
})

test_that("fit_tam warns and forces fill_missing to TRUE when mising", {
  (fit <- update(
    default_fit,
    catch_settings = list(sd_form = ~1),
    index_settings = list(sd_form = ~1, q_form = ~q_block, fill_missing = TRUE),
    silent = TRUE
  )) |>
    expect_warning(regexp = "catch_settings\\$fill_missing was NULL", fixed  = FALSE)
  expect_true(fit$dat$catch_settings$fill_missing)
  expect_true(fit$dat$index_settings$fill_missing)
})

test_that("update can add arguments absent from the original call", {
  cl <- update(
    default_fit,
    years = 1983:2020,
    silent = TRUE,
    evaluate = FALSE
  )
  expect_true(is.call(cl))
  expect_identical(cl$silent, TRUE)
  expect_identical(cl$years, quote(1983:2020))
})

## fit_retro ----

test_that("fit_retro runs peels and returns stacked outputs", {
  fit <- default_fit
  retros <- fit_retro(fit, folds = 1, progress = FALSE)
  expect_true(is.list(retros))
  expect_true(all(c("obs_pred","pop","fits") %in% names(retros)))
  # At least one retro fit kept (may drop if non-converged)
  if (length(retros$fits) > 0) {
    rf <- retros$fits[[1]]
    expect_true(is.list(rf$rep))
    expect_s3_class(rf$sdrep, "sdreport")
    expect_s3_class(rf, "tam_fit")
  }
})

test_that("fit_retro returns error when no fits converge", {
  fit <- default_fit
  suppressWarnings(fit_retro(fit, folds = 1, progress = FALSE, grad_tol = 1e-16)) |>
    expect_error("All folds failed convergence checks")
})

test_that("fit_retro inherits grad_tol stored on the fit when omitted", {
  fit <- default_fit
  fit$grad_tol <- 1e-16
  suppressWarnings(fit_retro(fit, folds = 1, progress = FALSE)) |>
    expect_error("All folds failed convergence checks")
})

test_that("fit_retro falls back to default grad_tol when fit has none", {
  fit <- default_fit
  fit$grad_tol <- NULL
  fit$sdrep$gradient.fixed[] <- 5e-4

  implicit <- fit_retro(fit, folds = 0, progress = FALSE)
  explicit <- fit_retro(fit, folds = 0, progress = FALSE, grad_tol = 1e-2)

  expect_identical(names(implicit$fits), names(explicit$fits))
  expect_identical(lapply(implicit$fits, `[[`, "is_converged"),
                   lapply(explicit$fits, `[[`, "is_converged"))
  expect_identical(names(implicit$obs_pred), names(explicit$obs_pred))
  expect_identical(names(implicit$pop), names(explicit$pop))
})

test_that("tam_fit summary and print methods provide structured output", {
  fit <- default_fit
  sum_fit <- summary(fit)
  expect_s3_class(sum_fit, "summary_tam_fit")
  expect_equal(sum_fit$convergence, check_tam(fit))
  expect_output(print(fit), "Numerical convergence: passed")
  expect_output(print(sum_fit), "Structural checks:")
  expect_output(print(fit), "Coefficients:")
  expect_output(print(sum_fit), "Terminal year")
})

test_that("fit_hindcasts runs peels with a one year projection", {
  fit <- default_fit
  hindcasts <- suppressWarnings(fit_hindcast(fit, folds = 4, progress = FALSE))
  # At least one hindcast fit kept (may drop if non-converged)
  if (length(hindcasts$fits) > 1) {
    hindcast_year <- as.numeric(names(hindcasts$fits[1]))
    modeled_years <- hindcasts$fits[[1]]$dat$years
    expect_equal(hindcast_year + 1, max(modeled_years))
  }
})
