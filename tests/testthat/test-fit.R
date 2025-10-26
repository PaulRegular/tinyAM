
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
      "fixed_par", "random_par", "grad_tol"),
    ignore.order = TRUE
  )

  # Optimizer status
  expect_true(is.finite(fit$opt$objective))
  expect_true(fit$opt$objective > 0)
  expect_equal(round(fit$opt$objective, 4), 989.6083)
  expect_true(is.list(fit$rep))
  expect_s3_class(fit$sdrep, "sdreport")

  # Reported matrices have expected dims
  expect_equal(dim(fit$rep$F), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$M), c(length(YEARS), length(AGES)))
  expect_equal(dim(fit$rep$Z), c(length(YEARS), length(AGES)))
  expect_equal(length(fit$rep$ssb), length(YEARS))
  expect_equal(fit$grad_tol, 1e-2)
})


test_that("fit_tam emits warning if the model does not converge", {
  bad_fit <- suppressWarnings({
    update(
      default_fit,
      N_settings = list(process = "iid", init_N0 = TRUE),
      F_settings = list(process = "iid"),
      M_settings = list(process = "iid", mu_form = NULL, assumption = ~I(0.3)),
      silent = TRUE
    )
  })
  expect_warning(check_convergence(bad_fit, quiet = TRUE),
                 regexp = "Model may not have converged", fixed = FALSE) # one of many warnings
})

test_that("fit_tam works when an survey does not provide an index for all ages", {
  obs <- cod_obs
  sub_ages <- 2:10
  obs$index <- obs$index[obs$index$age %in% sub_ages, ]
  fit <- update(default_fit, obs = obs, silent = TRUE)
  expect_equal(range(fit$obs_pred$index$age), range(sub_ages))
})


test_that("fit_tam objective is unaffected by projections", {
  fit <- update(
    default_fit,
    proj_settings = list(n_proj = 20, n_mean = 20, F_mult = 1),
    silent = TRUE
  )
  expect_equal(round(fit$opt$objective, 4), 989.6083)

  # "missing" random effects in projections = predictions
  is_proj <- fit$dat$obs_map$is_proj
  expect_equal(fit$rep$log_obs[is_proj], fit$rep$log_pred[is_proj])
})

test_that("fit_tam does not estimate missing values when fill_missing = FALSE", {
  fit <- update(
    default_fit,
    obs_settings = list(q_form = ~ q_block, sd_catch_form = ~1,
                        sd_index_form = ~1, fill_missing = FALSE),
    silent = TRUE
  )
  expect_false("missing" %in% fit$obj$env$.random)
})

test_that("fit_tam warns and forces fill_missing to TRUE when mising", {
  (fit <- update(default_fit, obs_settings = list(q_form = ~q_block, sd_catch_form = ~1, sd_index_form = ~1), silent = TRUE)) |>
    expect_warning(regexp = "fill_missing was NULL", fixed  = FALSE)
  expect_true(fit$dat$obs_settings$fill_missing)
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

test_that("fit_retro returns structured empty results when no fits converge", {
  fit <- default_fit
  retros <- fit_retro(fit, folds = 1, progress = FALSE, grad_tol = 0)

  expect_true(is.list(retros))
  expect_equal(retros$fits, list())
  expect_equal(retros$obs_pred, list())
  expect_equal(retros$pop, list())
  expect_s3_class(retros$mohns_rho, "data.frame")
  expect_equal(nrow(retros$mohns_rho), 0)
})

test_that("fit_retro inherits grad_tol stored on the fit when omitted", {
  fit <- default_fit
  fit$grad_tol <- 0.0123
  captured <- list()

  original_check <- tinyAM:::check_convergence

  with_mocked_bindings(
    check_convergence = function(fit, grad_tol, ...) {
      captured <<- c(captured, list(grad_tol))
      original_check(fit, grad_tol = grad_tol, ...)
    }, {
      fit_retro(fit, folds = 1, progress = FALSE)
    },
    .package = "tinyAM"
  )

  expect_length(captured, 1)
  expect_equal(captured[[1]], 0.0123)
})

test_that("fit_retro falls back to default grad_tol when fit has none", {
  fit <- default_fit
  fit$grad_tol <- NULL
  captured <- list()

  original_check <- tinyAM:::check_convergence

  with_mocked_bindings(
    check_convergence = function(fit, grad_tol, ...) {
      captured <<- c(captured, list(grad_tol))
      original_check(fit, grad_tol = grad_tol, ...)
    }, {
      fit_retro(fit, folds = 1, progress = FALSE)
    },
    .package = "tinyAM"
  )

  expect_length(captured, 1)
  expect_equal(captured[[1]], 1e-3)
})

test_that("tam_fit summary and print methods provide structured output", {
  fit <- default_fit
  sum_fit <- summary(fit)
  expect_s3_class(sum_fit, "summary_tam_fit")
  expect_output(print(fit), "Coefficients:")
  expect_output(print(sum_fit), "Terminal year")
})

test_that("fit_hindcasts runs peels with a one year projection", {
  fit <- default_fit
  hindcasts <- fit_hindcast(fit, folds = 1, progress = FALSE, grad_tol = 0.01)
  # At least one hindcast fit kept (may drop if non-converged)
  if (length(hindcasts$fits) > 0) {
    hindcast_year <- as.numeric(names(hindcasts$fits[1]))
    modeled_years <- hindcasts$fits[[1]]$dat$years
    expect_equal(hindcast_year + 1, max(modeled_years))
  }
})

test_that("hindcast empty fits surface placeholder RMSE", {
  fit <- default_fit
  hindcasts <- fit_retro(fit, folds = 1, hindcast = TRUE, progress = FALSE, grad_tol = 0)

  expect_true("hindcast_rmse" %in% names(hindcasts))
  expect_true(is.na(hindcasts$hindcast_rmse))
  expect_equal(hindcasts$fits, list())
  expect_equal(hindcasts$obs_pred, list())
  expect_equal(hindcasts$pop, list())
  expect_s3_class(hindcasts$mohns_rho, "data.frame")
  expect_equal(nrow(hindcasts$mohns_rho), 0)
  expect_s3_class(hindcasts$fixed_par, "data.frame")
  expect_equal(nrow(hindcasts$fixed_par), 0)
  expect_equal(hindcasts$random_par, list())
})


## check_convergence ----

make_conv_fit <- function(max_grad, pd_hess = TRUE) {
  fit <- default_fit
  fit$sdrep$gradient.fixed <- rep(max_grad, length(fit$sdrep$gradient.fixed))
  fit$sdrep$pdHess <- pd_hess
  fit
}

test_that("check_convergence returns TRUE and messages when all checks pass", {
  fit_ok <- make_conv_fit(max_grad = 5e-5, pd_hess = TRUE)
  expect_message(
    val <- check_convergence(fit_ok, grad_tol = 1e-3, quiet = FALSE),
    "Model converged"
  )
  expect_true(val)
})

test_that("check_convergence is silent on success when quiet = TRUE", {
  fit_ok <- make_conv_fit(max_grad = 5e-5, pd_hess = TRUE)
  expect_silent(check_convergence(fit_ok, grad_tol = 1e-3, quiet = TRUE))
})

test_that("check_convergence warns and returns FALSE if gradient too large", {
  fit_bad_grad <- make_conv_fit(max_grad = 1e-1, pd_hess = TRUE)
  expect_warning(
    val <- check_convergence(fit_bad_grad, grad_tol = 1e-3, quiet = TRUE),
    "Model may not have converged"
  )
  expect_false(val)
})

test_that("check_convergence warns and returns FALSE if Hessian not PD", {
  fit_bad_hess <- make_conv_fit(max_grad = 5e-5, pd_hess = FALSE)
  expect_warning(
    val <- check_convergence(fit_bad_hess, quiet = TRUE),
    "^Model may not have converged"
  )
  expect_false(val)
})

test_that("check_convergence accepts sdreport objects and sdrep lists", {
  fit_ok <- make_conv_fit(max_grad = 5e-5, pd_hess = TRUE)

  expect_true(check_convergence(fit_ok$sdrep, grad_tol = 1e-3, quiet = TRUE))

  sdrep_list <- list(
    gradient.fixed = fit_ok$sdrep$gradient.fixed,
    pdHess = fit_ok$sdrep$pdHess
  )

  expect_true(check_convergence(list(sdrep = sdrep_list), grad_tol = 1e-3, quiet = TRUE))
})

test_that("check_convergence errors when sdrep details are missing", {
  expect_error(
    check_convergence(list()),
    "must be either",
    class = "rlang_error"
  )

  expect_error(
    check_convergence(list(sdrep = list(pdHess = TRUE))),
    "must provide a gradient",
    class = "rlang_error"
  )

  expect_error(
    check_convergence(list(sdrep = list(gradient.fixed = 0))),
    "must provide a Hessian flag",
    class = "rlang_error"
  )
})

