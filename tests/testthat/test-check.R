make_check_fit <- function() {
  structure(list(
    opt = list(par = c(beta = 1), objective = 12, convergence = 0, message = "OK"),
    sdrep = list(par.fixed = c(beta = 1), gradient.fixed = 1e-5,
                 cov.fixed = matrix(.1), pdHess = TRUE),
    dat = list(), rep = list(log_pred = c(0, 1), sd_obs = c(.1, .2)),
    obs_pred = list(), grad_tol = 1e-2
  ), class = "tam_fit")
}

test_that("numerical diagnostics include optimizer status and use one tolerance", {
  fit <- make_check_fit()
  original <- fit
  checks <- check_tam(fit)
  expect_s3_class(checks, "tam_check")
  expect_true(checks$is_converged)
  expect_true(all(checks$numerical$status == "pass"))
  expect_identical(fit, original)
  expect_message(print(checks), "Numerical convergence: passed")
  fit$opt$convergence <- 1
  expect_false(check_tam(fit)$is_converged)
  expect_identical(check_tam(fit)$numerical$status[1], "fail")
  fit$opt$convergence <- 0
  fit$sdrep$gradient.fixed <- .005
  expect_true(check_tam(fit)$is_converged)
  expect_false(check_tam(fit, grad_tol = .001)$is_converged)
  fit$grad_tol <- NULL
  expect_equal(check_tam(fit)$grad_tol, .01)
})

test_that("unavailable and non-finite numerical information cannot pass", {
  fit <- make_check_fit()
  for (field in c("objective", "par")) {
    bad <- fit
    bad$opt[[field]][] <- Inf
    expect_false(check_tam(bad)$is_converged)
  }
  for (gradient in c(NA_real_, NaN, Inf, -Inf, .1)) {
    bad <- fit
    bad$sdrep$gradient.fixed <- gradient
    expect_false(check_tam(bad)$is_converged)
  }
  bad <- fit
  bad$sdrep$pdHess <- FALSE
  expect_false(check_tam(bad)$is_converged)
  bad <- fit
  bad$sdrep$cov.fixed[] <- NA_real_
  expect_false(check_tam(bad)$is_converged)
  bad <- fit
  bad$sdrep <- NULL
  bad$gradient <- 1e-5
  bad$sdreport_error <- "Uncertainty calculation failed."
  checks <- check_tam(bad)
  expect_false(checks$is_converged)
  expect_true(any(checks$numerical$status == "not_assessed"))
  expect_true(any(checks$numerical$detail == bad$sdreport_error))
  bad <- fit
  bad$rep$sd_obs[1] <- 0
  expect_false(check_tam(bad)$is_converged)
  bad <- fit
  bad$opt <- NULL
  expect_false(check_tam(bad)$is_converged)
})

test_that("bound-adjusted gradients retain raw and non-finite derivatives", {
  fit <- make_check_fit()
  fit$opt$par <- c(dq = 0)
  fit$sdrep$par.fixed <- fit$opt$par
  fit$sdrep$gradient.fixed <- 2
  checks <- check_tam(fit)
  expect_true(checks$is_converged)
  expect_equal(checks$max_gradient, 0)
  expect_equal(checks$raw_max_gradient, 2)
  expect_equal(fit$sdrep$gradient.fixed, 2)
  expect_true(checks$active_bounds)
  expect_true("active_bounds" %in% checks$advisories$issue)
  for (gradient in c(-2, Inf, NaN)) {
    fit$sdrep$gradient.fixed <- gradient
    expect_false(check_tam(fit)$is_converged)
  }
  fit$opt$par <- c(beta = 1)
  fit$bounds <- list(lower = -Inf, upper = 1)
  fit$sdrep$gradient.fixed <- -2
  expect_true(check_tam(fit)$is_converged)
  fit$sdrep$gradient.fixed <- 2
  expect_false(check_tam(fit)$is_converged)
})

test_that("structural checks distinguish redundancies from convergence", {
  fit <- make_check_fit()
  fit$dat <- list(obs = list(index = data.frame(obs = 1:3, year = 2000:2002)),
    index_settings = list(q_link = "log"), q_modmat = matrix(1, 3, 2))
  checks <- check_tam(fit)
  expect_true(checks$is_converged)
  expect_identical(checks$structural_status, "issues")
  expect_match(checks$structural$detail, "Rank 1 of 2")
  fit$dat$q_modmat <- cbind(1, c(1, 2, 3))
  expect_identical(check_tam(fit)$structural_status, "no_known_issues")
  fit$dat$q_modmat <- cbind(1, c(1, 2, 3) * 1e9)
  expect_identical(check_tam(fit)$structural_status, "no_known_issues")
  fit$dat$q_modmat <- matrix(1, 3, 2)
  fit$parameter_map <- list(log_q = factor(c(1, NA)))
  expect_identical(check_tam(fit)$structural_status, "no_known_issues")
  fit$dat$obs$index$obs[] <- NA_real_
  expect_identical(check_tam(fit)$structural_status, "not_assessed")
})

test_that("residual summaries exclude projections, zero and missing observations", {
  fit <- make_check_fit()
  fit$obs_pred$index <- data.frame(survey = "A", age = 2, year = 2000:2014,
    obs = c(rep(1, 12), 0, NA, 1), std_res = c(rep(1, 12), 9, 9, 9),
    is_proj = c(rep(FALSE, 14), TRUE))
  checks <- check_tam(fit)
  expect_equal(checks$residuals$n, 12)
  expect_equal(checks$residuals$mean, 1)
  expect_equal(checks$residuals$max_abs, 1)
  expect_true(checks$is_converged)
  expect_true("residual_location" %in% checks$advisories$issue)
  expect_true("residual_spread" %in% checks$advisories$issue)
  fit$obs_pred$index$year <- seq(2000, by = 2, length.out = 15)
  fit$obs_pred$index$std_res <- seq_len(15)
  expect_true(is.na(check_tam(fit)$residuals$lag1))
})

test_that("detailed curvature checks describe weak directions without refitting", {
  fit <- make_check_fit()
  fit$opt$par <- c(a = 1, b = 1)
  fit$sdrep$par.fixed <- fit$opt$par
  fit$sdrep$gradient.fixed <- c(0, 0)
  fit$sdrep$cov.fixed <- diag(c(1, 1000))
  expect_null(check_tam(fit)$curvature)
  checks <- check_tam(fit, detailed = TRUE)
  expect_true(checks$is_converged)
  expect_equal(checks$curvature$hessian_condition, 1000)
  expect_identical(checks$curvature$weak_direction$parameter[1], "b")
  fit$sdrep$cov.fixed <- matrix(c(1, .99, .99, 1), 2)
  checks <- check_tam(fit, detailed = TRUE)
  expect_true("parameter_correlation" %in% checks$advisories$issue)
  expect_equal(checks$curvature$correlations$correlation, .99)
})

test_that("check_tam validates its public arguments", {
  expect_error(check_tam(list()), "tam_fit")
  expect_error(check_tam(structure(list(), class = "tam_ref")), "tam_fit")
  for (tol in list(0, -1, NA_real_, c(.1, .2), "small")) {
    expect_error(check_tam(make_check_fit(), grad_tol = tol), "positive finite")
  }
  expect_error(check_tam(make_check_fit(), detailed = NA), "TRUE or FALSE")
  expect_false("check_convergence" %in% getNamespaceExports("tinyAM"))
})

test_that("uncertainty failures preserve estimates without inventing SEs", {
  fit <- default_fit
  estimates <- as.list(fit$sdrep, "Estimate")
  fit$sdrep <- NULL
  fit$sdreport_error <- "test failure"
  expect_false(check_tam(fit)$is_converged)
  tab <- tidy_par(fit)
  expect_true(all(is.na(tab$fixed$se)))
  expect_true(all(is.na(tab$fixed$lwr)))
  expect_equal(tinyAM:::.tam_parameter_summary(fit, "Estimate"), estimates)
  pop <- tidy_pop(fit)
  expect_equal(pop$ssb$est, unname(fit$rep$ssb))
  expect_true(all(is.na(pop$ssb$se)))
  expect_error(tinyAM:::.draw_fixed(fit), "successful uncertainty")
  expect_equal(tinyAM:::.draw_none(fit), estimates)
})

test_that("fit_tam retains a usable fit when sdreport fails", {
  local_mocked_bindings(sdreport = function(...) stop("Deliberate uncertainty failure"),
                        .package = "RTMB")
  expect_warning(fit <- update(default_fit, silent = TRUE,
    start_par = as.list(default_fit$sdrep, "Estimate")), "Numerical convergence")
  expect_s3_class(fit, "tam_fit")
  expect_null(fit[["sdrep"]])
  expect_match(fit$sdreport_error, "Deliberate uncertainty failure")
  expect_false(fit$is_converged)
  expect_true(all(is.finite(fit$pop$ssb$est)))
  expect_true(all(is.na(fit$pop$ssb$se)))
  expect_true(all(is.na(fit$fixed_par$se)))
})
