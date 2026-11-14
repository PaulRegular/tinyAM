pld_diagnostic_rows <- function() {
  data.frame(year = 2000:2005, age = 2L, survey = "RV", unit = "fish",
    obs = c(0, 2, NA, 0, 0, 0), pred = c(1, 1.5, 2, 0, 2, 1000),
    family = "pld_lognormal", is_observed = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
    is_structural = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    is_proj = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    log_p_zero = c(log(.5), log(.2), log(.3), 0, log(.1), -1000),
    log_p_positive = c(log(.5), log(.8), log(.7), -Inf, log(.9), 0),
    sd = .5, log_res = c(NA, 0, NA, NA, NA, NA), std_res = c(NA, 0, NA, NA, NA, NA))
}

test_that("conditional quantile residuals include sampled zeros and use stable tails", {
  d <- pld_diagnostic_rows()
  old_seed <- .Random.seed
  old_kind <- RNGkind()
  r <- .add_quantile_residuals(list(index = d), seed = 1L)
  expect_identical(.Random.seed, old_seed)
  expect_identical(RNGkind(), old_kind)
  expect_identical(r, .add_quantile_residuals(list(index = d), seed = 1L))
  rr <- .add_quantile_residuals(list(index = d), seed = 2L)
  expect_false(identical(r$index$quantile_res[1L], rr$index$quantile_res[1L]))
  expect_lt(r$index$quantile_res[1L], 0)
  expect_equal(r$index$quantile_res[2L], qnorm(.2 + .8 * .5))
  expect_identical(which(is.na(r$index$quantile_res)), 3:5)
  expect_true(is.finite(r$index$quantile_res[6L]))
  expect_lt(r$index$quantile_res[6L], -40)
  d$obs[1L] <- 1
  d$log_p_zero[1L] <- -1200
  d$log_p_positive[1L] <- 0
  d$std_res[1L] <- -42
  expect_equal(.add_quantile_residuals(list(index = d))$index$quantile_res[1L],
    -42, tolerance = 1e-10)
})

test_that("retained quantile residuals preserve keys, seeds, and the supplied fit", {
  d <- pld_diagnostic_rows()
  fit <- structure(list(dat = list(pld = list(series = list(list()))),
    obs_pred = list(catch = d[FALSE, ], index = d)), class = "tam_fit")
  fit2 <- add_residuals(fit, type = "quantile", residual_seed = 97L)
  expect_null(fit$obs_pred$index$quantile_res)
  res <- residuals(fit2, type = "quantile")
  expect_equal(res$residual, fit2$obs_pred$index$quantile_res)
  expect_identical(attr(res, "residual_seed"), 97L)
  expect_equal(res$year, d$year)
  expect_equal(residuals(fit2, type = "response")$residual[1L], -1)
  expect_equal(residuals(fit2, type = "response")$residual[4L], 0)
  expect_equal(residuals(fit2, type = "log")$residual[2L], 0)
  expect_error(residuals(fit2, type = "osa"), "all-lognormal")
  expect_error(add_osa_residuals(fit2), "all-lognormal")
  expect_error(add_residuals(fit2, residual_seed = .5), "nonnegative integer")
})

test_that("delta displays include zero residuals and maintain calibration support", {
  d <- pld_diagnostic_rows()
  d$model <- factor("delta")
  d$p_zero <- exp(d$log_p_zero)
  d$pred_positive_mean <- d$pred / exp(d$log_p_positive)
  d$pred_positive_mean[d$is_structural] <- NA_real_
  d <- .add_quantile_residuals(list(index = d))$index
  config <- plot_resid(d)$jsHooks$render[[1L]]$data
  expect_identical(config$view, 2L)
  expect_identical(config$selectors, c(0L, 2L))
  expect_identical(config$labels[config$selectors + 1L], c("Standardized", "Quantile"))
  expect_identical(config$views[[3L]]$cells[[1L]]$status$n, 3L)
  frequencies <- .zero_frequency_data(d)
  expect_equal(frequencies$n, 3L)
  expect_equal(frequencies$obs, 2 / 3)
  expect_equal(frequencies$pred, (.5 + .2) / 3)
  traces <- plotly::plotly_build(plot_obs_pred(d[!is.na(d$obs), ]))$x$data
  expect_true(any(vapply(traces, function(trace) any(trace$y == 0, na.rm = TRUE), logical(1))))
  expect_true(any(grepl("mean including zeros", unlist(lapply(traces, `[[`, "text")))))
  expect_match(.obs_family_hover(d)[1L], "Probability of zero: 0.5", fixed = TRUE)
})
