# A deterministic model at fixed parameter values isolates observation diagnostics
# from optimizer and latent-state approximation errors.
observation_cor_fixture <- function(cor_struct = "ar1") {
  data <- cod_obs
  data$index$obs[data$index$year == 1983 & data$index$age == 3] <- 0
  data$catch$obs[data$catch$year == 1983 & data$catch$age == 3] <- NA
  d <- prepare_tam(data, years = 1983:1987, ages = 2:5,
    N_settings = list(process = "off"),
    rec_settings = list(process = "off", mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~1),
    catch_settings = list(cor_struct = cor_struct),
    index_settings = list(cor_struct = cor_struct, q_form = ~1))
  p <- make_par(d)
  p$log_sd_catch[] <- log(.8)
  p$log_sd_index[] <- log(.6)
  if (!is.null(p$logit_rho_obs)) p$logit_rho_obs[] <- qlogis(.64)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  objective <- obj$fn(obj$par)
  fit <- structure(list(dat = d, obj = obj, rep = obj$report(),
    parameter_values = obj$env$last.par, parameter_map = list(),
    opt = list(par = obj$par, objective = objective),
    is_converged = TRUE), class = "tam_fit")
  fit$obs_pred <- tidy_obs_pred(fit)
  fit
}

observation_cor_covariance <- function(fit) {
  map <- fit$dat$obs_map[fit$dat$is_observed, ]
  group <- paste(map$type, map$year)
  sd <- fit$rep$sd_obs[fit$dat$is_observed]
  tcrossprod(sd) * .64^abs(outer(map$age, map$age, "-")) * outer(group, group, "==")
}

test_that("OSA residuals use the marginal covariance of the included observations", {
  fit <- observation_cor_fixture()
  keep <- fit$dat$is_observed
  map <- fit$dat$obs_map[keep, ]
  row <- function(type, age) which(map$type == type & map$year == 1983 & map$age == age)
  selected <- c(row("index", 4), row("catch", 4), row("index", 2))
  conditioned <- row("index", 5)
  residual <- fit$dat$observed - fit$rep$log_pred[keep]
  covariance <- observation_cor_covariance(fit)
  expected <- vapply(seq_along(selected), function(i) {
    current <- selected[i]
    earlier <- c(conditioned, head(selected, i - 1L))
    coef <- solve(covariance[earlier, earlier, drop = FALSE], covariance[earlier, current])
    (residual[current] - sum(coef * residual[earlier])) /
      sqrt(covariance[current, current] - sum(coef * covariance[earlier, current]))
  }, numeric(1))

  for (method in c("oneStepGaussianOffMode", "fullGaussian")) {
    tables <- tidy_obs_pred(fit, add_osa_res = TRUE, method = method,
      subset = selected, conditional = conditioned, trace = FALSE)
    actual <- c(tables$catch$osa_res, tables$index$osa_res)[keep]
    expect_equal(actual[selected], expected, tolerance = 1e-7)
    expect_true(all(is.na(actual[-selected])))
  }
  expect_error(residuals(fit, type = "osa", method = "cdf"),
    "cdf.*does not support correlated")
  expect_equal(c(fit$obs_pred$catch$std_res, fit$obs_pred$index$std_res)[keep],
    residual / fit$rep$sd_obs[keep])
})

test_that("conditional and process simulations both draw correlated observation errors", {
  fit <- observation_cor_fixture()
  keep <- fit$dat$is_observed
  simulate <- function(redraw_random) {
    withr::with_seed(4, tinyAM:::.sim_obs(fit, tinyAM:::.draw_none, redraw_random))
  }
  conditional <- simulate(FALSE)
  process <- simulate(TRUE)
  observations <- c(conditional$catch$obs, conditional$index$obs)
  expect_identical(is.finite(observations), keep)
  expect_equal(observations, c(process$catch$obs, process$index$obs))
  independent <- withr::with_seed(4,
    stats::rnorm(sum(keep), fit$rep$log_pred[keep], fit$rep$sd_obs[keep]))
  expect_gt(max(abs(log(observations[keep]) - independent)), .01)
})

test_that("correlation parameters count toward AIC and conditional fixed EDF", {
  correlated <- observation_cor_fixture()
  independent <- observation_cor_fixture("independent")
  expect_equal(AIC(correlated, independent)$df[1] - AIC(correlated, independent)$df[2], 2)
  cor_ic <- cAIC(correlated, details = TRUE)
  ind_ic <- cAIC(independent, details = TRUE)
  expect_equal(cor_ic$edf_fixed - ind_ic$edf_fixed, 2)
  expect_equal(cor_ic$EDF, cor_ic$edf_fixed)
  keep <- correlated$dat$is_observed
  covariance <- observation_cor_covariance(correlated)
  residual <- correlated$dat$observed - correlated$rep$log_pred[keep]
  expected <- -.5 * (length(residual) * log(2 * pi) +
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
    drop(crossprod(residual, solve(covariance, residual))))
  expect_equal(cor_ic$logLik, expected)
})
