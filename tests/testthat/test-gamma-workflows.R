gamma_workflow_inputs <- function() {
  withr::local_seed(516)
  years <- 2000:2014
  ages <- 2:6
  grid <- expand.grid(year = years, age = ages)
  F <- c(0, rep(.2, length(ages) - 1L))
  Z <- F + .3
  N <- matrix(1000, length(years), length(ages))
  N[1, -1L] <- 1000 * exp(-cumsum(utils::head(Z, -1L)))
  for (i in 2:length(years)) {
    N[i, -1L] <- N[i - 1L, -length(ages)] * exp(-utils::head(Z, -1L))
    N[i, length(ages)] <- N[i, length(ages)] + N[i - 1L, length(ages)] * exp(-utils::tail(Z, 1L))
  }
  catch_mean <- as.numeric(sweep(N, 2, F * -expm1(-Z) / Z, "*"))
  index_mean <- as.numeric(.4 * N)
  caught <- stats::rgamma(nrow(grid), shape = .25^-2, scale = catch_mean * .25^2)
  detected <- -expm1(-log(2) * index_mean / 240)
  positive <- stats::runif(nrow(grid)) < detected
  gg_scale <- .36
  Q <- -.20
  k <- Q^-2
  a <- gg_scale / Q
  correction <- lgamma(k + a) - lgamma(k) - a * log(k)
  index <- numeric(nrow(grid))
  index[positive] <- index_mean[positive] / detected[positive] *
    exp(-correction + a * (log(stats::rgamma(sum(positive), k)) - log(k)))
  caught[grid$year == min(years) & grid$age == 4L] <- NA_real_
  index[grid$year == min(years) & grid$age == 5L] <- NA_real_
  make_tam_data(
    catch = transform(grid, obs = caught, unit = "thousand fish"),
    index = transform(grid, obs = index, survey = "RV", samp_time = 0, unit = "thousand fish"),
    weight = transform(grid, obs = as.numeric(age), unit = "kg"),
    maturity = transform(grid, obs = as.numeric(age >= 4))
  )
}

gamma_workflow_settings <- function(data = gamma_workflow_inputs()) {
  list(data = data,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "off", mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~1, age_breaks = 3:7),
    M_settings = list(process = "off", mu_supplied = ~I(.3)),
    catch_settings = list(family = stats::Gamma(link = "log")),
    index_settings = list(family = list(RV = delta_gengamma(type = "poisson-link", y50 = 240, Q = -.20)),
      q_form = ~1, sd_form = ~0, sd_supplied = ~I(.36)))
}

test_that("Gamma and Delta generalized gamma work together through fit, simulation and forecasts", {
  settings <- gamma_workflow_settings()
  fit <- do.call(fit_tam, c(settings, list(silent = TRUE)))
  expect_true(fit$is_converged)
  expect_true(is.finite(fit$opt$objective))
  catch <- fit$obs_pred$catch
  index <- fit$obs_pred$index
  eligible <- catch$is_observed
  cv <- catch$sd[eligible]
  mean <- catch$pred[eligible]
  expect_identical(unique(catch$sd_scale), "CV")
  expect_equal(catch$pred_positive_median[eligible],
    stats::qgamma(.5, shape = cv^-2, scale = mean * cv^2))
  expect_equal(catch$quantile_res[eligible],
    stats::qnorm(stats::pgamma(catch$obs[eligible], shape = cv^-2, scale = mean * cv^2)))
  expect_equal(catch$std_res[eligible], (log(catch$obs[eligible]) -
    (log(mean) + digamma(cv^-2) - log(cv^-2))) / sqrt(trigamma(cv^-2)))

  expect_identical(unique(index$sd_scale), "Prentice scale")
  expect_equal(index$pred_mean, index$pred)
  expect_equal(index$p_positive * index$pred_positive_mean, index$pred_mean)
  shape <- subset(fit$fixed_par, par == "gengamma_Q")
  expect_identical(shape$coef, "RV")
  expect_true(shape$is_fixed)
  expect_equal(shape$est, -.20)
  expect_false(any(c("eta_gengamma_Q_catch", "eta_gengamma_Q_index") %in% fit$fixed_par$par))
  expect_false("gengamma_Q_report" %in% names(tidy_sdrep(fit)))
  expect_error(residuals(fit, type = "osa"), "all-lognormal")
  expect_true(all(is.finite(residuals(fit, type = "quantile")$residual[
    c(catch$is_observed, index$is_observed)])))
  expect_true(all(is.na(catch$quantile_res[catch$is_structural])))
  expect_true(is.finite(stats::AIC(fit)))
  expect_true(is.finite(cAIC(fit)))
  expect_no_warning(stats::AIC(fit, fit))

  per_fish <- settings
  for (type in c("catch", "index")) {
    per_fish$data[[type]]$obs <- 1000 * per_fish$data[[type]]$obs
    per_fish$data[[type]]$unit <- "fish"
  }
  per_fish$index_settings$family <- delta_gengamma(type = "poisson-link", y50 = 240000, Q = -.20)
  scaled <- do.call(prepare_tam, per_fish)
  par <- fit$obj$env$parList(par = fit$parameter_values)
  par$log_r0 <- par$log_r0 + log(1000)
  par$rec_beta <- par$rec_beta + log(1000)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, scaled), par, silent = TRUE)
  expect_equal(obj$fn(), fit$opt$objective, tolerance = 1e-10)

  for (redraw in c(FALSE, TRUE)) {
    sim <- sim_tam(fit, n = 2L, seed = 81, par_uncertainty = "none",
      redraw_random = redraw, progress = FALSE)
    expect_true(all(sim$catch$obs[sim$catch$age == 2L] == 0))
    expect_identical(is.na(sim$catch$obs), rep(is.na(catch$obs), 2L))
    expect_identical(is.na(sim$index$obs), rep(is.na(index$obs), 2L))
    expect_true(any(sim$index$obs == 0, na.rm = TRUE))
    expect_true(any(sim$index$obs > 0, na.rm = TRUE))
  }

  hind <- fit_hindcast(fit, folds = 1, proj_settings = list(F_mult = 1), progress = FALSE)
  index_forecast <- subset(hind$obs_pred$index, is_proj)
  expect_gt(nrow(index_forecast), 0)
  expect_equal(.forecast_log_pred(index_forecast), index_forecast$log_pred_positive)
  expect_true(all(is.finite(index_forecast$p_zero)))
  expect_equal(.occurrence_forecasts(transform(index_forecast, obs = 0))$brier,
    (1 - index_forecast$p_zero)^2)
})

test_that("Delta Gamma and estimated generalized gamma shape are available in fitted models", {
  settings <- gamma_workflow_settings()
  settings$catch_settings <- list(family = delta_gamma(type = "poisson-link", y50 = 80))
  settings$index_settings$family <- gengamma()
  fit <- do.call(fit_tam, c(settings, list(silent = TRUE)))
  expect_true(fit$is_converged)
  expect_true(all(is.finite(fit$obj$gr(fit$opt$par))))
  shape <- subset(fit$fixed_par, par == "gengamma_Q")
  expect_identical(shape$coef, "RV")
  expect_false(shape$is_fixed)
  expect_true(all(is.finite(unlist(shape[c("est", "lwr", "upr", "se")]))))
  expect_gt(shape$se, 0)
  expect_true(all(fit$obs_pred$index$p_zero == 0))
  structural <- subset(fit$obs_pred$catch, is_structural)
  expect_true(all(structural$p_zero == 1 & structural$pred_mean == 0))
  expect_true(all(is.na(structural$quantile_res)))
  catch <- subset(fit$obs_pred$catch, !is_structural)
  expect_equal(catch$p_positive * catch$pred_positive_mean, catch$pred)
  expect_true(all(catch$p_zero > 0 & catch$p_zero < 1))
})

test_that("family assignments protect observation comparisons and dispersion semantics", {
  settings <- gamma_workflow_settings()
  data <- settings$data
  auxiliary <- transform(data$index, survey = "Aux")
  data$index <- rbind(data$index, auxiliary)
  settings$data <- data
  settings$index_settings <- list(family = list(RV = gengamma()), q_form = ~survey)
  expect_error(do.call(prepare_tam, settings), "dispersion|different")
  settings$index_settings$sd_form <- ~ 0 + survey
  dat <- do.call(prepare_tam, settings)
  expect_identical(names(make_par(dat)$eta_gengamma_Q_index), "RV")
  settings$index_settings$cor_struct <- c(RV = "ar1")
  expect_error(do.call(prepare_tam, settings), "requires independent")
  settings$index_settings$cor_struct <- c(Aux = "ar1")
  mixed <- do.call(prepare_tam, settings)
  expect_true(.obs_has_correlation(mixed))
  tape <- RTMB::MakeADFun(function(p) nll_fun(p, mixed), make_par(mixed), silent = TRUE)
  expect_true(all(is.finite(c(tape$fn(), tape$gr()))))

  settings <- gamma_workflow_settings()
  settings$index_settings$family <- gengamma()
  settings$index_settings$sd_form <- ~factor(age)
  settings$index_settings$sd_supplied <- NULL
  expect_error(do.call(prepare_tam, settings), "constant")
  settings$index_settings$family <- gengamma(Q = -3)
  settings$index_settings$sd_form <- ~1
  expect_error(do.call(prepare_tam, settings), "supplied")
  settings$index_settings$sd_form <- ~0
  settings$index_settings$sd_supplied <- ~I(.5)
  expect_error(do.call(prepare_tam, settings), "finite.*mean|mean.*finite")

  settings <- gamma_workflow_settings()
  settings$catch_settings$cor_struct <- "ar1"
  expect_error(do.call(prepare_tam, settings), "requires independent")
})

test_that("new Delta quantile residuals stay finite at low encounter intensity", {
  for (kind in c("delta_gamma", "delta_gengamma")) {
    s <- .4
    # At this shape, the generalized gamma positive component is Gamma.
    d <- data.frame(family = kind, obs = c(0, .25, .25),
      is_observed = TRUE, is_structural = FALSE,
      log_p_zero = c(-exp(-800), -exp(-800), log(.4)),
      log_p_positive = c(-800, -800, log(.6)),
      sd = s, gengamma_Q = s, log_obs_positive_mean = 0)
    r <- .add_quantile_residuals(list(index = d), seed = 11L)$index$quantile_res
    expect_true(all(is.finite(r)))
    expect_equal(r[1], stats::qnorm(.residual_uniform(1L, 11L)))
    expect_gt(r[2], 35)
    expect_equal(r[3], stats::qnorm(.4 + .6 * stats::pgamma(.25, s^-2, scale = s^2)),
      tolerance = 1e-12)
  }
})

