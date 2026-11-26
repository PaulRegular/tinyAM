tweedie_reference_case <- function() {
  withr::local_seed(621)
  years <- 1983:2022
  ages <- 2:6
  grid <- expand.grid(year = years, age = ages)
  F <- c(0, rep(.26, length(ages) - 1L))
  Z <- F + .3
  N <- matrix(100, length(years), length(ages))
  N[1L, -1L] <- 100 * exp(-cumsum(utils::head(Z, -1L)))
  youngest <- seq_len(length(ages) - 1L)
  oldest <- length(ages)
  for (i in 2:length(years)) {
    N[i, -1L] <- N[i - 1L, youngest] * exp(-Z[youngest])
    N[i, oldest] <- N[i, oldest] + N[i - 1L, oldest] * exp(-Z[oldest])
  }
  mean <- as.numeric(N)
  catch_mean <- as.numeric(sweep(N, 2, F * -expm1(-Z) / Z, "*"))
  fished <- catch_mean > 0
  catch <- numeric(nrow(grid))
  catch[fished] <- tweedie::rtweedie(sum(fished), mu = catch_mean[fished], phi = .16, power = 1.52)
  draw_delta <- function(probability) {
    observed <- stats::runif(nrow(grid)) < probability
    index <- numeric(nrow(grid))
    index[observed] <- stats::rlnorm(sum(observed),
      log(mean[observed] / probability[observed]) - .3^2 / 2, .3)
    index
  }
  truth <- c(logit = 42, poisson = 54, power = 1.52, phi = 8)
  logit <- draw_delta(mean / (mean + truth[["logit"]]))
  poisson <- draw_delta(-expm1(-log(2) * mean / truth[["poisson"]]))
  tweedie_index <- tweedie::rtweedie(nrow(grid), mu = mean,
    phi = truth[["phi"]], power = truth[["power"]])
  index <- do.call(rbind, lapply(names(truth)[1:3], function(series) {
    value <- switch(series, logit = logit, poisson = poisson, power = tweedie_index)
    transform(grid, obs = value, survey = series, samp_time = 0, unit = "fish per tow")
  }))
  index$obs[index$survey == "power" & index$year == min(years) & index$age == 3] <- NA_real_
  data <- make_tam_data(catch = transform(grid, obs = catch, unit = "fish"), index = index,
    weight = transform(grid, obs = as.numeric(age), unit = "kg"),
    maturity = transform(grid, obs = as.numeric(age >= 4)))
  args <- list(data = data, years = years, ages = ages,
    N_settings = list(process = "off", init = "exp"),
    R_settings = list(process = "off", mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~1, age_breaks = 3:7),
    M_settings = list(process = "off", mu_supplied = ~I(.3)),
    catch_settings = list(family = tweedie(1.52), sd_form = ~0, sd_supplied = ~I(.4)),
    index_settings = list(q_form = ~0, sd_form = ~0 + survey,
      family = list(logit = delta_lognormal(), poisson = delta_lognormal("poisson-link"),
        power = tweedie())))
  list(args = args, truth = truth)
}

test_that("Tweedie preparation explains an unavailable or outdated optional package", {
  require_namespace <- base::requireNamespace
  namespace_version <- base::getNamespaceVersion
  family <- tweedie(power = 1.5)
  survey <- as.character(cod_obs$index$survey[1L])
  settings <- list(
    list(catch_settings = list(family = family)),
    list(catch_settings = list(family = list(Fishery = family))),
    list(index_settings = list(family = family)),
    list(index_settings = list(family = stats::setNames(list(family), survey))))
  args <- list(data = cod_obs, years = 1983:1985, ages = 2:5)
  for (available in c(FALSE, TRUE)) {
    testthat::with_mocked_bindings({
      expect_identical(tweedie(power = 1.5), family)
      default <- do.call(prepare_tam, args)
      expect_true(all(default$obs_map$family == "lognormal"))
      for (setting in settings) {
        error <- tryCatch(do.call(prepare_tam, c(args, setting)), error = identity)
        expect_s3_class(error, "error")
        expect_match(conditionMessage(error), "3[.]2[.]0")
        expect_match(conditionMessage(error), "install.packages\\('tweedie'\\)")
      }
    }, requireNamespace = function(package, ...) {
      if (identical(package, "tweedie")) available else require_namespace(package, ...)
    }, getNamespaceVersion = function(ns) {
      if (identical(ns, "tweedie")) "3.1.0" else namespace_version(ns)
    }, .package = "base")
  }
})

test_that("Tweedie adds zero mass and the log-observation density", {
  skip_if_not_installed("tweedie", minimum_version = "3.2.0")
  case <- tweedie_reference_case()
  args <- case$args
  args$index_settings$family$power <- tweedie(power = case$truth[["power"]])
  dat <- do.call(prepare_tam, args)
  expect_null(make_par(dat)$logit_tweedie_power_index)
  rows <- which(dat$obs_map$type == "index" & dat$obs_map$family == "tweedie" & dat$is_observed)
  without <- dat
  without$is_observed[rows] <- FALSE
  dat$.likelihood_component <- without$.likelihood_component <- "observation"
  par <- make_par(dat)
  tape <- function(d) RTMB::MakeADFun(function(p) nll_fun(p, d), par, silent = TRUE)
  full <- tape(dat)
  excluded <- tape(without)
  args$index_settings$family$power <- tweedie()
  estimated <- do.call(prepare_tam, args)
  expect_identical(names(make_par(estimated)$logit_tweedie_power_index), "power")

  report_dat <- dat
  report_dat$.likelihood_component <- NULL
  report <- tape(report_dat)$report()
  p <- report$tweedie_power[rows]
  expect_equal(p, rep(case$truth[["power"]], length(rows)))
  mu <- exp(report$log_pred[rows])
  phi <- report$sd_obs[rows]^2
  y <- dat$obs_values[rows]
  zeros <- y == 0
  want <- -mu[zeros]^(2 - p[zeros]) / (phi[zeros] * (2 - p[zeros]))
  positive <- tweedie::dtweedie(y[!zeros], mu = mu[!zeros], phi = phi[!zeros],
    power = p[1L], log = TRUE) + log(y[!zeros])
  expect_equal(full$fn() - excluded$fn(), -sum(want, positive), tolerance = 1e-8)
  expect_equal(report$log_p_zero[rows[zeros]], want, tolerance = 1e-12)
  expect_true(all(is.finite(full$gr())))
  scale <- 73
  density <- tweedie::dtweedie(y[!zeros] * scale, mu = mu[!zeros] * scale,
    phi = phi[!zeros] * scale^(2 - p[1L]), power = p[1L], log = TRUE) + log(y[!zeros] * scale)
  expect_equal(density, positive, tolerance = 1e-8)
})

test_that("Tweedie uses sampled zeros for dispersion and excludes structural zeros", {
  skip_if_not_installed("tweedie", minimum_version = "3.2.0")
  case <- tweedie_reference_case()
  args <- case$args
  args$catch_settings <- list(family = tweedie(1.5))
  args$data$catch$obs[args$data$catch$age == 3L] <- 0
  args$data$catch$obs[args$data$catch$year == 1983 & args$data$catch$age == 4L] <- NA_real_
  args$catch_settings$sd_form <- ~ 0 + I(as.numeric(age == 3L))
  dat <- do.call(prepare_tam, args)
  catch <- dat$obs_map$type == "catch"
  structural <- catch & dat$obs_map$age == 2L
  zeros <- catch & dat$obs_map$age == 3L
  missing <- catch & dat$obs_map$year == 1983 & dat$obs_map$age == 4L
  expect_false(any(dat$is_observed[structural | missing]))
  expect_true(all(dat$is_observed[zeros]))
  expect_true(all(dat$log_obs[zeros] == -Inf))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), make_par(dat), silent = TRUE)
  report <- obj$report()
  expect_true(all(report$log_p_zero[structural] == 0))
  expect_true(all(report$log_p_positive[structural] == -Inf))
  expect_true(all(is.finite(obj$gr())))
})

test_that("Tweedie positive medians remain available for rare encounters", {
  skip_if_not_installed("tweedie", minimum_version = "3.2.0")
  power <- 1.5
  log_mean <- log(1e-48)
  rate <- exp((2 - power) * log_mean) / (2 - power)
  log_positive <- log(-expm1(-rate))
  expect_identical(-expm1(log_positive - log(2)), 1)
  median <- .tweedie_log_positive_median(log_mean, 1, power, log_positive)
  expect_true(is.finite(median))
  expect_equal(.tweedie_log_cdf(median, log_mean, 1, power, lower.tail = FALSE),
    log_positive - log(2), tolerance = 1e-8)
})

test_that("mixed zero families recover parameters and retain forecast support", {
  skip_if_not_installed("tweedie", minimum_version = "3.2.0")
  case <- tweedie_reference_case()
  fit <- do.call(fit_tam, c(case$args, list(silent = TRUE)))
  expect_true(fit$is_converged)
  expect_false("tweedie_power_report" %in% names(tidy_sdrep(fit)))
  fixed_catch_power <- subset(tidy_par(fit)$fixed, par == "tweedie_power_catch")
  expect_true(fixed_catch_power$is_fixed)
  expect_equal(fixed_catch_power$est, 1.52)
  index <- subset(fit$obs_pred$index, survey == "power")
  expect_true(all(index$tweedie_power > 1 & index$tweedie_power < 2))
  expect_equal(unique(index$tweedie_power), case$truth[["power"]], tolerance = .16)
  power <- subset(tidy_par(fit)$fixed, par == "tweedie_power_index")
  expect_equal(power$est, unique(index$tweedie_power))
  expect_true(is.finite(power$se) & power$se > 0)
  detection <- subset(tidy_par(fit)$fixed, par == "y50_index")
  expect_equal(detection$est[match(c("logit", "poisson"), detection$coef)],
    unname(case$truth[c("logit", "poisson")]), tolerance = .30)
  expect_equal(index$pred_mean, index$pred)
  expect_equal(index$p_positive * index$pred_positive_mean, index$pred)
  expect_true(all(index$pred_positive_median > 0))
  half_tail <- tweedie::ptweedie(index$pred_positive_median, mu = index$pred,
    phi = index$sd^2, power = index$tweedie_power[1L], lower.tail = FALSE)
  expect_equal(half_tail, index$p_positive / 2, tolerance = 1e-6)
  expect_equal(index$std_res, (index$obs - index$pred) /
    (index$sd * index$pred^(index$tweedie_power / 2)))
  expect_true(all(is.finite(index$quantile_res[index$is_observed])))
  expect_true(all(is.na(index$quantile_res[is.na(index$obs)])))
  expect_error(residuals(fit, type = "osa"), "all-lognormal")

  simulations <- sim_tam(fit, n = 2L, seed = 81, par_uncertainty = "none", progress = FALSE)
  expect_true(all(simulations$catch$obs[simulations$catch$age == 2L] == 0))
  for (series in c("logit", "poisson", "power")) {
    values <- simulations$index$obs[simulations$index$survey == series]
    expect_true(any(values == 0, na.rm = TRUE))
    expect_true(any(values > 0, na.rm = TRUE))
  }
  hind <- fit_hindcast(fit, folds = 1, proj_settings = list(F_mult = 1), progress = FALSE)
  score <- subset(hind$occurrence_scores, level == "series" & type == "index")
  expect_setequal(score$series, c("logit", "poisson", "power"))
  forecasts <- subset(hind$forecast_pairs$index, survey == "power")
  expect_equal(.forecast_log_pred(forecasts), forecasts$log_pred_positive)
})
