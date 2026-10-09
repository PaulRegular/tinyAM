effect_data <- function() {
  d <- expand.grid(age = 1:3, year = 2000:2005, survey = c("A", "B"))
  d$survey <- factor(d$survey)
  d$obs <- 1
  d$is_proj <- FALSE
  d
}

effect_dat <- function(formula, data = effect_data()) {
  design <- tinyAM:::.parse_q_formula(formula, data)
  c(list(obs = list(index = data), index_settings = list(q_link = "log"),
         sd_index_modmat = model.matrix(~ 1, data)), design)
}

test_that("formula compilation preserves ordinary and mono designs", {
  d <- effect_data()
  expect_identical(tinyAM:::.parse_q_formula(~ survey * age, d)$q_modmat,
                   model.matrix(~ survey * age, d))
  design <- tinyAM:::.parse_q_formula(~ survey + mono(age, by = survey) + iid(year), d)
  baseline <- tinyAM:::.parse_q_formula(~ survey + mono(age, by = survey), d)
  expect_equal(design$q_modmat, baseline$q_modmat)
  expect_equal(design$q_mono_modmat, baseline$q_mono_modmat)
  expect_named(design$q_terms, "iid_year")
  expect_equal(tinyAM:::.parse_q_formula(~ (1 | survey), d)$q_terms,
               tinyAM:::.parse_q_formula(~ iid(survey), d)$q_terms)
  expect_equal(tinyAM:::.parse_q_formula(~ tinyAM::ar1(year), d)$q_terms,
               tinyAM:::.parse_q_formula(~ ar1(year), d)$q_terms)
})

test_that("term states are unique, named, ordered and preserve gaps", {
  d <- effect_data()
  d <- d[d$year != 2002, ]
  d$is_proj <- d$year == 2005
  dat <- effect_dat(~ survey + rw(year, by = survey), d)
  term <- dat$q_terms[[1]]
  expect_equal(term$groups[[1]]$levels, 2000:2005)
  expect_equal(term$groups[[1]]$states, paste("A", 2001:2005, sep = ":"))
  expect_equal(term$groups[[1]]$is_proj, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  p <- tinyAM:::.q_term_parameters(dat$q_terms)
  expect_equal(length(p[[term$parameter]]), 10)
  expect_false(anyDuplicated(names(p[[term$parameter]])) > 0)
  expect_equal(length(p[[term$sd_parameter]]), 1)
  expect_equal(tinyAM:::.q_random_parameters(dat), term$parameter)
  later <- effect_dat(~ survey + rw(year, by = survey), d[d$year != 2004, ])
  expect_identical(names(later$q_terms), names(dat$q_terms))
  ordered <- d
  ordered$step <- ordered(ordered$year, levels = 2000:2005)
  expect_equal(effect_dat(~ rw(step), ordered)$q_terms[[1]]$groups[[1]]$levels,
               as.character(2000:2005))
})

test_that("numeric by multiplies a shared process and categorical by separates groups", {
  d <- effect_data()
  d$multiplier <- ifelse(d$survey == "A", 1, -2)
  dat <- effect_dat(~ iid(year, by = multiplier, sd = .2), d)
  term <- dat$q_terms[[1]]
  expect_length(term$groups, 1)
  p <- tinyAM:::.q_term_parameters(dat$q_terms)
  p[[term$parameter]][] <- .3
  effect <- tinyAM:::.q_effects(p, dat)
  expect_equal(effect$contribution, .3 * d$multiplier)
  expect_equal(effect$nll, -sum(dnorm(rep(.3, 6), 0, .2, log = TRUE)), tolerance = 1e-10)
  d$multiplier <- factor(d$multiplier)
  expect_length(effect_dat(~ iid(year, by = multiplier), d)$q_terms[[1]]$groups, 2)
})

test_that("IID RW and AR1 densities use unique states and exact normal densities", {
  for (type in c("iid", "rw", "ar1")) {
    form <- switch(type, iid = ~ iid(year, sd = .2),
      rw = ~ rw(year, sd = .2), ar1 = ~ ar1(year, sd = .2, phi = .6))
    dat <- effect_dat(form)
    p <- tinyAM:::.q_term_parameters(dat$q_terms)
    nm <- dat$q_terms[[1]]$parameter
    p[[nm]][] <- seq(.1, by = .1, length.out = length(p[[nm]]))
    x <- p[[nm]]
    expected <- switch(type,
      iid = -sum(dnorm(x, 0, .2, log = TRUE)),
      rw = -sum(dnorm(diff(c(0, x)), 0, .2, log = TRUE)),
      ar1 = -dnorm(x[1], 0, .2 / sqrt(1 - .6^2), log = TRUE) -
        sum(dnorm(x[-1], .6 * x[-length(x)], .2, log = TRUE)))
    effect <- tinyAM:::.q_effects(p, dat)
    expect_equal(unname(effect$nll), unname(expected), tolerance = 1e-10)
    expect_equal(effect$contribution[1:3], rep(if (type == "rw") 0 else .1, 3))
    obj <- RTMB::MakeADFun(function(p) tinyAM:::.q_effects(p, dat)$nll, p, silent = TRUE)
    expect_equal(obj$fn(obj$par), unname(expected), tolerance = 1e-10)
    numerical <- vapply(seq_along(obj$par), function(i) {
      hi <- lo <- obj$par
      hi[i] <- hi[i] + 1e-5
      lo[i] <- lo[i] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, numeric(1))
    expect_equal(as.numeric(obj$gr(obj$par)), numerical, tolerance = 1e-6)
  }
})

test_that("logistic multiplies q with positive slope and group-specific curves", {
  dat <- effect_dat(~ survey + logistic(age, by = survey))
  p <- tinyAM:::.q_term_parameters(dat$q_terms)
  term <- dat$q_terms[[1]]
  p[[paste0("q_a50_", term$id)]][] <- c(2, 1)
  p[[paste0("log_q_slope_", term$id)]][] <- log(c(1, 2))
  effect <- tinyAM:::.q_effects(p, dat)
  d <- dat$obs$index
  expected <- plogis(ifelse(d$survey == "A", d$age - 2, 2 * (d$age - 1)))
  expect_equal(exp(effect$log_selectivity), expected)
  expect_equal(effect$nll, 0)
  expect_true(all(.8 * exp(effect$log_selectivity) < 1))
  obj <- RTMB::MakeADFun(function(p) sum(tinyAM:::.q_effects(p, dat)$log_selectivity^2),
                        p, silent = TRUE)
  expect_true(all(is.finite(obj$gr(obj$par))))
})

test_that("unsupported terms and exact variance aliases are rejected", {
  d <- effect_data()
  for (formula in list(~ rw(age / 2), ~ iid(year):survey, ~ (age | survey),
                       ~ iid(year) + iid(year))) {
    expect_error(effect_dat(formula, d))
  }
  expect_error(effect_dat(~ logistic(age, by = year)), "categorical")
  expect_error(effect_dat(~ ar1(survey)), "ordered")
  expect_error(effect_dat(~ iid(year, sd = 0)), "positive")
  expect_error(effect_dat(~ ar1(year, phi = 1)), "below")
  expect_error(tinyAM:::.check_q_terms(effect_dat(~ survey + iid(survey))), "saturated")
  d$unique_row <- seq_len(nrow(d))
  expect_error(tinyAM:::.check_q_terms(effect_dat(~ iid(unique_row), d)), "observation error")
  expect_no_error(tinyAM:::.check_q_terms(effect_dat(~ iid(unique_row, sd = .2), d)))
  expect_no_error(tinyAM:::.check_q_terms(effect_dat(~ survey + iid(year, by = survey), d)))
  d$year_copy <- d$year
  expect_error(tinyAM:::.check_q_terms(effect_dat(~ iid(year) + iid(year_copy), d)),
               "same IID variance")
  expect_error(tinyAM:::.check_q_terms(effect_dat(~ factor(age) + logistic(age), d)),
               "redundant curve")
  expect_no_error(tinyAM:::.check_q_terms(effect_dat(~ survey + logistic(age, by = survey), d)))
})

test_that("process simulation is reproducible and agrees with the density convention", {
  for (type in c("iid", "rw", "ar1")) {
    form <- switch(type, iid = ~ iid(year, sd = .2), rw = ~ rw(year, sd = .2),
                   ar1 = ~ ar1(year, sd = .2, phi = .6))
    dat <- effect_dat(form)
    p <- tinyAM:::.q_term_parameters(dat$q_terms)
    set.seed(703)
    a <- tinyAM:::.q_effects(p, dat, simulate = TRUE)
    set.seed(703)
    expect_equal(a, tinyAM:::.q_effects(p, dat, simulate = TRUE))
    p[names(a$parameters)] <- a$parameters
    expect_equal(tinyAM:::.q_effects(p, dat)$contribution, a$contribution)
    draws <- replicate(2000, tinyAM:::.q_effects(p, dat, simulate = TRUE)$parameters[[1]])
    innovations <- switch(type, iid = draws,
      rw = rbind(draws[1, ], draws[-1, ] - draws[-nrow(draws), ]),
      ar1 = draws[-1, ] - .6 * draws[-nrow(draws), ])
    expect_lt(abs(mean(innovations)), .006)
    expect_equal(sd(as.numeric(innovations)), .2, tolerance = .01)
    if (type == "ar1") expect_equal(sd(draws[1, ]), .2 / sqrt(1 - .6^2), tolerance = .02)
  }
})

test_that("forecast states integrate out and have the intended conditional means", {
  for (type in c("iid", "rw", "ar1")) {
    form <- switch(type, iid = ~ iid(year, sd = .2), rw = ~ rw(year, sd = .2),
                   ar1 = ~ ar1(year, sd = .2, phi = .6))
    d <- effect_data()
    d <- d[d$survey == "A", ]
    d$is_proj <- d$year >= 2004
    dat <- effect_dat(form, d)
    past <- effect_dat(form, d[!d$is_proj, ])
    observation <- c(.1, .2, -.1, .4)
    fit_effect <- function(dat) {
      p <- tinyAM:::.q_term_parameters(dat$q_terms)
      RTMB::MakeADFun(function(p) {
        e <- tinyAM:::.q_effects(p, dat)
        e$nll - sum(RTMB::dnorm(observation,
          e$contribution[c(1, 4, 7, 10)], .1, log = TRUE))
      }, p, random = tinyAM:::.q_random_parameters(dat), silent = TRUE)
    }
    obj <- fit_effect(dat)
    historical <- fit_effect(past)
    expect_equal(obj$fn(obj$par), historical$fn(historical$par), tolerance = 1e-10)
    p <- obj$env$parList()
    e <- tinyAM:::.q_effects(p, dat)$contribution[c(10, 13, 16)]
    expect_equal(unname(e[-1]), switch(type, iid = c(0, 0),
      rw = rep(unname(e[1]), 2), ar1 = unname(e[1]) * c(.6, .6^2)), tolerance = 1e-7)
  }
})

test_that("new q parameters are registered and simulations reuse the returned states", {
  dat <- make_test_dat(index_settings = list(sd_form = ~1,
    q_form = ~ q_block + iid(year, sd = .2), fill_missing = FALSE))
  p <- make_par(dat)
  term <- dat$q_terms[[1]]
  expect_true(term$parameter %in% names(p))
  expect_false(any(startsWith(names(p), "log_sd_q_")))
  set.seed(307)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  expect_true(term$parameter %in% names(simulated))
  p[names(simulated)[names(simulated) %in% names(p)]] <- simulated[names(simulated) %in% names(p)]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  report <- obj$report()
  expected <- drop(dat$q_modmat %*% p$log_q) + tinyAM:::.q_effects(p, dat)$contribution
  expect_equal(report$log_q_obs, expected)
  error <- simulated$log_obs[dat$is_observed] - report$log_pred[dat$is_observed]
  expect_lt(abs(mean(error / report$sd_obs[dat$is_observed])), .1)
  expect_equal(sd(error / report$sd_obs[dat$is_observed]), 1, tolerance = .1)
})

test_that("formula processes fit and warm starts align named states", {
  fit <- update(default_fit, silent = TRUE,
    index_settings = list(q_form = ~ q_block + ar1(year, sd = .15, phi = .6),
                          sd_form = ~1, fill_missing = FALSE),
    start_par = as.list(default_fit$sdrep, "Estimate"))
  term <- fit$dat$q_terms[[1]]
  expect_true(term$parameter %in% fit$obj$env$.random)
  expect_true(is.finite(fit$opt$objective))
  expect_true(all(is.finite(fit$rep$log_q_obs)))
  estimates <- as.list(fit$sdrep, "Estimate")
  shorter <- prepare_tam(cod_obs, years = 1983:2020, ages = AGES,
    index_settings = fit$dat$index_settings)
  base <- make_par(shorter)
  merged <- tinyAM:::.merge_start_par(base, estimates)
  labels <- names(base[[term$parameter]])
  expect_equal(merged[[term$parameter]], estimates[[term$parameter]][labels])
  expect_true(all(grepl("198[3-9]|199[0-9]|200[0-9]|201[0-9]|2020", labels)))
  set.seed(205)
  expect_no_error(tinyAM:::.sim_obs(fit, tinyAM:::.draw_none, redraw_random = TRUE))
  effects <- fit$formula_effects$levels[[term$id]]
  expect_equal(effects$est, unname(estimates[[term$parameter]]))
  expect_true(all(effects$se_scale == "reported"))
  expect_equal(effects$coordinate, YEARS)
  expect_true(all(is.finite(fit$obs_pred$index$q_se)))
  expect_equal(fit$obs_pred$index$q_se,
    as.numeric(as.list(fit$sdrep, "Std. Error", report = TRUE)$q_link_prediction))
  tabs <- tidy_tam(model_list = list(process = fit, ordinary = default_fit), interval = .9)
  expect_equal(unique(tabs$formula_effects$levels[[term$id]]$model), "process")
  interval <- subset(tabs$obs_pred$index, model == "process")
  expect_true(all(interval$q_upr - interval$q_lwr <
    fit$obs_pred$index$q_upr - fit$obs_pred$index$q_lwr))
})

test_that("RW increments and numeric-by contributions retain correct signed uncertainty", {
  dat <- make_test_dat(index_settings = list(sd_form = ~1,
    q_form = ~ q_block + rw(year, by = age, sd = .2), fill_missing = FALSE))
  p <- make_par(dat)
  term <- dat$q_terms[[1]]
  p[[term$parameter]][] <- -.2
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  # Supply known local curvature without fitting the population again.
  summaries <- p
  errors <- lapply(p, function(x) { x[] <- .1; x })
  fit <- structure(list(dat = dat, obj = obj, rep = obj$report(),
    parameter_values = obj$env$last.par), class = "tam_fit")
  local_mocked_bindings(.tam_parameter_summary = function(fit, what) {
    if (what == "Estimate") summaries else errors
  }, .package = "tinyAM")
  effects <- tinyAM:::.tidy_q_effects(fit)
  levels <- effects$levels[[term$id]]
  expect_equal(levels$est, c(0, rep(-.2, length(YEARS) - 1L)))
  expect_equal(levels$se, c(0, rep(.1, length(YEARS) - 1L)))
  expect_true(levels$is_fixed[1])
  increments <- effects$increments[[term$id]]
  expect_equal(increments$est, c(-.2, rep(0, length(YEARS) - 2L)))
  expect_true(all(is.na(increments$se)))
  contribution <- effects$contributions[[term$id]]
  expect_equal(contribution$est[contribution$year > min(YEARS)],
    -.2 * contribution$age[contribution$year > min(YEARS)])
  expect_equal(contribution$se[contribution$year > min(YEARS)],
    .1 * contribution$age[contribution$year > min(YEARS)])
  expect_false("eta_q_increments" %in% names(tidy_rep(fit)))
  plot <- plotly::plotly_build(tinyAM:::.plot_formula_effect(levels))
  expect_equal(plot$x$layout$yaxis$type, "linear")
  expect_lt(plot$x$layout$yaxis$range[1], -.2)
})

test_that("logit catchability intervals include logistic curves and remain bounded", {
  dat <- make_test_dat(index_settings = list(sd_form = ~1,
    q_form = ~ logistic(age), q_link = "logit", fill_missing = FALSE))
  p <- make_par(dat)
  obj <- RTMB::MakeADFun(function(p) {
    nll_fun(p, dat)
    Reduce(`+`, lapply(p, function(x) sum(x^2))) / (2 * .1^2)
  }, p, silent = TRUE)
  fit <- structure(list(dat = dat, obj = obj, rep = obj$report(),
    sdrep = RTMB::sdreport(obj, getReportCovariance = FALSE)), class = "tam_fit")
  tab <- tidy_obs_pred(fit)$index
  expect_true(all(tab$q_lwr > 0 & tab$q_upr < 1))
  expect_true(all(tab$q_se_scale == "logit"))
  z <- qnorm(.975)
  values <- as.list(fit$sdrep, "Estimate", report = TRUE)$q_link_prediction
  ses <- as.list(fit$sdrep, "Std. Error", report = TRUE)$q_link_prediction
  expect_equal(tab$q_lwr, as.numeric(plogis(values - z * ses)))
  expect_equal(tab$q, as.numeric(plogis(values)))
})
