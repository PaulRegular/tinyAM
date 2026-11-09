sd_test_dat <- function(type = "iid", component = "catch", ...) {
  form <- switch(type, iid = ~ iid(age, sd = .3), rw = ~ rw(age, sd = .3),
    ar1 = ~ ar1(age, sd = .3, phi = .6))
  args <- list(data = cod_obs, years = 1983:1989, ages = 2:8,
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(q_form = ~ q_block, sd_form = ~ 1))
  args[[paste0(component, "_settings")]]$sd_form <- form
  do.call(prepare_tam, utils::modifyList(args, list(...)))
}

test_that("SD formulas use compact Gaussian states and unchanged densities", {
  for (component in c("catch", "index")) {
    for (type in c("iid", "rw", "ar1")) {
      dat <- sd_test_dat(type, component)
      p <- make_par(dat)
      term <- dat[[paste0("sd_", component, "_terms")]][[1]]
      x <- seq(-.2, .3, length.out = length(p[[term$parameter]]))
      p[[term$parameter]][] <- x
      e <- tinyAM:::.sd_effects(p, dat, component)
      expected <- switch(type,
        iid = -sum(dnorm(x, 0, .3, log = TRUE)),
        rw = -sum(dnorm(c(x[1], diff(x)), 0, .3, log = TRUE)),
        ar1 = -dnorm(x[1], 0, .3 / sqrt(1 - .6^2), log = TRUE) -
          sum(dnorm(x[-1], .6 * x[-length(x)], .3, log = TRUE))
      )
      expect_equal(unname(e$nll), expected, tolerance = 1e-10)
      expect_length(x, if (type == "rw") 6L else 7L)
      surface <- if (type == "rw") c(0, x) else x
      expect_equal(e$contribution, surface[match(dat$obs[[component]]$age, dat$ages)])
      expect_true(term$parameter %in% tinyAM:::.formula_random_parameters(dat))
    }
  }
})

test_that("fixed formulas and supplied SD offsets preserve their scale", {
  dat <- sd_test_dat(catch_settings = list(
    sd_form = ~ 0 + rw(age, sd = .3),
    sd_supplied = ~ I(.2)
  ))
  p <- make_par(dat)
  term <- dat$sd_catch_terms[[1]]
  p[[term$parameter]][] <- .1
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  expect_equal(
    obj$report()$sd_obs[dat$obs_map$type == "catch"],
    .2 * exp(ifelse(dat$obs$catch$age == 2, 0, .1))
  )
  ordinary <- sd_test_dat(catch_settings = list(sd_form = ~ age + I(age^2)))
  expect_length(ordinary$sd_catch_terms, 0)
  expect_equal(ordinary$sd_catch_modmat, model.matrix(~ age + I(age^2), ordinary$obs$catch))
  expect_warning(sd_test_dat(catch_settings = list(
    sd_form = ~ iid(age, sd = .3),
    sd_supplied = ~ I(.2)
  )), "Dropping intercept")
})

test_that("random SD changes spread and its density, not the observation mean", {
  dat <- sd_test_dat()
  p <- make_par(dat)
  term <- dat$sd_catch_terms[[1]]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  before <- obj$fn(obj$par)
  initial <- obj$report()
  modified <- obj$par
  indices <- names(modified) == term$parameter
  modified[indices] <- seq(-.3, .3, length.out = sum(indices))
  after <- obj$fn(modified)
  final <- obj$report()
  expect_equal(final$log_pred, initial$log_pred)
  rows <- dat$is_observed
  observation_change <- sum(dnorm(dat$observed, initial$log_pred[rows], initial$sd_obs[rows], log = TRUE)) -
    sum(dnorm(dat$observed, final$log_pred[rows], final$sd_obs[rows], log = TRUE))
  process_change <- sum(dnorm(rep(0, sum(indices)), 0, .3, log = TRUE)) -
    sum(dnorm(modified[indices], 0, .3, log = TRUE))
  expect_equal(after - before, observation_change + process_change, tolerance = 1e-8)
})

test_that("SD simulation uses returned states and matching observation rows", {
  dat <- sd_test_dat("rw", index_settings = list(q_form = ~ q_block,
    sd_form = ~ ar1(age, sd = .3, phi = .6)))
  p <- make_par(dat)
  p$log_sd_catch[] <- log(.2)
  p$log_sd_index[] <- log(.4)
  p$log_sd_f <- p$log_sd_r <- log(.1)
  set.seed(754)
  draw <- nll_fun(p, dat, simulate = TRUE)
  set.seed(754)
  expect_equal(draw, nll_fun(p, dat, simulate = TRUE))
  for (term in tinyAM:::.formula_terms(dat)) expect_true(term$parameter %in% names(draw))
  p[intersect(names(p), names(draw))] <- draw[intersect(names(p), names(draw))]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  report <- obj$report()
  for (component in c("catch", "index")) {
    expected <- exp(p[[paste0("log_sd_", component)]] +
      tinyAM:::.sd_effects(p, dat, component)$contribution)
    expect_equal(report$sd_obs[dat$obs_map$type == component], expected)
  }
  residual <- (draw$log_obs[dat$is_observed] - report$log_pred[dat$is_observed]) /
    report$sd_obs[dat$is_observed]
  expect_lt(abs(mean(residual)) * sqrt(length(residual)), 3)
  expect_equal(sd(residual), 1, tolerance = .2)
  analytic <- obj$gr(obj$par)
  target <- grep("eta_sd_", names(obj$par))[1]
  hi <- lo <- obj$par
  hi[target] <- hi[target] + 1e-5
  lo[target] <- lo[target] - 1e-5
  expect_true(all(is.finite(analytic)))
  expect_equal(analytic[target], (obj$fn(hi) - obj$fn(lo)) / 2e-5, tolerance = 1e-6)
})

test_that("grouping, warm starts, forecasts and uncertainty retain SD semantics", {
  dat <- sd_test_dat("ar1", "index", index_settings = list(q_form = ~ q_block,
    sd_form = ~ ar1(year, by = survey, sd = .3, phi = .6)),
    proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1))
  p <- make_par(dat)
  term <- dat$sd_index_terms[[1]]
  p[[term$parameter]][] <- seq_along(p[[term$parameter]]) / 100
  shorter <- sd_test_dat("ar1", "index",
    years = 1983:1987,
    index_settings = dat$index_settings
  )
  merged <- tinyAM:::.merge_start_par(make_par(shorter), p)
  expect_equal(merged[[term$parameter]], p[[term$parameter]][names(merged[[term$parameter]])])
  proxy <- list(
    dat = dat, parameter_values = p, rep = list(),
    obj = list(env = list(parList = function(par) par))
  )
  tables <- tinyAM:::.tidy_formula_effects(proxy)
  tab <- tables$levels[[term$id]]
  expect_true(all(tab$component == "sd_index"))
  expect_equal(tab$est, unname(p[[term$parameter]]))
  expect_true(all(is.na(tab$se)))
  expect_true(any(tab$is_proj))
  expect_match(tinyAM:::.formula_effect_label(tab), "Survey SD AR1")
  # Future states have proper densities, so their integral is one.
  objective <- function(d) {
    pp <- make_par(d)[term$parameter]
    RTMB::MakeADFun(function(pp) {
      e <- tinyAM:::.sd_effects(pp, d, "index")
      rows <- !d$obs$index$is_proj
      e$nll - sum(RTMB::dnorm(rep(.2, sum(rows)), 0,
        exp(log(.3) + e$contribution[rows]),
        log = TRUE
      ))
    }, pp, random = term$parameter, silent = TRUE)
  }
  historical <- sd_test_dat("ar1", "index", index_settings = dat$index_settings)
  a <- objective(dat)
  b <- objective(historical)
  expect_equal(a$fn(a$par), b$fn(b$par), tolerance = 1e-8)
})

test_that("SD designs reject redundancies but do not alias mean and scale variance", {
  expect_error(sd_test_dat(catch_settings = list(sd_form = ~ logistic(age))), "not catchability curves")
  expect_error(sd_test_dat(catch_settings = list(sd_form = ~ mono(age))), "not catchability curves")
  expect_error(sd_test_dat(catch_settings = list(sd_form = ~ factor(age) + iid(age))), "saturated fixed")
  expect_error(sd_test_dat(catch_settings = list(sd_form = ~ rw(age) + ar1(age))), "one RW/AR1")
  obs <- cod_obs
  obs$catch$age_copy <- obs$catch$age
  expect_error(
    sd_test_dat(data = obs, catch_settings = list(sd_form = ~ iid(age) + iid(age_copy))),
    "same IID variance"
  )
  # Unlike a Gaussian random mean, one random log-SD per row is a scale mixture.
  obs$catch$age_group <- factor(obs$catch$age)
  dat <- sd_test_dat(data = obs, catch_settings = list(sd_form = ~ iid(year, by = age_group)))
  expect_true("SD_effect_replication" %in% tinyAM:::.sd_process_advisories(dat)$issue)
  well <- sd_test_dat()
  expect_equal(nrow(tinyAM:::.sd_process_advisories(well)), 0L)
  thin <- sd_test_dat(ages = 2:5, catch_settings = list(sd_form = ~ ar1(age)))
  expect_true("SD_effect_support" %in% tinyAM:::.sd_process_advisories(thin)$issue)
  fragile <- sd_test_dat(catch_settings = list(sd_form = ~ rw(age)))
  advisory <- tinyAM:::.sd_process_advisories(fragile)
  expect_true("catch_SD_RW" %in% advisory$issue)
  expect_match(advisory$detail[advisory$issue == "catch_SD_RW"], "when F is also estimated")
  expect_false("catch_SD_RW" %in% tinyAM:::.sd_process_advisories(sd_test_dat("rw"))$issue)
  expect_false("catch_SD_RW" %in% tinyAM:::.sd_process_advisories(
    sd_test_dat("rw", "index", index_settings = list(sd_form = ~ rw(age)))
  )$issue)
})

test_that("SD random intercepts and by specifications use the same formula rules", {
  d <- expand.grid(age = 1:6, year = 2000:2004, survey = c("early", "late"))
  d$obs <- 1
  d$is_proj <- FALSE
  d$multiplier <- ifelse(d$survey == "early", -2, 0)
  a <- tinyAM:::.parse_sd_formula(~ (1 | age), d, "catch")
  b <- tinyAM:::.parse_sd_formula(~ iid(age), d, "catch")
  expect_equal(a, b)
  compiled <- tinyAM:::.parse_sd_formula(~ rw(age, by = survey, sd = .2), d, "index")
  expect_length(compiled$terms[[1]]$groups, 2L)
  compiled <- tinyAM:::.parse_sd_formula(~ iid(age, by = multiplier, sd = .2), d, "index")
  dat <- list(obs = list(index = d), sd_index_terms = compiled$terms)
  p <- tinyAM:::.q_term_parameters(compiled$terms)
  p[[compiled$terms[[1]]$parameter]][] <- .1
  expect_equal(tinyAM:::.sd_effects(p, dat, "index")$contribution, .1 * d$multiplier)
})

test_that("smooth observation SDs change spread with Gaussian penalties", {
  expect_warning(
    dat <- sd_test_dat(
      catch_settings = list(sd_form = ~ s(age, k = 4), sd_supplied = ~ I(.2)),
      index_settings = list(sd_form = ~ t2(age, year, k = c(4, 4))),
      proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1)
    ),
    "Dropping intercept"
  )
  expect_false("(Intercept)" %in% colnames(dat$sd_catch_modmat))
  expect_length(dat$sd_catch_terms, 1L)
  expect_gt(length(dat$sd_index_terms), 1L)
  p <- make_par(dat)
  p$log_sd_index["(Intercept)"] <- log(.4)
  p$log_sd_f <- p$log_sd_r <- log(.1)
  objective <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  initial <- objective$report()

  for (term in c(dat$sd_catch_terms, dat$sd_index_terms)) {
    p[[term$parameter]][] <- seq(-.05, .06, length.out = length(p[[term$parameter]]))
    expect_true(term$parameter %in% tinyAM:::.formula_random_parameters(dat))
  }
  p$log_sd_catch[] <- -.08
  p$log_sd_index[names(p$log_sd_index) != "(Intercept)"] <- -.03
  modified <- objective$par
  for (name in unique(names(modified))) {
    modified[names(modified) == name] <- p[[name]]
  }
  expect_true(is.finite(objective$fn(modified)))
  expect_true(all(is.finite(objective$gr(modified))))
  final <- objective$report()
  expect_equal(final$log_pred, initial$log_pred)

  for (component in c("catch", "index")) {
    field <- paste0("sd_", component)
    terms <- dat[[paste0(field, "_terms")]]
    fitted <- p[[paste0("log_", field)]]
    z <- Reduce("+", lapply(terms, function(term) {
      drop(term$design %*% p[[term$parameter]])
    }))
    penalty <- sum(vapply(terms, function(term) {
      -sum(dnorm(p[[term$parameter]], sd = exp(p[[term$sd_parameter]]), log = TRUE))
    }, numeric(1)))
    effects <- tinyAM:::.sd_effects(p, dat, component)
    expect_equal(effects$nll, penalty, tolerance = 1e-10)
    expect_equal(effects$contribution, z)
    smooth <- terms[[1L]]$smooth_fixed
    smooth_value <- drop(smooth %*% fitted[colnames(smooth)]) + z
    expect_equal(as.numeric(final[[paste0("eta_", field, "_smooth")]]), smooth_value)
    offset <- if (component == "catch") .2 else .4
    expect_equal(final$sd_obs[dat$obs_map$type == component], offset * exp(smooth_value))
  }

  index <- dat$obs$index
  key <- function(d) paste(d$survey, d$age)
  terminal <- which(index$year == max(index$year[!index$is_proj]))
  matched <- terminal[match(key(index[index$is_proj, ]), key(index[terminal, ]))]
  expect_true(any(index$is_proj))
  expect_false(anyNA(matched))
  expect_equal(as.numeric(final$eta_sd_index_smooth[index$is_proj, ]),
    as.numeric(final$eta_sd_index_smooth[matched, ]))

  set.seed(912)
  draw <- nll_fun(p, dat, simulate = TRUE)
  set.seed(912)
  expect_equal(draw, nll_fun(p, dat, simulate = TRUE))
  for (component in c("catch", "index")) {
    terms <- dat[[paste0("sd_", component, "_terms")]]
    expect_true(all(vapply(terms, function(term) {
      !is.null(draw[[term$parameter]]) && any(draw[[term$parameter]] != p[[term$parameter]])
    }, logical(1))))
  }

  proxy <- list(dat = dat, parameter_values = p, rep = final,
    obj = list(env = list(parList = function(par) par)))
  effects <- tinyAM:::.tidy_formula_effects(proxy)$smooths
  expect_length(effects, 2L)
  for (tab in effects) {
    expect_identical(unique(tab$unit), "log SD contribution")
    expect_equal(tab$est, as.numeric(final[[paste0("eta_", tab$component[1L], "_smooth")]]))
  }
})

test_that("survey SD warm starts preserve smooths across basis sign changes", {
  dat <- sd_test_dat(index_settings = list(sd_form = ~ s(age, k = 4)))
  old <- dat
  term <- dat$sd_index_terms[[1L]]
  columns <- colnames(term$smooth_fixed)
  old$sd_index_modmat[, columns] <- -old$sd_index_modmat[, columns]
  reversed <- term
  reversed$smooth_fixed <- -term$smooth_fixed
  reversed$design <- -term$design
  reversed$basis_signature$design$fixed <- -term$basis_signature$design$fixed
  reversed$basis_signature$design$random[[1L]] <- -term$basis_signature$design$random[[1L]]
  old$sd_index_terms[[1L]] <- reversed
  p <- make_par(old)
  p$log_sd_index[columns] <- -.03
  p[[term$parameter]][] <- seq_along(p[[term$parameter]]) / 50
  p[[term$sd_parameter]][] <- log(.3)
  source <- structure(list(dat = old, parameter_values = p,
    obj = list(env = list(parList = function(par) par))), class = "tam_fit")
  start <- tinyAM:::.start_parameters(source, dat)
  merged <- tinyAM:::.merge_start_par(make_par(dat), start)
  before <- tinyAM:::.sd_effects(p, old, "index")
  after <- tinyAM:::.sd_effects(merged, dat, "index")
  expect_equal(drop(dat$sd_index_modmat %*% merged$log_sd_index) + after$contribution,
    drop(old$sd_index_modmat %*% p$log_sd_index) + before$contribution)
  expect_equal(after$nll, before$nll)
  expect_equal(merged[[term$sd_parameter]], p[[term$sd_parameter]])
})
