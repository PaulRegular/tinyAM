sd_test_dat <- function(type = "iid", component = "catch", ...) {
  form <- switch(type, iid = ~ iid(age, sd = .3), rw = ~ rw(age, sd = .3),
    ar1 = ~ ar1(age, sd = .3, phi = .6))
  args <- list(data = cod_obs, years = 1983:1989, ages = 2:8,
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ q_block, sd_form = ~ 1, fill_missing = FALSE))
  args[[paste0(component, "_settings")]]$sd_form <- form
  do.call(prepare_tam, utils::modifyList(args, list(...)))
}

test_that("SD formulas use compact Gaussian states and unchanged densities", {
  for (component in c("catch", "index")) for (type in c("iid", "rw", "ar1")) {
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
        sum(dnorm(x[-1], .6 * x[-length(x)], .3, log = TRUE)))
    expect_equal(unname(e$nll), expected, tolerance = 1e-10)
    expect_length(x, if (type == "rw") 6L else 7L)
    surface <- if (type == "rw") c(0, x) else x
    expect_equal(e$contribution, surface[match(dat$obs[[component]]$age, dat$ages)])
    expect_true(term$parameter %in% tinyAM:::.formula_random_parameters(dat))
  }
})

test_that("fixed formulas and supplied SD offsets preserve their scale", {
  dat <- sd_test_dat(catch_settings = list(sd_form = ~ 0 + rw(age, sd = .3),
    sd_supplied = ~ I(.2)))
  p <- make_par(dat)
  term <- dat$sd_catch_terms[[1]]
  p[[term$parameter]][] <- .1
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  expect_equal(obj$report()$sd_obs[dat$obs_map$type == "catch"],
    .2 * exp(ifelse(dat$obs$catch$age == 2, 0, .1)))
  ordinary <- sd_test_dat(catch_settings = list(sd_form = ~ age + I(age^2)))
  expect_length(ordinary$sd_catch_terms, 0)
  expect_equal(ordinary$sd_catch_modmat, model.matrix(~ age + I(age^2), ordinary$obs$catch))
  expect_warning(sd_test_dat(catch_settings = list(sd_form = ~ iid(age, sd = .3),
    sd_supplied = ~ I(.2))), "Dropping intercept")
})

test_that("SD simulation uses returned states and matching observation rows", {
  dat <- sd_test_dat("rw", index_settings = list(q_form = ~ q_block,
    sd_form = ~ ar1(age, sd = .3, phi = .6), fill_missing = FALSE))
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
    sd_form = ~ ar1(year, by = survey, sd = .3, phi = .6), fill_missing = FALSE),
    proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1))
  p <- make_par(dat)
  term <- dat$sd_index_terms[[1]]
  p[[term$parameter]][] <- seq_along(p[[term$parameter]]) / 100
  shorter <- sd_test_dat("ar1", "index", years = 1983:1987,
    index_settings = dat$index_settings)
  merged <- tinyAM:::.merge_start_par(make_par(shorter), p)
  expect_equal(merged[[term$parameter]], p[[term$parameter]][names(merged[[term$parameter]])])
  proxy <- list(dat = dat, parameter_values = p, rep = list(),
    obj = list(env = list(parList = function(par) par)))
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
        exp(log(.3) + e$contribution[rows]), log = TRUE))
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
  expect_error(sd_test_dat(data = obs, catch_settings = list(sd_form = ~ iid(age) + iid(age_copy))),
    "same IID variance")
  # Unlike a Gaussian random mean, one random log-SD per row is a scale mixture.
  obs$catch$age_group <- factor(obs$catch$age)
  dat <- sd_test_dat(data = obs, catch_settings = list(sd_form = ~ iid(year, by = age_group)))
  expect_true("SD_effect_replication" %in% tinyAM:::.sd_process_advisories(dat)$issue)
  well <- sd_test_dat()
  expect_equal(nrow(tinyAM:::.sd_process_advisories(well)), 0L)
  thin <- sd_test_dat(ages = 2:5, catch_settings = list(sd_form = ~ ar1(age)))
  expect_true("SD_effect_support" %in% tinyAM:::.sd_process_advisories(thin)$issue)
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
