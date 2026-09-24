mono_values <- function(design, beta, theta) {
  drop(design$q_modmat %*% beta + design$q_mono_modmat %*% exp(theta))
}

test_that("monotonic q uses cumulative positive log-q steps", {
  d <- data.frame(x = 1:4)
  design <- tinyAM:::.parse_q_formula(~mono(x), d)
  expect_equal(unname(mono_values(design, -2, log(c(.1, .2, .3)))),
               c(-2, -1.9, -1.7, -1.4))
  set.seed(823)
  for (i in 1:10) {
    log_q <- mono_values(design, -2, runif(3, -5, 1))
    expect_true(all(diff(log_q) > 0))
    expect_true(all(diff(exp(log_q)) > 0))
  }
})

test_that("numeric order, declared factor order and pooled plateaus are retained", {
  for (x in list(c(10, 2, 5, 5), factor(c("old", "young", "mid", "mid"),
                  levels = c("young", "mid", "old")),
                ordered(c("old", "young", "mid", "mid"),
                  levels = c("young", "mid", "old")))) {
    design <- tinyAM:::.parse_q_formula(~ mono(x), data.frame(x = x))
    expect_equal(unname(mono_values(design, -2, log(c(.1, .2)))), c(-1.7, -2, -1.9, -1.9))
  }
})

test_that("surveys have separate baselines and steps with different level counts", {
  d <- data.frame(survey = factor(c("B", "A", "A", "B", "A", "B", "A")),
                  x = factor(c(4, 2, 4, 1, 1, 3, 3), levels = 1:4))
  design <- tinyAM:::.parse_q_formula(~ survey + mono(x, by = survey), d)
  expect_equal(ncol(design$q_modmat), 2L)
  expect_equal(ncol(design$q_mono_modmat), 5L)
  expect_equal(design$q_mono_steps$by_level, c(rep("A", 3), rep("B", 2)))
  theta <- log(c(.1, .2, .3, .4, .5))
  q <- mono_values(design, c(-2, 1), theta)
  expect_equal(unname(q), c(-.1, -1.9, -1.4, -1, -2, -.6, -1.7))
  theta[1] <- log(.8)
  changed <- mono_values(design, c(-2, 1), theta)
  expect_equal(changed[d$survey == "B"], q[d$survey == "B"])
  expect_gt(changed[2], q[2])
  shared <- tinyAM:::.parse_q_formula(~ mono(x, by = survey), d)
  expect_equal(ncol(shared$q_modmat), 1L)
  expect_equal(shared$q_mono_modmat, design$q_mono_modmat)
  no_intercept <- tinyAM:::.parse_q_formula(~ mono(x) - 1, d)
  expect_equal(ncol(no_intercept$q_modmat), 0L)
})

test_that("malformed or ambiguous monotonic formulas fail clearly", {
  d <- data.frame(x = 1:4, survey = c("A", "A", "B", "B"))
  for (f in list(~mono(), ~mono(missing), ~mono(x, by = missing), ~mono(log(x)))) {
    expect_error(tinyAM:::.parse_q_formula(f, d), "existing index column")
  }
  expect_error(tinyAM:::.parse_q_formula(~mono(x, direction = "increasing"), d), "Invalid mono")
  expect_error(tinyAM:::.parse_q_formula(~mono(x, extra = 1), d), "Invalid mono")
  for (f in list(~x + mono(x), ~I(x^2) + mono(x), ~mono(x) + mono(x, by = survey))) {
    expect_error(tinyAM:::.parse_q_formula(f, d), "ordinary, duplicate")
  }
  for (f in list(~survey * mono(x), ~I(mono(x)), ~mono(x):survey, ~1 - mono(x))) {
    expect_error(tinyAM:::.parse_q_formula(f, d), "additive term")
  }
  expect_error(tinyAM:::.parse_q_formula(~mono(x), data.frame(x = 1)), "at least two")
  d$survey[2] <- "B"
  expect_error(tinyAM:::.parse_q_formula(~mono(x, by = survey), d), "at least two")
  expect_error(tinyAM:::.parse_q_formula(~mono(survey), d), "numeric or a factor")
  d$x[1] <- NA
  expect_error(tinyAM:::.parse_q_formula(~mono(x), d), "missing")
  expect_error(mono(1:4), "only supported")
})

test_that("multiple surveys map independent monotonic curves into the likelihood", {
  obs <- cod_obs
  a <- b <- obs$index
  a$survey <- "A"
  b$survey <- "B"
  a$q_block <- pmin(a$age, 5)
  b$q_block <- pmin(b$age, 4)
  obs$index <- rbind(a, b)
  obs$index$survey <- factor(obs$index$survey)
  dat <- make_dat(obs = obs, years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~survey + mono(q_block, by = survey),
                          sd_form = ~1, fill_missing = TRUE))
  par <- make_par(dat)
  par$log_q[] <- c(-2, 1)
  par$log_dq[] <- log(c(.1, .2, .3, .4, .5))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
  index <- dat$obs$index
  expected <- ifelse(index$survey == "A", c(-2, -1.9, -1.7, -1.4)[index$q_block - 1],
                     c(-1, -.6, -.1)[index$q_block - 1])
  expect_equal(unname(obj$report()$log_q_obs), expected)
})

test_that("ordinary q matrices, parameters, predictions and simulation remain unchanged", {
  ordinary <- data.frame(mono = 1:4, x = 4:1)
  for (f in list(~mono, x ~ mono)) {
    expect_identical(tinyAM:::.parse_q_formula(f, ordinary)$q_modmat,
                     stats::model.matrix(f, ordinary))
  }
  for (f in list(~q_block, ~1)) {
    dat <- make_test_dat(years = 2000:2005, ages = 2:6,
                         index_settings = list(q_form = f, sd_form = ~1, fill_missing = TRUE))
    old <- dat
    old$q_modmat <- stats::model.matrix(f, data = dat$obs$index)
    expect_identical(dat$q_modmat, old$q_modmat)
    par <- make_par(dat)
    expect_false("log_dq" %in% names(par))
    expect_equal(length(par$log_q), ncol(old$q_modmat))
    par$log_q[] <- seq_along(par$log_q) / 10
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
    expect_equal(obj$report()$log_q_obs, drop(old$q_modmat %*% par$log_q), ignore_attr = TRUE)
    expect_identical(nll_fun(par, dat), nll_fun(par, old))
    set.seed(412); actual <- nll_fun(par, dat, simulate = TRUE)
    set.seed(412); expected <- nll_fun(par, old, simulate = TRUE)
    expect_identical(actual, expected)
  }
})

test_that("monotonic q integrates with RTMB, simulation and observation predictions", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
                       F_settings = list(process = "iid", mu_form = ~1),
                       index_settings = list(q_form = ~mono(q_block), sd_form = ~1, fill_missing = TRUE))
  par <- make_par(dat)
  expect_equal(unname(par$log_dq), rep(log(.05), ncol(dat$q_mono_modmat)))
  expect_identical(names(par$log_dq), dat$q_mono_steps$coef)
  par$log_q[] <- -2
  par$log_dq[] <- log(seq_along(par$log_dq) / 10)
  par$log_mu_f[] <- log(.3)
  par$log_sd_f <- log(.1)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
  expected_q <- mono_values(dat, par$log_q, par$log_dq)
  expect_equal(obj$report()$log_q_obs, expected_q, ignore_attr = TRUE)
  calls <- list()
  testthat::local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    calls[[length(calls) + 1L]] <<- mean
    rep_len(mean, n)
  }, .package = "stats")
  set.seed(624)
  sim <- nll_fun(par, dat, simulate = TRUE)
  par[intersect(names(par), names(sim))] <- sim[intersect(names(par), names(sim))]
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expect_equal(report$log_q_obs, expected_q, ignore_attr = TRUE)
  expect_equal(tail(calls, 2)[[1]], report$log_pred[dat$is_observed])
  ii <- dat$obs_map$type == "index"
  rows <- cbind(match(dat$obs$index$year, dat$years), match(dat$obs$index$age, dat$ages))
  expect_equal(unname(report$log_pred[ii]), unname(expected_q + log(report$N)[rows] -
                 report$Z[rows] * dat$obs$index$samp_time))
})

test_that("constant and monotonic q fit and retain useful tidy outputs", {
  for (f in list(~1, ~mono(q_block))) {
    fit <- update(default_fit, index_settings = list(q_form = f, sd_form = ~1, fill_missing = TRUE),
                  start_par = as.list(default_fit$sdrep, "Estimate"), silent = TRUE)
    expect_true(is.finite(fit$opt$objective))
    expect_true(all(is.finite(fit$obj$gr(fit$opt$par))))
    expect_equal(fit$obs_pred$index$q, unname(exp(fit$rep$log_q_obs)))
    if (length(as.list(fit$sdrep, "Estimate")$log_dq)) {
      steps <- subset(tidy_par(fit)$fixed, par == "dq")
      expect_equal(steps$coef, fit$dat$q_mono_steps$coef)
      expect_equal(steps$est, unname(exp(as.list(fit$sdrep, "Estimate")$log_dq)))
      by_block <- tapply(fit$obs_pred$index$q, fit$obs_pred$index$q_block, unique)
      expect_true(all(lengths(by_block) == 1L))
      expect_true(all(diff(unlist(by_block)) >= 0))
      expect_false("log_dq" %in% names(fit$random_par))
    } else {
      expect_length(unique(fit$obs_pred$index$q), 1L)
    }
  }
})
