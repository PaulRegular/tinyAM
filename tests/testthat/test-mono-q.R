mono_values <- function(design, beta, dq) {
  drop(design$q_modmat %*% beta + design$q_mono_modmat %*% dq)
}

test_that("monotonic q uses cumulative non-negative log-q steps", {
  d <- data.frame(x = 1:4)
  design <- tinyAM:::.parse_q_formula(~mono(x), d)
  expect_equal(unname(mono_values(design, -2, c(.1, .2, .3))),
               c(-2, -1.9, -1.7, -1.4))
  expect_equal(unname(mono_values(design, -2, c(0, .2, 0))),
               c(-2, -2, -1.8, -1.8))
  set.seed(823)
  for (i in 1:10) {
    log_q <- mono_values(design, -2, runif(3, 0, 1))
    expect_true(all(diff(log_q) > 0))
    expect_true(all(diff(exp(log_q)) > 0))
  }
})

test_that("a fitted plateau has finite direct-scale increments and SEs", {
  design <- tinyAM:::.parse_q_formula(~mono(x), data.frame(x = 1:4))
  par <- list(log_q = c("(Intercept)" = -2),
              dq = setNames(rep(.05, 3), colnames(design$q_mono_modmat)))
  # Decreasing observations put the constrained optimum exactly on all bounds.
  y <- c(-1.7, -1.9, -2.1, -2.3)
  obj <- RTMB::MakeADFun(function(p) {
    pred <- design$q_modmat %*% p$log_q + design$q_mono_modmat %*% p$dq
    -sum(RTMB::dnorm(y, pred, .1, log = TRUE))
  }, par, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = c(-Inf, 0, 0, 0))
  obj$fn(opt$par)
  sdrep <- RTMB::sdreport(obj)
  estimates <- as.list(sdrep, "Estimate")
  expect_equal(unname(estimates$dq), rep(0, 3))
  expect_equal(unname(exp(mono_values(design, estimates$log_q, estimates$dq))),
               rep(exp(-2), 4), tolerance = 1e-6)
  fit <- structure(list(obj = obj, dat = list(), sdrep = sdrep), class = "tam_fit")
  steps <- subset(tidy_par(fit)$fixed, par == "dq")
  expect_equal(steps$est, rep(0, 3))
  expect_equal(steps$se, unname(as.list(sdrep, "Std. Error")$dq))
  expect_true(all(is.finite(steps$se) & steps$se > 0 & steps$se < 1))
  expect_true(check_convergence(fit))
})

test_that("convergence diagnostics respect only valid active-bound gradients", {
  fit <- list(sdrep = list(par.fixed = c(dq = 0), gradient.fixed = 2, pdHess = TRUE))
  expect_true(check_convergence(fit))
  expect_identical(fit$sdrep$gradient.fixed, 2)
  fit$sdrep$gradient.fixed <- Inf
  expect_warning(expect_false(check_convergence(fit)), "may not have converged")
  fit$sdrep$gradient.fixed <- -2
  expect_warning(expect_false(check_convergence(fit)), "may not have converged")
  fit$sdrep$par.fixed[] <- .1
  fit$sdrep$gradient.fixed <- 2
  expect_warning(expect_false(check_convergence(fit)), "may not have converged")
})

test_that("numeric order, declared factor order and pooled plateaus are retained", {
  for (x in list(c(10, 2, 5, 5), factor(c("old", "young", "mid", "mid"),
                  levels = c("young", "mid", "old")),
                ordered(c("old", "young", "mid", "mid"),
                  levels = c("young", "mid", "old")))) {
    design <- tinyAM:::.parse_q_formula(~ mono(x), data.frame(x = x))
    expect_equal(unname(mono_values(design, -2, c(.1, .2))), c(-1.7, -2, -1.9, -1.9))
  }
})

test_that("surveys have separate baselines and steps with different level counts", {
  d <- data.frame(survey = factor(c("B", "A", "A", "B", "A", "B", "A")),
                  x = factor(c(4, 2, 4, 1, 1, 3, 3), levels = 1:4))
  design <- tinyAM:::.parse_q_formula(~ survey + mono(x, by = survey), d)
  expect_equal(ncol(design$q_modmat), 2L)
  expect_equal(ncol(design$q_mono_modmat), 5L)
  expect_equal(design$q_mono_steps$by_level, c(rep("A", 3), rep("B", 2)))
  dq <- c(.1, .2, .3, .4, .5)
  q <- mono_values(design, c(-2, 1), dq)
  expect_equal(unname(q), c(-.1, -1.9, -1.4, -1, -2, -.6, -1.7))
  dq[1] <- .8
  changed <- mono_values(design, c(-2, 1), dq)
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
  par$dq[] <- c(.1, .2, .3, .4, .5)
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
    expect_false("dq" %in% names(par))
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
  expect_equal(unname(par$dq), rep(.05, ncol(dat$q_mono_modmat)))
  expect_identical(names(par$dq), dat$q_mono_steps$coef)
  par$log_q[] <- -2
  par$dq[] <- seq_along(par$dq) / 10
  par$log_mu_f[] <- log(.3)
  par$log_sd_f <- log(.1)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
  expected_q <- mono_values(dat, par$log_q, par$dq)
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

test_that("simulation projects negative uncertainty draws onto feasible dq", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~mono(q_block), sd_form = ~1, fill_missing = TRUE))
  par <- make_par(dat)
  par$dq[] <- rep_len(c(-.2, .3, 0), length(par$dq))
  feasible <- par
  feasible$dq[] <- pmax(par$dq, 0)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), feasible, silent = TRUE)
  fit <- list(dat = dat, obj = obj)
  for (redraw in c(FALSE, TRUE)) {
    set.seed(812)
    actual <- tinyAM:::.sim_obs(fit, function(fit) par, redraw_random = redraw)
    set.seed(812)
    expected <- tinyAM:::.sim_obs(fit, function(fit) feasible, redraw_random = redraw)
    expect_identical(actual, expected)
  }
})

test_that("negative monotonic warm starts fail clearly", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~mono(q_block), sd_form = ~1, fill_missing = TRUE))
  par <- make_par(dat)
  par$dq[1] <- -.1
  expect_error(fit_tam(obs = cod_obs, years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~mono(q_block), sd_form = ~1, fill_missing = TRUE),
    start_par = par, silent = TRUE), "finite and non-negative")
})

test_that("constant and monotonic q fit and retain useful tidy outputs", {
  for (f in list(~1, ~mono(q_block))) {
    fit <- update(default_fit, index_settings = list(q_form = f, sd_form = ~1, fill_missing = TRUE),
                  start_par = as.list(default_fit$sdrep, "Estimate"), silent = TRUE)
    expect_true(is.finite(fit$opt$objective))
    expect_true(all(is.finite(fit$obj$gr(fit$opt$par))))
    expect_equal(fit$obs_pred$index$q, unname(exp(fit$rep$log_q_obs)))
    if (length(as.list(fit$sdrep, "Estimate")$dq)) {
      steps <- subset(tidy_par(fit)$fixed, par == "dq")
      expect_equal(steps$coef, fit$dat$q_mono_steps$coef)
      expect_equal(steps$est, unname(as.list(fit$sdrep, "Estimate")$dq))
      expect_equal(steps$se, unname(as.list(fit$sdrep, "Std. Error")$dq))
      expect_true(all(steps$est >= 0))
      expect_true(any(steps$est == 0))
      expect_true(all(is.finite(steps$se)))
      by_block <- tapply(fit$obs_pred$index$q, fit$obs_pred$index$q_block, unique)
      expect_true(all(lengths(by_block) == 1L))
      expect_true(all(diff(unlist(by_block)) >= 0))
      expect_false("dq" %in% names(fit$random_par))
    } else {
      expect_length(unique(fit$obs_pred$index$q), 1L)
    }
  }
})
