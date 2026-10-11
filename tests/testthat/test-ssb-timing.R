timed_obs <- function(time = .25) {
  data <- cod_obs
  data$maturity$spawn_time <- time
  data
}

test_that("spawning-time survival uses the supplied fraction of F and M", {
  for (time in c(0, .25, 1)) {
    dat <- make_test_dat(years = 1983:1990, ages = 2:8, data = timed_obs(time))
    par <- make_par(dat)
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    # Shortened plus-group biology is constructed inside the population path.
    population <- tinyAM:::.population_states(par, dat, report$Z)
    expected <- report$N * population$W * population$P * exp(-time * report$Z)
    expect_equal(report$ssb_mat, expected)
    expect_equal(unname(report$ssb), unname(rowSums(expected)))
    expect_equal(report$biomass, rowSums(report$N * population$W), ignore_attr = TRUE)
  }
  for (time in list(-.1, 1.1, NA_real_, Inf, "fall")) {
    expect_error(make_test_dat(data = timed_obs(time)), "spawn_time")
  }
})

test_that("simulated BH/Ricker parents use the same spawning-time population", {
  for (curve in c("bh", "ricker")) {
    dat <- make_test_dat(years = 1983:1995, ages = 2:8,
      N_settings = list(process = "iid", rec_form = as.formula(paste0("~ ", curve, "(ssb) + iid(year, sd = .2)"))),
      data = timed_obs(.25),
      proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1))
    par <- make_par(dat)
    set.seed(731)
    simulated <- nll_fun(par, dat, simulate = TRUE)
    shared <- intersect(names(par), names(simulated))
    par[shared] <- simulated[shared]
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    population <- tinyAM:::.population_states(par, dat, report$Z)
    i <- dat$rec$eligible
    parent <- log(report$ssb[i - dat$rec$curve$lag])
    expect_equal(unname(report$rec_log_parent), unname(parent))
    expect_equal(population$log_mu_R[i], tinyAM:::.rec_log_curve(parent, par, curve), ignore_attr = TRUE)
    expect_true(all(is.finite(report$ssb)))
  }
})

test_that("curve initialization anchors the configured spawning-time SSB", {
  for (curve in c("bh", "ricker")) for (time in c(0, .25, 1)) {
    dat <- make_test_dat(years = 1983:1995, ages = 2:8,
      N_settings = list(process = "off", rec_form = as.formula(paste0("~ ", curve, "(ssb) + iid(year)"))),
      data = timed_obs(time))
    par <- make_par(dat)
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    parent <- dat$rec$eligible[1L] - dat$rec$curve$lag
    S <- report$ssb[parent]
    expect_equal(unname(exp(par$log_sr_beta)), unname(1 / S))
    expect_equal(unname(exp(tinyAM:::.rec_log_curve(log(S), par, curve))), exp(par$log_r0))
  }
})

test_that("spawning timing survives updates and retrospective/hindcast folds", {
  fit <- suppressWarnings(update(default_fit, data = timed_obs(.25), silent = TRUE))
  expect_true(all(fit$pop$ssb$spawn_time == .25))
  expect_true(all(is.finite(fit$pop$ssb$se)))
  refit <- suppressWarnings(update(fit, years = head(YEARS, -1), silent = TRUE,
    start_par = as.list(fit$sdrep, "Estimate"), proj_settings = list(n_proj = 1, n_mean = 1, F_mult = 1)))
  expect_identical(unique(as.vector(refit$dat$spawn_time)), .25)
  expect_true(all(refit$pop$ssb$spawn_time == .25))
  retro <- fit_retro(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  hind <- fit_hindcast(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  expect_true(all(vapply(retro$fits, function(x) all(x$dat$spawn_time == .25), logical(1))))
  expect_true(all(vapply(hind$fits, function(x) all(x$dat$spawn_time == .25), logical(1))))
})

test_that("default SSB timing preserves likelihood and population states", {
  dat <- make_test_dat()
  zero <- make_test_dat(data = timed_obs(0))
  expect_identical(make_par(dat), make_par(zero))
  expect_identical(nll_fun(make_par(dat), dat), nll_fun(make_par(zero), zero))
})

test_that("age/year timing is aligned and contributes exact survival", {
  data <- cod_obs
  before <- data
  data$maturity$spawn_time <- .1 + .01 * (data$maturity$age - 2) + .001 * (data$maturity$year - 1954)
  input <- data
  dat <- make_test_dat(data = data, years = 1983:1990, ages = 2:max(data$maturity$age))
  expect_identical(data, input)
  par <- make_par(dat)
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expected <- report$N * dat$W * dat$P * exp(-dat$spawn_time * report$Z)
  expect_equal(report$ssb_mat, expected)
  expect_equal(report$ssb, rowSums(expected), ignore_attr = TRUE)
  expect_identical(data$catch, before$catch)
  expect_match(tinyAM:::.ssb_timing_label(dat), "Age-specific")
})

test_that("hidden ages retain their own timing and project terminal values", {
  data <- timed_obs()
  data$maturity$spawn_time <- .01 * data$maturity$age + .001 * (data$maturity$year - 1954)
  dat <- make_test_dat(data = data, years = 1983:1990, ages = 2:8,
    proj_settings = list(n_proj = 2))
  par <- make_par(dat)
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expected <- rowSums(report$N_plus * dat$W_plus_input * dat$P_plus_input *
    exp(-dat$spawn_time_plus * report$Z[, "8"]))
  expect_equal(unname(report$ssb_mat[, "8"]), unname(expected))
  averaged <- rowSums(report$N_plus * dat$W_plus_input * dat$P_plus_input) *
    exp(-rowMeans(dat$spawn_time_plus) * report$Z[, "8"])
  expect_gt(max(abs(expected - averaged)), 1e-6)
  for (matrix in list(dat$spawn_time, dat$spawn_time_plus)) {
    expect_equal(unname(matrix[9, ]), unname(matrix[8, ]))
    expect_equal(matrix[9, ], matrix[10, ])
  }
  expect_equal(sum(report$ssb_mat[1, ]), report$ssb[1], ignore_attr = TRUE)
})

test_that("timing is validated even at immature and hidden biological ages", {
  for (value in list(NA_real_, NaN, Inf, -.1, 1.1, "fall", factor("fall"))) {
    data <- timed_obs()
    data$maturity$spawn_time <- value
    expect_error(check_obs(data), "spawn_time")
    expect_error(make_test_dat(data = data, ages = 2:8), "spawn_time")
  }
  data <- timed_obs()
  data$maturity$spawn_time[data$maturity$age == max(data$maturity$age)] <- NA
  expect_error(make_test_dat(data = data, ages = 2:8), "spawn_time")
})

test_that("zero-lag parent SSB uses age-specific timing without circularity", {
  data <- timed_obs()
  data$maturity$obs[data$maturity$age == 2] <- 0
  data$maturity$spawn_time <- data$maturity$age / 20
  dat <- make_test_dat(data = data, years = 1983:1995, ages = 2:8,
    N_settings = list(process = "off", rec_form = ~ bh(ssb, lag = 0) + iid(year)))
  par <- make_par(dat)
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expect_equal(report$rec_log_parent, log(report$ssb[dat$rec$eligible]), ignore_attr = TRUE)
  data$maturity$obs[data$maturity$age == 2] <- .1
  expect_error(make_test_dat(data = data, N_settings = list(rec_form = ~ bh(ssb, lag = 0) + iid(year))), "circular")
})

test_that("variable timing is reported honestly and retains derived uncertainty", {
  data <- timed_obs()
  data$maturity$spawn_time <- data$maturity$age / 30
  fit <- suppressWarnings(update(default_fit, data = data, silent = TRUE,
    start_par = as.list(default_fit$sdrep, "Estimate")))
  expect_true(all(is.na(fit$pop$ssb$spawn_time)))
  expect_true(all(is.finite(fit$pop$ssb$se)))
  expect_true(all(grepl("Age-specific", fit$pop$ssb$timing)))
  timing <- tinyAM:::.ssb_timing_table(fit$dat)
  expect_equal(timing$spawn_time, timing$age / 30)
  for (curve in c("bh", "ricker")) {
    dat <- make_test_dat(data = data, years = 1983:1995, ages = 2:8,
      N_settings = list(rec_form = as.formula(paste0("~ ", curve, "(ssb) + iid(year, sd = .2)"))))
    par <- make_par(dat)
    set.seed(897)
    simulated <- nll_fun(par, dat, simulate = TRUE)
    par[intersect(names(par), names(simulated))] <- simulated[intersect(names(par), names(simulated))]
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    expect_equal(report$rec_log_parent, log(report$ssb[dat$rec$eligible - dat$rec$curve$lag]), ignore_attr = TRUE)
    anchor <- make_par(dat)
    anchored <- RTMB::MakeADFun(function(p) nll_fun(p, dat), anchor, silent = TRUE)$report()
    parent <- dat$rec$eligible[1] - dat$rec$curve$lag
    expect_equal(exp(anchor$log_sr_beta), 1 / anchored$ssb[parent], ignore_attr = TRUE)
  }
})
