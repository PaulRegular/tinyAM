test_that("spawning-time survival uses the supplied fraction of F and M", {
  for (time in c(0, .25, 1)) {
    dat <- make_test_dat(years = 1983:1990, ages = 2:8, ssb_settings = list(spawn_time = time))
    par <- make_par(dat)
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    # Shortened plus-group biology is constructed inside the population path.
    population <- tinyAM:::.population_states(par, dat, report$Z)
    expected <- report$N * population$W * population$P * exp(-time * report$Z)
    expect_equal(report$ssb_mat, expected)
    expect_equal(unname(report$ssb), unname(rowSums(expected)))
    expect_equal(report$biomass, rowSums(report$N * population$W), ignore_attr = TRUE)
  }
  for (time in list(-.1, 1.1, NA_real_, Inf, c(.1, .2), "fall")) {
    expect_error(make_test_dat(ssb_settings = list(spawn_time = time)), "spawn_time")
  }
})

test_that("simulated BH/Ricker parents use the same spawning-time population", {
  for (curve in c("bh", "ricker")) {
    dat <- make_test_dat(years = 1983:1995, ages = 2:8,
      N_settings = list(process = "iid", rec_form = as.formula(paste0("~ ", curve, "(ssb) + iid(year, sd = .2)"))),
      ssb_settings = list(spawn_time = .25),
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
      ssb_settings = list(spawn_time = time))
    par <- make_par(dat)
    report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    parent <- dat$rec$eligible[1L] - dat$rec$curve$lag
    S <- report$ssb[parent]
    expect_equal(unname(exp(par$log_sr_beta)), unname(1 / S))
    expect_equal(unname(exp(tinyAM:::.rec_log_curve(log(S), par, curve))), exp(par$log_r0))
  }
})

test_that("spawning timing survives updates and retrospective/hindcast folds", {
  fit <- suppressWarnings(update(default_fit, ssb_settings = list(spawn_time = .25), silent = TRUE))
  expect_true(all(fit$pop$ssb$spawn_time == .25))
  expect_true(all(is.finite(fit$pop$ssb$se)))
  refit <- suppressWarnings(update(fit, years = head(YEARS, -1), silent = TRUE,
    start_par = as.list(fit$sdrep, "Estimate"), proj_settings = list(n_proj = 1, n_mean = 1, F_mult = 1)))
  expect_identical(refit$dat$ssb_settings, fit$dat$ssb_settings)
  expect_true(all(refit$pop$ssb$spawn_time == .25))
  retro <- fit_retro(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  hind <- fit_hindcast(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  expect_true(all(vapply(retro$fits, function(x) identical(x$dat$ssb_settings, fit$dat$ssb_settings), logical(1))))
  expect_true(all(vapply(hind$fits, function(x) identical(x$dat$ssb_settings, fit$dat$ssb_settings), logical(1))))
})

test_that("default SSB timing preserves likelihood and population states", {
  dat <- make_test_dat()
  zero <- make_test_dat(ssb_settings = list(spawn_time = 0))
  expect_identical(make_par(dat), make_par(zero))
  expect_identical(nll_fun(make_par(dat), dat), nll_fun(make_par(zero), zero))
})
