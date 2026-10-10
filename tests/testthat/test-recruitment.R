test_that("default recruitment retains the original random-walk density and layout", {
  d <- make_test_dat()
  p <- make_par(d)
  r <- c(0, seq_len(length(d$years) - 1L) / 10)
  expect_equal(d$rec$type, "rw")
  expect_equal(ncol(d$rec$matrix), 0L)
  expect_equal(names(p$log_r), as.character(d$years[-1L]))
  expect_false("rec_beta" %in% names(p))
  expect_equal(tinyAM:::.rec_nll(r, tinyAM:::.rec_mean(p, d), p, d),
               -sum(dnorm(diff(r), 0, exp(p$log_sd_r), log = TRUE)))
})

test_that("recruitment IID and AR1 use normalized residual densities", {
  for (type in c("iid", "ar1")) {
    form <- as.formula(paste0("~ ", type, "(year, sd = 0.3", if (type == "ar1") ", phi = 0.6" else "", ")"))
    d <- make_test_dat(N_settings = list(process = "off", rec_form = form))
    p <- make_par(d)
    p$rec_beta[] <- 2
    r <- seq_len(length(d$years)) / 10
    u <- r[-1L] - 2
    expected <- if (type == "iid") -sum(dnorm(u, 0, .3, log = TRUE)) else
      -dnorm(u[1], 0, .3 / sqrt(1 - .6^2), log = TRUE) - sum(dnorm(u[-1], .6 * head(u, -1), .3, log = TRUE))
    expect_equal(unname(tinyAM:::.rec_nll(r, tinyAM:::.rec_mean(p, d), p, d)), expected)
    expect_false("log_sd_r" %in% names(p))
    expect_true(is.finite(nll_fun(p, d)))
  }
})

test_that("recruitment covariates come only from youngest-age maturity rows", {
  obs <- tinyAM::cod_obs
  obs$maturity$temp <- ifelse(obs$maturity$age == min(AGES), sin(obs$maturity$year), NA)
  d <- make_test_dat(data = obs, N_settings = list(process = "off", rec_form = ~ temp + rw(year)))
  p <- make_par(d)
  p$rec_beta[] <- .4
  r <- seq_len(length(d$years)) / 10
  expect_equal(tinyAM:::.rec_nll(r, tinyAM:::.rec_mean(p, d), p, d),
               -sum(dnorm(diff(r) - .4 * diff(d$rec$data$temp), 0, 1, log = TRUE)))
  expect_error(make_test_dat(N_settings = list(rec_form = ~ absent + iid(year))), "covariates")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ iid(year) + ar1(year))), "exactly one")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ rw(age))), "year without")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ factor(year) + iid(year))), "saturated")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ year * rw(year))), "additive")
})

test_that("Gaussian recruitment simulation keeps the anchor and residual convention", {
  for (type in c("iid", "ar1", "rw")) {
    form <- as.formula(paste0("~ ", type, "(year, sd = 0.2", if (type == "ar1") ", phi = 0.5" else "", ")"))
    d <- make_test_dat(N_settings = list(process = "off", rec_form = form))
    p <- make_par(d)
    p$log_r0 <- 1
    if (!is.null(p$rec_beta)) p$rec_beta[] <- 2
    set.seed(912)
    draws <- replicate(500, tinyAM:::.simulate_rec(p, d))
    expected <- if (type == "rw") 1 else 2
    expect_lt(max(abs(rowMeans(draws) - expected)), if (type == "rw") .15 else .05)
    expect_equal(p$log_r0, 1)
  }
})

test_that("stock-recruit curves and derivatives match their definitions", {
  p <- list(log_sr_alpha = log(4), log_sr_beta = log(.02))
  S <- c(1, 30, 100)
  for (type in c("bh", "ricker")) {
    expected <- if (type == "bh") 4 * S / (1 + .02 * S) else 4 * S * exp(-.02 * S)
    expect_equal(exp(tinyAM:::.rec_log_curve(log(S), p, type)), expected)
    obj <- RTMB::MakeADFun(function(p) sum(tinyAM:::.rec_log_curve(log(S), p, type)), p, silent = TRUE)
    expected_grad <- c(3, if (type == "bh") -sum(.02 * S / (1 + .02 * S)) else -sum(.02 * S))
    expect_equal(unname(obj$gr(obj$par)), matrix(expected_grad, 1), tolerance = 1e-10)
  }
})

test_that("stock-recruit boundaries and lags align by year", {
  for (lag in c(1L, 2L, 4L)) {
    d <- make_test_dat(years = 1983:1995, ages = 2:8,
      N_settings = list(process = "off", rec_form = as.formula(substitute(~ bh(ssb, lag = L) + iid(year), list(L = lag)))))
    p <- make_par(d)
    expect_equal(names(p$log_r), as.character(d$years[seq.int(max(2L, lag + 1L), 13L)]))
    expect_equal(length(p$log_r_init), max(0L, lag - 1L))
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    expect_true(is.finite(obj$fn()))
    shortened <- make_test_dat(years = 1983:1990, ages = 2:8,
      N_settings = d$N_settings)
    p$log_r[] <- seq_along(p$log_r)
    start <- tinyAM:::.merge_start_par(make_par(shortened), p)
    expect_equal(start$log_r, p$log_r[names(start$log_r)])
  }
})

test_that("simulated recruitment uses parent SSB from the same population", {
  for (type in c("bh", "ricker")) for (n_process in c("off", "iid")) {
    form <- as.formula(paste0("~ ", type, "(ssb) + ar1(year, sd = 0.2, phi = 0.4)"))
    d <- make_test_dat(years = 1983:1995, ages = 2:8,
      N_settings = list(process = n_process, rec_form = form),
      proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1))
    p <- make_par(d)
    set.seed(713)
    simulated <- nll_fun(p, d, simulate = TRUE)
    p[names(simulated)[names(simulated) %in% names(p)]] <- simulated[names(simulated) %in% names(p)]
    report <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$report()
    states <- tinyAM:::.population_states(p, d, report$Z)
    i <- d$rec$eligible
    expected <- tinyAM:::.rec_log_curve(log(report$ssb[i - 2L]), p, type)
    expect_equal(unname(states$log_mu_R[i]), unname(expected))
    expect_equal(states$log_N, log(report$N))
    expect_equal(states$W * states$P * report$N, report$ssb_mat)
    expect_equal(simulated$log_r_init, p$log_r_init)
  }
})

test_that("zero-lag models reject circular SSB and unavailable parent data", {
  expect_error(make_test_dat(N_settings = list(rec_form = ~ bh(ssb, lag = 0) + iid(year))), "circular")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ bh(ssb) + rw(year))), "not RW")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ bh(ssb) + ricker(ssb) + iid(year))), "one stock-recruit")
  expect_error(make_test_dat(N_settings = list(rec_form = ~ bh(ssb, lag = -1) + iid(year))), "non-negative")
  expect_error(make_test_dat(years = 1983:1985, N_settings = list(rec_form = ~ bh(ssb, lag = 3) + iid(year))), "No historical")
  obs <- cod_obs
  obs$maturity$obs[obs$maturity$age == 2] <- 0
  d <- make_test_dat(data = obs, years = 1983:1995, ages = 2:8,
    N_settings = list(process = "off", rec_form = ~ bh(ssb, lag = 0) + iid(year)))
  p <- make_par(d)
  expect_true(is.finite(nll_fun(p, d)))
  set.seed(94)
  expect_true(all(is.finite(nll_fun(p, d, simulate = TRUE)$log_r)))
})
