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
