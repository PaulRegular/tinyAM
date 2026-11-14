test_that("PLD preserves the mean and has stable derivatives across detection levels", {
  y50 <- 20
  sd <- .5
  mu <- y50 * c(.001, .1, 1, 2, 10)
  p <- .pld_positive(log(mu), log(y50), sd)
  probability <- -expm1(-log(2) * mu / y50)
  expect_equal(exp(p$log_p_positive), probability, tolerance = 1e-12)
  expect_equal(exp(-exp(p$log_intensity))[3:4], c(.5, .25), tolerance = 1e-14)
  positive_mean <- exp(p$log_median + sd^2 / 2)
  expect_equal(positive_mean * probability, mu, tolerance = 1e-11)
  expect_equal(probability * (positive_mean^2 * exp(sd^2)) - mu^2,
    mu^2 * (exp(sd^2) / probability - 1), tolerance = 1e-11)
  scaled <- .pld_positive(log(mu * 1000), log(y50 * 1000), sd)
  expect_equal(scaled$log_median - log(1000), p$log_median, tolerance = 1e-12)
  expect_equal(scaled$log_p_positive, p$log_p_positive, tolerance = 1e-12)

  nll <- function(p) {
    result <- .pld_positive(p$log_mu, log(y50), sd)
    -sum(result$log_p_positive + RTMB::dnorm(2, result$log_median, sd, log = TRUE))
  }
  tape <- RTMB::MakeADFun(nll, list(log_mu = 0), silent = TRUE)
  for (x in c(-1000, -50, 0, 50, 1000)) {
    expect_true(all(is.finite(c(tape$fn(x), tape$gr(x), tape$he(x)))))
  }
  small <- .pld_positive(-1000, log(y50), sd)
  expect_equal(small$log_median, log(y50 / log(2)) - sd^2 / 2, tolerance = 1e-13)
})

test_that("only supported PLD zeros join the likelihood", {
  d <- cod_obs
  d$catch$obs[d$catch$age == 2] <- 0
  d$catch$obs[d$catch$age == 3 & d$catch$year == 1991] <- NA_real_
  d$catch$obs[d$catch$age == 4 & d$catch$year == 1991] <- 0
  args <- list(data = d, years = 1991:1994,
    F_settings = list(age_breaks = 3:21),
    catch_settings = list(family = pld_lognormal()))
  dat <- do.call(prepare_tam, args)
  catch <- dat$obs_map$type == "catch"
  structural <- catch & dat$obs_map$age == 2
  ordinary <- catch & dat$obs_map$age == 4 & dat$obs_map$year == 1991
  missing <- catch & dat$obs_map$age == 3 & dat$obs_map$year == 1991
  expect_true(all(dat$obs_map$is_structural[structural]))
  expect_false(any(dat$is_observed[structural | missing]))
  expect_true(dat$is_observed[ordinary])
  expect_identical(dat$log_obs[ordinary], -Inf)

  make_tape <- function(dat) RTMB::MakeADFun(function(p) nll_fun(p, dat),
    make_par(dat), silent = TRUE)
  obj <- make_tape(dat)
  report <- obj$report()
  expect_equal(report$log_p_zero[structural], rep(0, sum(structural)))
  expect_true(all(is.finite(obj$gr())))
  obs <- dat
  obs$.likelihood_component <- "observation"
  gaussian <- dat$pld$lognormal
  positive <- dat$pld$positive
  zero <- dat$pld$zero
  expected <- -sum(stats::dnorm(dat$log_obs[gaussian],
    report$log_pred[gaussian], report$sd_obs[gaussian], log = TRUE)) -
    sum(report$log_p_zero[zero]) -
    sum(report$log_p_positive[positive] + stats::dnorm(dat$log_obs[positive],
      report$log_pred_positive[positive], report$sd_obs[positive], log = TRUE))
  expect_equal(make_tape(obs)$fn(), expected, tolerance = 1e-10)
  sim <- nll_fun(make_par(dat), dat, simulate = TRUE)
  expect_true(all(sim$log_obs[structural] == -Inf))
  expect_true(is.na(sim$log_obs[missing]))

  args$ages <- 2:14
  expect_error(do.call(prepare_tam, args), "Supply pooled.*plus age")
  args$ages <- NULL
  args$catch_settings$cor_struct <- "ar1"
  expect_error(do.call(prepare_tam, args), "requires independent")
  args$catch_settings$cor_struct <- NULL
  args$index_settings <- list(cor_struct = "ar1")
  mixed <- do.call(prepare_tam, args)
  expect_identical(mixed$obs_correlation$groups, "index:Fall 2J3KL")
  expect_true(all(is.finite(make_tape(mixed)$gr())))
})

test_that("PLD series have separate estimated or fixed scales with clear support", {
  d <- cod_obs
  d$index$survey <- "A"
  b <- d$index
  b$survey <- "B"
  b$obs <- ifelse(is.na(b$obs), NA, pmax(1, b$obs))
  d$index <- rbind(d$index, b)
  family <- list(A = pld_lognormal(), B = pld_lognormal(100))
  args <- list(data = d, years = 1991:1998,
    index_settings = list(family = family, q_form = ~ survey))
  dat <- do.call(prepare_tam, args)
  par <- make_par(dat)
  expect_identical(names(par$log_y50_index), "A")
  expect_null(par$log_y50_catch)
  expect_equal(dat$pld$series[[2L]]$log_start, log(100))
  expect_error(do.call(prepare_tam, within(args, index_settings$family <- pld_lognormal())),
    "no eligible sampled zeros")
  d$index$obs[d$index$survey == "B"] <- 0
  args$data <- d
  expect_error(do.call(prepare_tam, args), "needs at least one positive")
  expect_error(pld_lognormal(0), "positive finite")
  expect_error(lognormal("identity"), "link")
})
