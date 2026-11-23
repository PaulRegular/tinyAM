test_that("Poisson-link delta preserves the mean and has stable derivatives", {
  y50 <- 20
  sd <- .5
  mu <- y50 * c(.001, .1, 1, 2, 10)
  p <- .poisson_link_positive(log(mu), log(y50), sd)
  probability <- -expm1(-log(2) * mu / y50)
  expect_equal(exp(p$log_p_positive), probability, tolerance = 1e-12)
  expect_equal(exp(-exp(p$log_intensity))[3:4], c(.5, .25), tolerance = 1e-14)
  positive_mean <- exp(p$log_median + sd^2 / 2)
  expect_equal(positive_mean * probability, mu, tolerance = 1e-11)
  expect_equal(probability * (positive_mean^2 * exp(sd^2)) - mu^2,
    mu^2 * (exp(sd^2) / probability - 1), tolerance = 1e-11)
  scaled <- .poisson_link_positive(log(mu * 1000), log(y50 * 1000), sd)
  expect_equal(scaled$log_median - log(1000), p$log_median, tolerance = 1e-12)
  expect_equal(scaled$log_p_positive, p$log_p_positive, tolerance = 1e-12)

  nll <- function(p) {
    result <- .poisson_link_positive(p$log_mu, log(y50), sd)
    -sum(result$log_p_positive + RTMB::dnorm(2, result$log_median, sd, log = TRUE))
  }
  tape <- RTMB::MakeADFun(nll, list(log_mu = 0), silent = TRUE)
  for (x in c(-1000, -50, 0, 50, 1000)) {
    expect_true(all(is.finite(c(tape$fn(x), tape$gr(x), tape$he(x)))))
  }
  small <- .poisson_link_positive(-1000, log(y50), sd)
  expect_equal(small$log_median, log(y50 / log(2)) - sd^2 / 2, tolerance = 1e-13)
})

test_that("only supported delta zeros join the likelihood", {
  d <- cod_obs
  d$catch$obs[d$catch$age == 2] <- 0
  d$catch$obs[d$catch$age == 3 & d$catch$year == 1991] <- NA_real_
  d$catch$obs[d$catch$age == 4 & d$catch$year == 1991] <- 0
  args <- list(data = d, years = 1991:1994,
    F_settings = list(age_breaks = 3:21),
    catch_settings = list(family = delta_lognormal(type = "poisson-link")))
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
  gaussian <- dat$delta$lognormal
  positive <- dat$delta$positive
  zero <- dat$delta$zero
  expected <- -sum(stats::dnorm(dat$log_obs[gaussian],
    report$log_pred[gaussian], report$sd_obs[gaussian], log = TRUE)) -
    sum(report$log_p_zero[zero]) -
    sum(report$log_p_positive[positive] + stats::dnorm(dat$log_obs[positive],
      report$log_pred_positive[positive], report$sd_obs[positive], log = TRUE))
  expect_equal(make_tape(obs)$fn(), expected, tolerance = 1e-10)
  sim <- nll_fun(make_par(dat), dat, simulate = TRUE)
  expect_true(all(sim$log_obs[structural] == -Inf))
  expect_true(is.na(sim$log_obs[missing]))
  args$catch_settings$cor_struct <- "ar1"
  expect_error(do.call(prepare_tam, args), "requires independent")
  args$catch_settings$cor_struct <- NULL
  args$index_settings <- list(cor_struct = "ar1")
  mixed <- do.call(prepare_tam, args)
  expect_identical(mixed$obs_correlation$groups, "index:Fall 2J3KL")
  expect_true(all(is.finite(make_tape(mixed)$gr())))
})

test_that("delta zeros do not support observation SD coefficients", {
  d <- cod_obs
  d$index$obs[d$index$age == 2] <- 0
  # Use a single numerical indicator so this checks the unsupported direction
  # independently of other ages that also have limited positive observations.
  settings <- list(family = delta_lognormal(type = "poisson-link"), sd_form = ~ 0 + I(as.numeric(age == 2)))
  expect_error(prepare_tam(d, years = 1991:1998, index_settings = settings),
    "index_settings.*rank-deficient.*eligible observations")
  settings$sd_form <- ~ 0
  settings$sd_supplied <- ~ I(0.3)
  dat <- prepare_tam(d, years = 1991:1998, index_settings = settings)
  expect_equal(ncol(dat$sd_index_modmat), 0L)
})

test_that("zero generalized gamma shape matches delta lognormal for both types", {
  make_tape <- function(family) {
    dat <- prepare_tam(cod_obs, years = 1991:1994,
      catch_settings = list(family = family, sd_form = ~ 1))
    RTMB::MakeADFun(function(p) nll_fun(p, dat), make_par(dat), silent = TRUE)
  }
  for (type in c("logit-link", "poisson-link")) {
    normal <- make_tape(delta_lognormal(type = type))
    zero_shape <- make_tape(delta_gengamma(type = type, Q = 0))
    expect_identical(zero_shape$par, normal$par)
    expect_equal(zero_shape$fn(), normal$fn(), tolerance = 1e-12)
    expect_equal(zero_shape$gr(), normal$gr(), tolerance = 1e-12)
  }
})

test_that("delta series have separate estimated or fixed scales with clear support", {
  d <- cod_obs
  d$index$survey <- "A"
  b <- d$index
  b$survey <- "B"
  b$obs <- ifelse(is.na(b$obs), NA, pmax(1, b$obs))
  d$index <- rbind(d$index, b)
  family <- list(A = delta_lognormal(), B = delta_lognormal("poisson-link", y50 = 100))
  args <- list(data = d, years = 1991:1998,
    index_settings = list(family = family, q_form = ~ survey))
  dat <- do.call(prepare_tam, args)
  par <- make_par(dat)
  expect_identical(names(par$log_y50_index), "A")
  expect_null(par$log_y50_catch)
  expect_identical(unique(dat$obs_map$delta_type[dat$obs_map$survey == "A" &
    dat$obs_map$type == "index"]), "logit-link")
  expect_equal(dat$delta$series[[2L]]$log_start, log(100))
  expect_error(do.call(prepare_tam, within(args, index_settings$family <- delta_lognormal(type = "poisson-link"))),
    "no eligible sampled zeros")
  d$index$obs[d$index$survey == "B"] <- 0
  args$data <- d
  expect_error(do.call(prepare_tam, args), "needs at least one positive")
  expect_error(delta_lognormal(type = "poisson-link", y50 = 0), "positive finite")
})

test_that("logit-link delta respects its half-detection scale and mean", {
  y50 <- 20
  sd <- 0.5
  mean <- y50 * c(.001, .1, 1, 2, 10)
  p <- .delta_positive(log(mean), log(y50), sd, "logit-link")
  expect_equal(exp(p$log_p_zero), y50 / (mean + y50), tolerance = 1e-13)
  expect_equal(exp(p$log_p_positive + p$log_mean), mean, tolerance = 1e-13)
  expect_equal(exp(p$log_mean), mean + y50, tolerance = 1e-13)
  scaled <- .delta_positive(log(mean * 1000), log(y50 * 1000), sd, "logit-link")
  expect_equal(scaled$log_p_positive, p$log_p_positive, tolerance = 1e-13)
  expect_equal(scaled$log_mean - log(1000), p$log_mean, tolerance = 1e-13)

  nll <- function(par) {
    pred <- .delta_positive(par$log_mu, log(y50), sd, "logit-link")
    -pred$log_p_zero - pred$log_p_positive -
      RTMB::dnorm(2, pred$log_median, sd, log = TRUE)
  }
  tape <- RTMB::MakeADFun(nll, list(log_mu = 0), silent = TRUE)
  for (x in c(-1000, -50, 0, 50, 1000)) {
    expect_true(all(is.finite(c(tape$fn(x), tape$gr(x), tape$he(x)))))
  }
})

test_that("renamed Poisson-link and ordinary families preserve numeric results", {
  # Captured on 6d71bf37 from the same fixed input and starting parameters.
  reference <- dget(test_path("fixtures", "observation-family-reference.R"))
  cases <- list(delta_log = delta_lognormal("poisson-link"),
    delta_gamma = delta_gamma("poisson-link"),
    delta_gg = delta_gengamma("poisson-link", Q = .23), normal = lognormal(),
    gamma = stats::Gamma(link = "log"), gg = gengamma(Q = .23))
  d <- cod_obs
  d$catch$obs[d$catch$age == 2] <- 0
  d$catch$obs[d$catch$age == 4 & d$catch$year == 1991] <- 0
  d$catch$obs[d$catch$age == 3 & d$catch$year == 1991] <- NA_real_

  for (name in names(cases)) {
    dat <- prepare_tam(d, years = 1991:1994, ages = 2:8,
      F_settings = list(age_breaks = 3:9), catch_settings = list(family = cases[[name]]))
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), make_par(dat), silent = TRUE)
    grad <- as.numeric(obj$gr())
    report <- obj$report()
    selected <- which(dat$obs_map$type == "catch" & dat$obs_map$year == 1991)
    ref <- reference[[name]]
    expect_equal(obj$fn(), ref$objective, tolerance = 1e-9, info = name)
    expect_equal(c(sum(grad), sqrt(sum(grad^2)), sum(grad * seq_along(grad))),
      ref$gradient, tolerance = 1e-9, info = name)
    for (field in setdiff(names(ref), c("objective", "gradient"))) {
      expect_equal(report[[field]][selected], ref[[field]], tolerance = 1e-9,
        info = paste(name, field))
    }
  }
})
