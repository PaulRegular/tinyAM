test_that("Gamma specs use the standard log-link family object", {
  series <- c("RV", "Sentinel")
  direct <- .obs_families(stats::Gamma(link = "log"), series, "index_settings$family")
  expect_identical(vapply(direct, `[[`, character(1), "name"),
    c(RV = "gamma", Sentinel = "gamma"))
  mixed <- .obs_families(list(RV = pld_gengamma(y50 = 21, Q = -.2)),
    series, "index_settings$family")
  expect_identical(mixed$RV$name, "pld_gengamma")
  expect_identical(mixed$RV$y50, 21)
  expect_identical(mixed$RV$Q, -.2)
  expect_identical(mixed$Sentinel, lognormal())
  expect_error(.obs_families(stats::Gamma(), series, "index_settings$family"), "Gamma")
  expect_error(.obs_families(list(RV = stats::gaussian()), series,
    "index_settings$family"), "Gamma")
  expect_error(pld_gamma(y50 = 0), "positive finite")
  expect_error(gengamma(Q = Inf), "one finite number")
  expect_error(pld_gengamma(Q = c(.1, .2)), "one finite number")
  expect_error(gengamma(link = "identity"), "link")
})

test_that("PLD constructors expose the log link without changing positional arguments", {
  expect_identical(pld_lognormal(20), pld_lognormal(y50 = 20, link = "log"))
  expect_identical(pld_gamma(40), pld_gamma(y50 = 40, link = "log"))
  expect_identical(pld_gengamma(60, .2), pld_gengamma(y50 = 60, Q = .2, link = "log"))

  for (constructor in list(pld_lognormal, pld_gamma, pld_gengamma)) {
    expect_identical(constructor()$link, "log")
    for (link in list("cloglog", NA_character_, c("log", "log"), NULL)) {
      expect_error(constructor(link = link), "link = 'log'.*overall assessment prediction")
    }
  }
})

test_that("Gamma density uses a log observation measure and CV dispersion", {
  y <- c(.2, 2, 10, 120, 4)
  mean <- c(.8, 3, 5, 50, 2)
  cv <- c(.35, .25, .8, .65, 1.8)
  expected <- stats::dgamma(y, shape = cv^-2, scale = mean * cv^2, log = TRUE) + log(y)
  expect_equal(.gamma_log_density(log(y), log(mean), cv), expected, tolerance = 1e-13)
  expect_equal(.gamma_log_density(log(y) - 300, log(mean) - 300, cv),
    expected, tolerance = 1e-12)
  expect_equal(.positive_log_cdf("gamma", log(y), log(mean), cv),
    stats::pgamma(y, shape = cv^-2, scale = mean * cv^2, log.p = TRUE), tolerance = 1e-13)
  expect_equal(.positive_log_summary("pld_gamma", log(mean), cv)$median_log,
    log(stats::qgamma(.5, shape = cv^-2, scale = mean * cv^2)), tolerance = 1e-13)
  expect_true(is.finite(.positive_log_cdf("gamma", -1000, 0, .5)))
})

test_that("generalized gamma agrees with independent Prentice calculations", {
  y <- log(c(.2, 1, 2, 10, 120))
  m <- log(c(.8, 3, 5, 30, 50))
  sd <- .45
  test_ref <- function(Q) {
    k <- Q^-2
    center <- m - lgamma(k + sd / Q) + lgamma(k) + sd / Q * log(k)
    lg <- log(k) + Q * (y - center) / sd
    density <- stats::dgamma(exp(lg), k, log = TRUE) + lg + log(abs(Q) / sd)
    expect_equal(.gengamma_log_density(y, m, sd, Q), density, tolerance = 1e-11)
    for (lower in c(TRUE, FALSE)) {
      expected <- stats::pgamma(exp(lg), k, lower.tail = xor(lower, Q < 0), log.p = TRUE)
      expect_equal(.positive_log_cdf("gengamma", y, m, sd, Q, lower.tail = lower),
        expected, tolerance = 1e-11)
    }
    expect_equal(.positive_log_summary("gengamma", m, sd, Q)$mean_log,
      center + sd / Q * (digamma(k) - log(k)), tolerance = 1e-12)
  }
  for (Q in c(-1.8, -.7, -.15, .1, .75, 2)) test_ref(Q)
  expect_equal(.gengamma_log_density(y, m, sd, 0),
    stats::dnorm(y, m - sd^2 / 2, sd, log = TRUE))
  expect_equal(.positive_log_cdf("pld_gengamma", y, m, sd, Q = 0),
    .positive_log_cdf("pld_lognormal", y, m, sd), tolerance = 1e-14)
  expect_equal(.gengamma_log_density(y, m, sd, sd),
    .gamma_log_density(y, m, sd), tolerance = 1e-12)
  expect_equal(.positive_log_summary("pld_gengamma", m, sd, sd),
    .positive_log_summary("pld_gamma", m, sd), tolerance = 1e-12)
})

test_that("the Prentice tape retains the signed shape derivatives through zero", {
  y <- log(c(.9, 1.3, 2.5, 6, 10))
  m <- log(c(1, 1.15, 3, 4.5, 8))
  sd <- .6
  fn <- function(Q) sum(.gengamma_log_density(y, m, sd, Q))
  tape <- RTMB::MakeTape(fn, 0)
  w <- (y - m + sd^2 / 2) / sd
  w_prime <- -1/2 - sd^2/6
  w_second <- sd/2 + sd^3/6
  first <- sum(-w * w_prime - w^3/6)
  second <- sum(-1/6 - w_prime^2 - w * w_second - w^2 * w_prime - w^4/12)
  expect_equal(as.numeric(tape$jacobian(0)), first, tolerance = 1e-12)
  expect_equal(as.numeric(tape$jacfun()$jacobian(0)), second, tolerance = 1e-12)
  for (Q in c(-.4, -1e-8, 1e-8, .6)) {
    expect_equal(tape(Q), fn(Q), tolerance = 1e-13)
    expect_true(all(is.finite(tape$jacobian(Q))))
  }
  expect_equal(as.numeric(tape$jacfun()$jacobian(-1e-8)), second, tolerance = 1e-6)
  expect_equal(as.numeric(tape$jacfun()$jacobian(1e-8)), second, tolerance = 1e-6)
})

test_that("estimated shapes obey the finite mean domain and start at zero", {
  dat <- list(gengamma = list(group = c(0L, 1L, 1L, 2L), series = list(
    list(type = "index", name = "RV", fixed = FALSE, Q = NULL),
    list(type = "catch", name = "Fishery", fixed = TRUE, Q = -.7))))
  expect_identical(.gengamma_parameters(dat, "index"), c(RV = 0))
  expect_length(.gengamma_parameters(dat, "catch"), 0L)
  shape <- function(eta) .gengamma_shape(dat, list(eta_gengamma_Q_index = c(RV = eta)),
    c(.2, .4, .4, .6))
  expect_equal(shape(0), c(0, 0, 0, -.7))
  expect_lt(shape(-500)[2L], 0)
  expect_gt(1 + .4 * shape(-500)[2L], 0)
  expect_gt(shape(500)[2L], 0)
  expect_true(is.finite(.gengamma_log_density(0, log(2), .4, shape(-500)[2L])))
})

test_that("Gamma and generalized gamma simulation preserve conditional means", {
  withr::local_seed(137)
  n <- 20000L
  m <- 3
  for (case in list(list(f = "gamma", Q = 0), list(f = "gengamma", Q = -.35),
      list(f = "pld_gengamma", Q = .55))) {
    log_draw <- .draw_positive_log_obs(case$f, rep(log(m), n), .42, case$Q)
    summary <- .positive_log_summary(case$f, log(m), .42, case$Q)
    expect_true(all(is.finite(log_draw)))
    expect_lt(abs(mean(exp(log_draw)) - m), .045)
    expect_lt(abs(mean(log_draw) - summary$mean_log), .012)
    expect_lt(abs(stats::sd(log_draw) - summary$sd_log), .012)
  }
})
