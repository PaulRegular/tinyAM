test_that("M AR1 projection states already satisfy their conditional forecasts", {
  d <- prepare_tam(cod_obs, years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off"), M_settings = list(process = "ar1"),
    proj_settings = list(n_proj = 2))
  p <- make_par(d)
  phi <- .6
  p$logit_phi_m <- qlogis(c(.2, phi))
  p$log_m[] <- log(.2) + .4
  p$log_m[as.character(2006:2007), ] <- log(.2) + .4 * phi^(1:2)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  expect_equal(unname(log(r$M[as.character(2006:2007), "3"]) - log(.2)), .4 * phi^(1:2))
  # Conditional mode of future log-M under the existing separable density.
  x <- matrix(c(.4, .4 * phi, .4 * phi^2), 3, 1)
  den <- function(v) -dprocess_ar1(matrix(c(.4, v), 3, 1), phi = c(0, phi), sd = .3)
  fitted <- stats::optim(c(.1, .1), den, method = "BFGS", control = list(reltol = 1e-12))$par
  expect_equal(fitted, x[-1, 1], tolerance = 1e-5)
})


test_that("simulated AR1 M forecasts retain innovation variance and dependence", {
  d <- prepare_tam(cod_obs, years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off"), M_settings = list(process = "ar1"),
    proj_settings = list(n_proj = 2))
  p <- make_par(d)
  p$logit_phi_m <- qlogis(c(.25, .6))
  p$log_sd_m <- log(.15)
  p$log_sd_f <- p$log_sd_r <- log(.05)
  set.seed(123)
  draws <- replicate(400, {
    x <- nll_fun(p, d, simulate = TRUE)$log_m[, 1] - log(.2)
    c(innovation = x[length(x)] - .6 * x[length(x) - 1], previous = x[length(x) - 1])
  })
  expect_lt(abs(mean(draws[1, ])), .03)
  expected <- .15 / sqrt(1 - .25^2)
  expect_lt(abs(stats::sd(draws[1, ]) - expected), 4 * expected / sqrt(2 * (ncol(draws) - 1)))
  expect_lt(abs(stats::cor(draws[1, ], draws[2, ])), .15)
})
