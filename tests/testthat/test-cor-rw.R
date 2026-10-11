test_that("correlated RW uses marginal increment SDs and the full Gaussian density", {
  x <- matrix(c(.1, .3, -.2, .4, .2, .8, -.1, .5, .7), 3, 3)
  sd <- c(.2, .4, .7)
  for (rho in c(-.6, 0, .7)) {
    covariance <- outer(sd, sd) * rho^abs(outer(1:3, 1:3, "-"))
    increments <- x[-1, ] - x[-3, ]
    expected <- sum(apply(increments, 1, function(z) -.5 * (
      3 * log(2 * pi) + as.numeric(determinant(covariance)$modulus) +
        sum(z * solve(covariance, z)))))
    expect_equal(tinyAM:::.dprocess_cor_rw(x, sd, rho), expected, tolerance = 1e-10)
    obj <- RTMB::MakeADFun(function(p) -tinyAM:::.dprocess_cor_rw(p$x, exp(p$sd), tanh(p$rho)),
      list(x = x, sd = log(sd), rho = atanh(rho)), silent = TRUE)
    expect_equal(obj$fn(), -expected, tolerance = 1e-10)
    expect_true(all(is.finite(obj$gr())))
  }
  expect_equal(tinyAM:::.dprocess_cor_rw(x, .3, 0), dprocess_rw(x, .3))
  expect_equal(tinyAM:::.dprocess_cor_rw(x[1, , drop = FALSE], sd, .6), 0)
})

test_that("correlated RW simulation preserves starts and increment covariance", {
  x <- matrix(0, 10001, 3, dimnames = list(NULL, c("2", "3", "4")))
  x[1, ] <- c(-1, 2, 3)
  sd <- c(.2, .4, .7)
  rho <- -.5
  set.seed(162)
  draw <- tinyAM:::.rprocess_cor_rw(x, sd, rho)
  expect_identical(draw[1, ], x[1, ])
  expect_identical(dimnames(draw), dimnames(x))
  covariance <- outer(sd, sd) * rho^abs(outer(1:3, 1:3, "-"))
  expect_lt(max(abs(cov(draw[-1, ] - draw[-nrow(draw), ]) - covariance)), .015)
})

test_that("correlated F RW is generative and retains the independent limit", {
  dat <- make_test_dat(F_settings = list(process = "cor_rw"))
  par <- make_par(dat)
  expect_true("atanh_rho_f" %in% names(par))
  independent <- dat
  independent$F_settings$process <- "rw"
  expect_equal(nll_fun(par, dat), nll_fun(par[setdiff(names(par), "atanh_rho_f")], independent))
  par$atanh_rho_f <- atanh(.4)
  set.seed(29)
  draw <- nll_fun(par, dat, simulate = TRUE)
  par[intersect(names(par), names(draw))] <- draw[intersect(names(par), names(draw))]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  expect_equal(obj$report()$F[!dat$is_proj, ], exp(draw$log_f))
  expect_true(is.finite(obj$fn()))
})

test_that("correlated F supports fitted uncertainty, warm starts and folds", {
  fit <- suppressWarnings(update(default_fit, F_settings = list(process = "cor_rw"), silent = TRUE))
  expect_true(is.finite(fit$opt$objective))
  expect_true("rho_f" %in% fit$fixed_par$par)
  refit <- suppressWarnings(update(fit, years = head(YEARS, -1),
    start_par = as.list(fit$sdrep, "Estimate"),
    proj_settings = list(n_proj = 1, n_mean = 1, F_mult = 1), silent = TRUE))
  expect_true(is.finite(refit$opt$objective))
  expect_equal(refit$rep$F[nrow(refit$rep$F), ], refit$rep$F[nrow(refit$rep$F) - 1, ])
  expect_error(make_test_dat(M_settings = list(process = "cor_rw")), "arg")
})
