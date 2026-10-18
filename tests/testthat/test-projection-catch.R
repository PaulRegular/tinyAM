projection_dat <- function(...) {
  prepare_tam(cod_obs, years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off"), ...)
}
projection_par <- function(dat) {
  par <- make_par(dat)
  par$log_r0 <- log(1e5)
  par$log_r[] <- log(1e5)
  par$log_f[] <- log(.15)
  par
}

test_that("catch controls switch modes and align projection years", {
  a <- .merge_settings(NULL, list(n_proj = 2, catch = c(`2007` = 80, `2006` = 100)), "proj_settings")
  expect_null(a$F_mult)
  d <- projection_dat(proj_settings = a)
  expect_equal(unname(d$proj_settings$catch), c(100, 80))
  b <- .merge_settings(a, list(F_mult = .8), "proj_settings")
  expect_null(b$catch)
  expect_equal(b$F_mult, .8)
  b <- .merge_settings(a, list(catch = NULL), "proj_settings")
  expect_equal(b$F_mult, 1)
  expect_error(.merge_settings(NULL, list(n_proj = 1, catch = 1, F_mult = 1), "proj_settings"), "not both")
  for (bad in list(-1, NA_real_, Inf, "100", numeric())) {
    expect_error(projection_dat(proj_settings = list(n_proj = 1, catch = bad)), "finite")
  }
  expect_error(projection_dat(proj_settings = list(n_proj = 2, catch = c(`2008` = 1, `2009` = 2))), "every projection year")
})

test_that("catch solver agrees with independent equations and derivatives", {
  N <- c(1000, 800, 600); W <- c(1, 2, 3); F <- c(.1, .3, .5); M <- c(.2, .3, .4)
  target <- 300
  independent <- function(k, N) sum(N * W * k * F / (k * F + M) * (1 - exp(-k * F - M))) - target
  k <- .solve_catch(target, N, W, F, M, 2000)
  expect_equal(k, stats::uniroot(independent, c(0, 5), N = N, tol = 1e-11)$root, tolerance = 1e-8)
  obj <- RTMB::MakeADFun(function(p) .solve_catch(target, exp(p$logN), W, F, M, 2000), list(logN = log(N)), silent = TRUE)
  step <- 1e-5
  finite <- vapply(seq_along(N), function(i) {
    a <- b <- log(N); a[i] <- a[i] + step; b[i] <- b[i] - step
    (.solve_catch(target, exp(a), W, F, M, 2000) - .solve_catch(target, exp(b), W, F, M, 2000)) / (2 * step)
  }, numeric(1))
  expect_equal(as.numeric(obj$gr()), finite, tolerance = 1e-6)
  expect_equal(.solve_catch(0, N, W, F, M, 2000), 0)
  expect_error(.solve_catch(sum(N * W), N, W, F, M, 2000), "2000.*infeasible")
  expect_equal(.baranov_fraction(c(0, .2)), c(1, (1 - exp(-.2)) / .2))
  obj <- RTMB::MakeADFun(function(p) {
    k <- .solve_catch(target, exp(p$logN), W, F, M, 2000)
    RTMB::ADREPORT(k)
    sum((p$logN - log(N))^2) / 2
  }, list(logN = log(N)), silent = TRUE)
  sr <- RTMB::sdreport(obj)
  expect_true(all(is.finite(sr$sd) & sr$sd > 0))
})

test_that("catch projections use chronological abundance and retained plus biology", {
  d <- projection_dat(proj_settings = list(n_proj = 3, catch = c(100, 80, 60)))
  p <- projection_par(d)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  expect_equal(unname(r$total_yield_pred[d$is_proj]), c(100, 80, 60), tolerance = 1e-7)
  expect_true(all(is.finite(obj$gr())))
  for (y in which(d$is_proj)) {
    expect_equal(as.numeric(r$F[y, ] / r$F[max(which(!d$is_proj)), ]), rep(r$F[y, 1] / r$F[max(which(!d$is_proj)), 1], length(d$ages)))
    N <- r$N[y - 1, ]
    expect_equal(unname(r$N[y, 2]), unname(N[1] * exp(-r$Z[y - 1, 1])))
    expect_equal(r$N[y, length(d$ages)], sum(tail(N * exp(-r$Z[y - 1, ]), 2)))
  }
  old <- projection_dat(proj_settings = list(n_proj = 3))
  ordinary <- RTMB::MakeADFun(function(p) nll_fun(p, old), p, silent = TRUE)
  expect_equal(obj$fn(), ordinary$fn())
  expect_equal(obj$gr(), ordinary$gr())
  zero <- projection_dat(proj_settings = list(n_proj = 1, catch = 0))
  op <- RTMB::MakeADFun(function(p) nll_fun(p, zero), projection_par(zero), silent = TRUE)
  rz <- op$report()
  expect_equal(as.numeric(rz$F[zero$is_proj, ]), rep(0, length(zero$ages)))
  expect_equal(unname(tail(rz$total_yield_pred, 1)), 0)
  expect_equal(as.numeric(rz$S[zero$is_proj, ]), as.numeric(rz$S[sum(!zero$is_proj), ]))
  expect_true(all(is.finite(op$gr())))
})



test_that("simulated population and catch projections use the same realization", {
  d <- projection_dat(proj_settings = list(n_proj = 2, catch = c(100, 80)))
  p <- projection_par(d)
  set.seed(42)
  simulated <- nll_fun(p, d, simulate = TRUE)
  p$log_r[] <- simulated$log_r
  p$log_f[] <- simulated$log_f
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  expect_equal(unname(r$total_yield_pred[d$is_proj]), c(100, 80), tolerance = 1e-7)
  expect_true(all(is.na(simulated$log_obs[d$obs_map$is_proj])))
  # Complete hidden-age weights are abundance-weighted before solving catch.
  expect_equal(r$W[d$is_proj, ncol(r$W)],
    rowSums(r$N_plus[d$is_proj, , drop = FALSE] * d$W_plus_input[d$is_proj, , drop = FALSE]) /
      r$N[d$is_proj, ncol(r$N)], tolerance = 1e-8)
})


test_that("a fitted zero-catch projection has exact zero F and safe uncertainty", {
  skip_if_not(exists("default_fit"), "Fitted model fixture unavailable")
  fitted <- suppressWarnings(update(default_fit,
    proj_settings = list(n_proj = 1, catch = 0),
    start_par = as.list(default_fit$sdrep, "Estimate"), silent = TRUE))
  expect_equal(fitted$opt$objective, default_fit$opt$objective, tolerance = 1e-6)
  expect_equal(unname(tail(fitted$rep$total_yield_pred, 1)), 0)
  expect_equal(tail(fitted$pop$F_bar$est, 1), 0)
  expect_equal(tail(fitted$pop$F_bar$se, 1), 0)
  expect_equal(tail(fitted$pop$F_bar$se_scale, 1), "reported")
  expect_true(all(is.finite(fitted$gradient)))
})


test_that("catch projections reject spawning-time zero-lag circularity", {
  data <- cod_obs
  data$maturity$obs[data$maturity$age == 2] <- 0
  data$maturity$spawn_time <- .25
  expect_error(prepare_tam(data, years = 2000:2008, ages = 2:6,
    rec_settings = list(process = "iid", mu_form = ~ bh(ssb, lag = 0)),
    proj_settings = list(n_proj = 1, catch = 1)), "zero-lag.*circular")
  data$maturity$spawn_time <- 0
  expect_no_error(prepare_tam(data, years = 2000:2008, ages = 2:6,
    rec_settings = list(process = "iid", mu_form = ~ bh(ssb, lag = 0)),
    proj_settings = list(n_proj = 1, catch = 1)))
})
