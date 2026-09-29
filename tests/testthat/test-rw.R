test_that("RW density contains only independent temporal increments", {
  x <- matrix(c(2, 2.5, 2.3, -4, -3, -2.8), 3, 2)
  increments <- x[-1, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  expected <- sum(dnorm(increments, 0, .4, log = TRUE))
  expect_equal(dprocess_rw(x, .4), expected)
  expect_equal(dprocess_rw(sweep(x, 2, c(100, -50), "+"), .4), expected)
  expect_equal(dprocess_rw(x, .4),
               sum(vapply(seq_len(ncol(x)), function(a) dprocess_rw(x[, a, drop = FALSE], .4), numeric(1))))
  expect_equal(dprocess_rw(x[1, , drop = FALSE], .4), 0)
  obj <- RTMB::MakeADFun(function(p) -dprocess_rw(p$x, .4), list(x = x), silent = TRUE)
  expect_equal(obj$fn(obj$par), -expected)
  expect_equal(as.numeric(obj$gr(obj$par)[c(1, 4)]), -increments[1, ] / .4^2)
})

test_that("RW N and M fits support warm starts, projections and shorter year ranges", {
  for (process in c("N", "M")) {
    fit <- update(default_fit, silent = TRUE,
      N_settings = list(process = if (process == "N") "rw" else "off", init = "exp"),
      F_settings = list(process = "rw", mu_form = ~ 1),
      M_settings = list(process = if (process == "M") "rw" else "off",
                        mu_supplied = ~ I(.3), age_breaks = c(3, 14)),
      start_par = as.list(default_fit$sdrep, "Estimate"))
    expect_true(is.finite(fit$opt$objective))
    expect_true(all(is.finite(fit$rep$N)))
    expect_false("log_mu_f" %in% names(fit$obj$par))
    expect_false(any(grepl("logit_phi", names(fit$obj$par))))
    projected <- update(fit, years = head(YEARS, -1),
      proj_settings = list(n_proj = 2, n_mean = 1, F_mult = 1),
      start_par = as.list(fit$sdrep, "Estimate"), silent = TRUE)
    states <- as.list(projected$sdrep, "Estimate")
    expect_true(is.finite(projected$opt$objective))
    expect_equal(rownames(states$log_f), as.character(head(YEARS, -1)))
    if (process == "N") expect_equal(rownames(states$log_n), as.character(projected$dat$years[-1]))
    if (process == "M") expect_equal(rownames(states$log_m), as.character(projected$dat$M_settings$years))
    set.seed(271)
    sim <- sim_tam(projected, n = 1, par_uncertainty = "none", redraw_random = TRUE, progress = FALSE)
    expect_true(all(is.finite(sim$N$est)))
  }
})

test_that("RW simulation retains starting rows and accumulates Normal increments", {
  x <- matrix(99, 6, 3, dimnames = list(2000:2005, 2:4))
  x[1, ] <- c(-2, 0, 3)
  set.seed(429)
  expected <- t(replicate(nrow(x) - 1L, rnorm(ncol(x), 0, .2)))
  set.seed(429)
  draw <- rprocess_rw(x, .2)
  expect_identical(draw[1, ], x[1, ])
  expect_identical(dimnames(draw), dimnames(x))
  expect_equal(unname(draw[-1, ] - draw[-nrow(draw), ]), expected)
  expect_identical(rprocess_rw(x[1, , drop = FALSE], .2), x[1, , drop = FALSE])
  expect_equal(dim(rprocess_rw(matrix(1, 4, 1), .2)), c(4L, 1L))
})

test_that("renamed AR1 helpers retain their separable Gaussian mathematics", {
  x <- matrix(c(.2, -.1, .3, -.4, .5, .6), 3, 2)
  phi <- c(.3, .6)
  variance <- .4^2 / prod(1 - phi^2)
  covariance <- variance * kronecker(toeplitz(phi[1]^(0:1)), toeplitz(phi[2]^(0:2)))
  expected <- -.5 * (length(x) * log(2 * pi) + as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
                       drop(t(c(x)) %*% solve(covariance, c(x))))
  expect_equal(dprocess_ar1(x, phi, .4), expected, tolerance = 1e-10)
  set.seed(49)
  expected_draw <- matrix(MASS::mvrnorm(1L, rep(0, length(x)), covariance), 3, 2)
  set.seed(49)
  expect_identical(rprocess_ar1(3, 2, phi, .4), expected_draw)
})

test_that("RW is accepted for N/F/M without AR parameters; retired process is rejected", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    N_settings = list(process = "rw", init = "exp"),
    F_settings = list(process = "rw"), M_settings = list(process = "rw", mu_supplied = ~ I(.3)))
  expect_false(any(grepl("logit_phi", names(make_par(dat)))))
  expect_false(any(grepl("logit_phi", names(dat))))
  expect_identical(make_test_dat()$F_settings$process, "rw")
  for (setting in c("N_settings", "F_settings", "M_settings")) {
    args <- setNames(list(list(process = "approx_rw", mu_supplied = ~ I(.3))), setting)
    expect_error(do.call(make_test_dat, args), "arg")
  }
})

test_that("likelihood and simulation use RW residuals and preserve model boundaries", {
  density <- dprocess_rw
  simulator <- rprocess_rw
  densities <- draws <- list()
  normal <- stats::rnorm
  normal_calls <- list()
  local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    normal_calls[[length(normal_calls) + 1L]] <<- list(mean = mean, sd = sd)
    normal(n, mean, sd)
  }, .package = "stats")
  local_mocked_bindings(
    dprocess_rw = function(x, sd = 1) {
      densities[[length(densities) + 1L]] <<- x
      density(x, sd)
    },
    rprocess_rw = function(x, sd = 1) {
      out <- simulator(x, sd)
      draws[[length(draws) + 1L]] <<- list(start = x[1, ], result = out)
      out
    }, .package = "tinyAM")
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    N_settings = list(process = "rw", init = "free"),
    F_settings = list(process = "rw", mu_form = ~ 1 + I(year - 2000)),
    M_settings = list(process = "rw", mu_form = ~ 0 + I(year - 2000),
      mu_supplied = ~ I(.3), first_dev_year = 2002, age_breaks = c(3, 6)),
    proj_settings = list(n_proj = 2, n_mean = 1, F_mult = c(.8, 1.2)))
  par <- make_par(dat)
  par$log_mu_f[] <- c(-2, .1)
  par$mu_m[] <- .05
  par$log_f[] <- -1
  par$log_m[] <- -1.5
  par$log_n[] <- 2
  par$log_sd_f <- par$log_sd_m <- par$log_sd_n <- log(.1)
  expect_true(is.finite(nll_fun(par, dat)))
  expect_equal(lapply(densities, dim), list(dim(par$log_n), dim(par$log_m), dim(par$log_f)))
  set.seed(744)
  sims <- nll_fun(par, dat, simulate = TRUE)
  expect_equal(sims$log_f[1, ], par$log_f[1, ])
  expect_equal(sims$log_m[1, ], par$log_m[1, ])
  expect_equal(sims$log_n[1, ], par$log_n[1, ])
  expect_equal(sims$log_n0, par$log_n0)
  simulated_densities <- tail(densities, 3) # N, M, F
  expect_equal(unname(simulated_densities[[1]]), unname(draws[[3]]$result))
  expect_equal(unname(simulated_densities[[2]]), unname(draws[[2]]$result))
  expect_equal(unname(simulated_densities[[3]]), unname(draws[[1]]$result))
  observation_calls <- tail(normal_calls, 2)
  shared <- intersect(names(par), names(sims))
  par[shared] <- sims[shared]
  rep <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expect_equal(unname(rep$F[dat$is_proj, ]),
    unname(sweep(exp(sims$log_f[rep(nrow(sims$log_f), 2), , drop = FALSE]), 1, c(.8, 1.2), "*")))
  expect_equal(rep$M[as.character(2000:2001), ], rep$mu_M[as.character(2000:2001), ])
  expect_equal(rep$M[, "2"], rep$mu_M[, "2"])
  expect_equal(observation_calls[[1]]$mean, rep$log_pred[dat$is_observed])
  expect_equal(observation_calls[[1]]$sd, rep$sd_obs[dat$is_observed])
  expect_equal(observation_calls[[2]]$mean, rep$log_pred[dat$fill_missing_map])
})
