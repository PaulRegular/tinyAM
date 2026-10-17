n0_report <- function(par, dat) {
  RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
}

test_that("recruitment has a fixed first-year anchor and subsequent RW states", {
  dat <- make_test_dat(years = 2000:2005, N_settings = list(process = "off", init = "exp"))
  par <- make_par(dat)
  par$log_r0 <- log(200)
  increments <- c(0.2, -0.3, 0.1, 0.4, -0.1)
  par$log_r[] <- par$log_r0 + cumsum(increments)
  rep <- n0_report(par, dat)
  expect_equal(unname(log(rep$N[1, 1])), par$log_r0)
  expect_identical(names(par$log_r), as.character(dat$years[-1]))
  expect_equal(as.numeric(diff(log(rep$recruitment))), increments)
  expect_equal(names(rep$recruitment), as.character(dat$years))
  no_obs <- dat
  no_obs$observed <- numeric()
  no_obs$is_observed[] <- FALSE
  constant <- par
  constant$log_r[] <- constant$log_r0
  expect_equal(
    nll_fun(par, no_obs) - nll_fun(constant, no_obs),
    -sum(dnorm(increments, 0, 1, log = TRUE) - dnorm(0, 0, 1, log = TRUE))
  )
  constant$log_r0 <- constant$log_r0 + 10
  constant$log_r[] <- constant$log_r0
  expect_equal(nll_fun(constant, no_obs), nll_fun(par, no_obs) - sum(increments^2) / 2)
  next_increment <- 0L
  testthat::local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    if (n == 1L) {
      next_increment <<- next_increment + 1L
      return(mean + increments[next_increment])
    }
    rep_len(mean, n)
  }, .package = "stats")
  sim <- nll_fun(par, dat, simulate = TRUE)
  expect_equal(sim$log_r, par$log_r)
  expect_identical(names(sim$log_r), names(par$log_r))
})

test_that("each initializer is independent of the subsequent N process", {
  for (init in c("exp", "free", "latent")) {
    initial <- list()
    for (process in c("off", "iid", "rw", "ar1")) {
      dat <- make_test_dat(
        years = 2000:2003,
        N_settings = list(process = process, init = init)
      )
      par <- make_par(dat)
      par$log_r0 <- log(80)
      par$log_f[1, ] <- log(seq(0.1, 0.4, length.out = length(dat$ages)))
      if (init != "exp") par$log_n0[] <- seq(-0.2, 0.3, length.out = length(par$log_n0))
      rep <- n0_report(par, dat)
      initial[[process]] <- rep$N[1, ]
      expect_equal(unname(log(rep$N[1, 1])), par$log_r0)
      if (init == "exp") {
        expect_equal(unname(log(rep$N[1, -1])), par$log_r0 - cumsum(unname(rep$Z[1, -ncol(rep$Z)])))
      } else {
        expect_equal(unname(log(rep$N[1, -1])), unname(par$log_n0))
      }
      if (init == "exp") expect_false(any(c("log_n0", "log_sd_n0") %in% names(par)))
      if (process != "off") {
        expect_equal(dim(par$log_n), c(length(dat$years) - 1L, length(dat$ages) - 1L))
        expect_equal(rownames(par$log_n), as.character(dat$years[-1]))
      }
    }
    for (process in names(initial)[-1]) expect_equal(initial[[process]], initial$off)
  }
})

test_that("free initialization spans arbitrary initial older-age abundance without a penalty", {
  dat <- make_test_dat(
    years = 2000:2003,
    N_settings = list(process = "off", init = "free")
  )
  par <- make_par(dat)
  target <- log(c(100, 400, 50, 250, 80, 30, 90, 12, 17, 5, 10, 3))
  par$log_n0[] <- target
  expect_length(par$log_n0, length(dat$ages) - 1L)
  expect_equal(log(unname(n0_report(par, dat)$N[1, -1])), target)
  expect_false("log_sd_n0" %in% names(par))
  # With observations removed and N deterministic, no likelihood term depends
  # on the initial older ages. Free states must add no penalty.
  dat$observed <- numeric()
  dat$is_observed[] <- FALSE
  nll <- nll_fun(par, dat)
  par$log_n0[] <- 0
  expect_equal(nll_fun(par, dat), nll)
})

test_that("latent initialization has the free joint likelihood without an initial density", {
  dat <- make_test_dat(years = 2000:2003, N_settings = list(process = "iid", init = "latent"))
  par <- make_par(dat)
  par$log_n0[] <- seq(-.4, .7, length.out = length(par$log_n0))
  free <- dat
  free$N_settings$init <- "free"
  expect_equal(nll_fun(par, dat), nll_fun(par, free))
  expect_equal(n0_report(par, dat)$N, n0_report(par, free)$N)
  expect_false("log_sd_n0" %in% names(par))
})


test_that("latent initial abundance is retained under process simulation", {
  dat <- make_test_dat(years = 2000:2003,
    N_settings = list(process = "iid", init = "latent"),
    F_settings = list(process = "iid", mu_form = ~1))
  par <- make_par(dat)
  par$log_n0[] <- log(seq(100, 20, length.out = length(par$log_n0)))
  set.seed(734)
  sim <- nll_fun(par, dat, simulate = TRUE)
  expect_identical(sim$log_n0, par$log_n0)
  expect_false(isTRUE(all.equal(sim$log_f, par$log_f)))
  par[intersect(names(sim), names(par))] <- sim[intersect(names(sim), names(par))]
  expect_equal(unname(log(n0_report(par, dat)$N[1, -1])), unname(par$log_n0))
})


test_that("N states represent only genuine cohort transitions with a plus group", {
  dat <- make_test_dat(years = 2000:2003, N_settings = list(process = "iid", init = "exp"))
  par <- make_par(dat)
  deterministic <- dat
  deterministic$N_settings$process <- "off"
  baseline <- n0_report(par, deterministic)
  par$log_n[] <- log(baseline$N[-1, -1])
  par$log_sd_n <- log(0.3)
  expect_equal(n0_report(par, dat)$N, baseline$N)
  expect_equal(
    nll_fun(par, dat) - nll_fun(par, deterministic),
    -length(par$log_n) * dnorm(0, 0, 0.3, log = TRUE)
  )
  testthat::local_mocked_bindings(
    rnorm = function(n, mean = 0, sd = 1) rep_len(mean, n) + 0.2,
    .package = "stats"
  )
  sims <- nll_fun(par, dat, simulate = TRUE)
  par[intersect(names(sims), names(par))] <- sims[intersect(names(sims), names(par))]
  rep <- n0_report(par, dat)
  expect_equal(unname(log(rep$N[-1, -1])), unname(sims$log_n))
  for (y in 2:length(dat$years)) {
    survivors <- rep$N[y - 1, ] * exp(-rep$Z[y - 1, ])
    expected <- survivors[-length(survivors)]
    expected[length(expected)] <- tail(expected, 1) + tail(survivors, 1)
    expect_equal(unname(log(rep$N[y, -1]) - log(expected)), rep(0.2, length(expected)))
  }
})

test_that("initialization validates options and supports two ages", {
  expect_error(make_test_dat(ages = 2, N_settings = list(init = "latent")), "at least two modeled ages")
  expect_error(make_test_dat(N_settings = list(init = "random")), "arg")
  for (init in c("exp", "free", "latent", "mono", "capped")) {
    dat <- make_test_dat(ages = 2:3, years = 2000:2003, N_settings = list(process = "off", init = init))
    par <- make_par(dat)
    expect_true(is.finite(nll_fun(par, dat)))
    expect_equal(unname(n0_report(par, dat)$N[1, 1]), exp(par$log_r0))
    expect_false("log_sd_n0" %in% names(par))
  }
})


test_that("fits, summaries, updates and conditional simulations honor N0 roles", {
  for (init in c("free", "latent")) {
    fit <- update(default_fit, N_settings = list(process = "iid", init = init), silent = TRUE)
    par <- as.list(fit$sdrep, "Estimate")
    expect_true(is.finite(fit$opt$objective))
    expect_identical("log_n0" %in% fit$obj$env$.random, init == "latent")
    expect_false(any(c("log_r0", "log_sd_n0") %in% fit$obj$env$.random))
    expect_equal(fit$random_par$log_r$year, fit$dat$years[-1])
    expect_equal(fit$pop$recruitment$year, fit$dat$years)
    tab <- if (init == "latent") fit$random_par$log_n0 else subset(fit$fixed_par, par == "n0")
    expect_equal(tab$age, fit$dat$ages[-1])
    expect_equal(unname(tab$est), unname(exp(par$log_n0)))
    for (nm in "r0") {
      expect_equal(fit$fixed_par$est[fit$fixed_par$par == nm], exp(par[[paste0("log_", nm)]]))
    }
    # A controlled random draw isolates the redraw switch without refitting.
    par$log_sd_f <- par$log_sd_n <- par$log_sd_r <- log(0.05)
    fit$dat$logit_phi_f <- qlogis(c(0, 0))
    for (redraw in c(FALSE, TRUE)) {
      set.seed(18)
      sim <- tinyAM:::.sim_obs(fit, function(fit) par, redraw_random = redraw)
      n0 <- subset(sim$N, year == min(fit$dat$years))$est
      expect_equal(log(n0[-1]), unname(par$log_n0))
    }
    dat <- make_test_dat(
      years = 1990:2020, ages = 2:20,
      N_settings = list(process = "iid", init = init)
    )
    start <- tinyAM:::.merge_start_par(make_par(dat), par)
    expect_equal(names(start$log_r), as.character(dat$years[-1]))
    expect_equal(rownames(start$log_n), as.character(dat$years[-1]))
    expect_equal(start$log_n0[names(par$log_n0)], par$log_n0)
  }
})

test_that("warm starts align single-element state vectors by year or age", {
  par <- list(log_r0 = 0, log_r = c(`2001` = 0), log_n0 = c(`4` = 0))
  start <- list(log_r0 = 2, log_r = c(`2002` = 3), log_n0 = c(`5` = 4))
  merged <- tinyAM:::.merge_start_par(par, start)
  expect_equal(merged$log_r0, 2)
  expect_identical(merged$log_r, par$log_r)
  expect_identical(merged$log_n0, par$log_n0)
})
