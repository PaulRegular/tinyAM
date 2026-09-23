core_report <- function(par, dat) {
  RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
}

test_that("M states start at their mean, including age blocks and delayed deviations", {
  for (process in c("iid", "ar1")) {
    dat <- make_test_dat(
      years = 2000:2005, ages = 2:6,
      M_settings = list(process = process, mu_form = ~ 0 + I(year - 2000),
                        mu_supplied = ~ I(0.2 + 0.1 * (age >= 4)),
                        age_breaks = c(2, 4, 6), first_dev_year = 2002)
    )
    par <- make_par(dat)
    rep <- core_report(par, dat)
    expect_equal(rep$M, rep$mu_M)
    expect_equal(unname(par$log_m),
                 unname(log(rep$mu_M[rownames(par$log_m), dat$M_settings$age_block_start])))
    expect_true(is.finite(nll_fun(par, dat)))

    par$mu_m[] <- 0.1
    rep <- core_report(par, dat)
    par$log_m[] <- log(rep$mu_M[rownames(par$log_m), dat$M_settings$age_block_start])
    expect_equal(core_report(par, dat)$M, rep$mu_M)
  }
})

test_that("simulated absolute log M is centred on its mean surface", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    F_settings = list(process = "iid"),
    M_settings = list(process = "iid", mu_form = ~ 0 + I(year - 2000),
                      mu_supplied = ~ I(0.3), age_breaks = c(3, 6)))
  par <- make_par(dat)
  par$mu_m[] <- 0.1
  par$log_sd_m <- log(0.2)
  rep <- core_report(par, dat)
  mu <- log(rep$mu_M[rownames(par$log_m), dat$M_settings$age_block_start, drop = FALSE])
  set.seed(704)
  draws <- replicate(200, nll_fun(par, dat, simulate = TRUE)$log_m)
  expect_equal(as.numeric(apply(draws, c(1, 2), mean)), as.numeric(mu), tolerance = 0.04)
})

test_that("simulation uses returned states, recursive cohorts and matching observation SDs", {
  # Fixed innovations make the generative equations and row alignment exact checks.
  calls <- list()
  testthat::local_mocked_bindings(
    rnorm = function(n, mean = 0, sd = 1) {
      calls[[length(calls) + 1L]] <<- list(n = n, mean = mean, sd = sd)
      rep_len(mean, n) + rep_len(sd, n)
    }, .package = "stats")
  testthat::local_mocked_bindings(
    rprocess_2d = function(ny, na, phi = c(0, 0), sd = 1) matrix(0.1, ny, na),
    .package = "tinyAM")

  for (n_process in c("off", "iid", "ar1")) {
    dat <- make_test_dat(years = 2000:2005, ages = 2:6,
      N_settings = list(process = n_process, init_N0 = TRUE),
      F_settings = list(process = "iid"),
      M_settings = list(process = "iid", mu_supplied = ~ I(0.3)),
      proj_settings = list(n_proj = 2, n_mean = 1, F_mult = c(0.8, 1.2)))
    # Missing rows precede observed rows; every row has a distinct SD.
    dat$log_sd_catch_supplied <- log(seq_len(nrow(dat$obs$catch)) / 100)
    dat$log_sd_index_supplied <- log(seq_len(nrow(dat$obs$index)) / 50)
    par <- make_par(dat)
    par$log_sd_r <- log(0.2)
    before <- core_report(par, dat)
    sims <- nll_fun(par, dat, simulate = TRUE)
    obs_calls <- tail(calls, 2)
    par[intersect(names(par), names(sims))] <- sims[intersect(names(par), names(sims))]
    rep <- core_report(par, dat)
    expect_equal(as.numeric(diff(sims$log_r)), rep(0.2, length(dat$years) - 1))
    expect_false(isTRUE(all.equal(rep$log_pred, before$log_pred)))
    expect_equal(obs_calls[[1]]$mean, rep$log_pred[dat$is_observed])
    expect_equal(obs_calls[[1]]$sd, rep$sd_obs[dat$is_observed])
    expect_equal(sims$log_obs[dat$is_observed],
                 unname(rep$log_pred[dat$is_observed] + rep$sd_obs[dat$is_observed]))
    expect_equal(sims$missing,
                 unname(rep$log_pred[dat$fill_missing_map] + rep$sd_obs[dat$fill_missing_map]))
    for (y in 2:nrow(rep$N)) {
      survivors <- rep$N[y - 1, ] * exp(-rep$Z[y - 1, ])
      expected <- survivors[-length(survivors)]
      expected[length(expected)] <- expected[length(expected)] + tail(survivors, 1)
      if (n_process != "off") expected <- expected * exp(0.1)
      expect_equal(unname(rep$N[y, -1]), unname(expected))
    }
  }
})

test_that("IID and AR1 M models fit with absolute mortality states", {
  for (process in c("iid", "ar1")) {
    fit <- update(default_fit,
      N_settings = list(process = "off", init_N0 = TRUE),
      M_settings = list(process = process, mu_supplied = ~ I(0.3), age_breaks = c(3, 14)),
      silent = TRUE)
    expect_true(is.finite(fit$opt$objective))
    expect_true(fit$is_converged)
    par <- as.list(fit$sdrep, "Estimate")
    expect_equal(unname(fit$rep$M[rownames(par$log_m), "3"]),
                 unname(exp(par$log_m[, 1])))
  }
})

test_that("generative predictions remain finite when abundance underflows", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off", init_N0 = TRUE),
    F_settings = list(process = "iid", mu_form = ~ 1))
  par <- make_par(dat)
  par$log_mu_f[] <- log(2000)
  par$log_sd_f <- log(0.01)
  set.seed(510)
  sims <- nll_fun(par, dat, simulate = TRUE)
  expect_true(all(is.finite(sims$log_obs)))
  par[intersect(names(par), names(sims))] <- sims[intersect(names(par), names(sims))]
  report <- core_report(par, dat)
  expect_true(any(report$N == 0))
  expect_true(all(is.finite(report$log_pred)))
})

test_that("retros and every hindcast fold use observed terminal years", {
  fit <- update(default_fit,
    proj_settings = list(n_proj = 3, n_mean = 3, F_mult = c(0.5, 0.7, 0.9)),
    silent = TRUE)
  original <- fit
  original_par <- fit$obj$env$last.par.best
  terminal <- max(fit$dat$years[!fit$dat$is_proj])
  retro <- fit_retro(fit, folds = 1, progress = FALSE)
  expect_equal(as.integer(names(retro$fits)), (terminal - 1):terminal)
  expect_equal(vapply(retro$fits, tinyAM:::.terminal_year, numeric(1)),
               setNames(as.numeric((terminal - 1):terminal), names(retro$fits)))
  hindcast <- fit_hindcast(fit, folds = 1, progress = FALSE)
  expect_equal(names(hindcast$fits), names(retro$fits))
  for (nm in names(hindcast$fits)) {
    fold <- hindcast$fits[[nm]]
    expect_equal(sum(fold$dat$is_proj), 1)
    expect_equal(max(fold$dat$years), as.numeric(nm) + 1)
    expect_equal(fold$dat$proj_settings$n_mean, 1)
    expect_equal(unname(fold$dat$proj_settings$F_mult), 1)
    expect_equal(unname(fold$rep$F[fold$dat$is_proj, ]),
                 unname(fold$rep$F[as.character(nm), ]))
  }
  expect_identical(fit, original)
  expect_identical(fit$obj$env$last.par.best, original_par)
})

test_that("sim_tam keeps observations consistent with redrawn and retained states", {
  fit <- default_fit
  fit$dat$logit_phi_f <- qlogis(c(0, 0))
  par <- as.list(default_fit$sdrep, "Estimate")
  par$log_sd_f <- par$log_sd_n <- par$log_sd_r <- log(0.1)
  par$log_sd_catch[] <- log(1e-8)
  par$log_sd_index[] <- log(1e-8)
  for (redraw in c(FALSE, TRUE)) {
    set.seed(62)
    sim <- tinyAM:::.sim_obs(fit, function(fit) par, redraw_random = redraw)
    catch <- merge(sim$catch, sim$N[, c("year", "age", "est")], by = c("year", "age"))
    names(catch)[names(catch) == "est"] <- "N"
    catch <- merge(catch, sim$F[, c("year", "age", "est")], by = c("year", "age"))
    names(catch)[names(catch) == "est"] <- "F"
    catch <- merge(catch, sim$Z[, c("year", "age", "est")], by = c("year", "age"))
    expected <- with(catch, log(N * F / est * (1 - exp(-est))))
    expect_equal(log(catch$obs), expected, tolerance = 1e-6)
  }
})
