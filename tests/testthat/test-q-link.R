test_that("q_link defaults to log and rejects invalid links", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6)
  explicit <- make_test_dat(years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~q_block, sd_form = ~1, q_link = "log", fill_missing = TRUE))
  expect_identical(dat$index_settings$q_link, "log")
  expect_identical(make_par(dat), make_par(explicit))
  expect_identical(nll_fun(make_par(dat), dat), nll_fun(make_par(explicit), explicit))
  set.seed(182); a <- nll_fun(make_par(dat), dat, simulate = TRUE)
  set.seed(182); b <- nll_fun(make_par(explicit), explicit, simulate = TRUE)
  expect_identical(a, b)
  for (link in list("identity", "logi", c("log", "logit"), NA_character_, 1)) {
    expect_error(make_test_dat(years = 2000:2005, ages = 2:6,
      index_settings = list(q_form = ~1, sd_form = ~1, q_link = link, fill_missing = TRUE)), "q_link")
  }
})

test_that("logit bounds the full formula predictor without changing its design", {
  obs <- cod_obs
  obs$index$x <- (obs$index$year - 2000) / 5
  obs$index$group <- factor(ifelse(obs$index$age <= 3, "young", "older"))
  settings <- list(q_form = ~group + x, sd_form = ~1, q_link = "logit", fill_missing = TRUE)
  dat <- make_test_dat(data = obs, years = 2000:2005, ages = 2:6,
                       index_settings = settings)
  expect_identical(dat$q_modmat, model.matrix(~group + x, dat$obs$index))
  par <- make_par(dat)
  expect_false("log_q" %in% names(par))
  expect_equal(par$logit_q, setNames(rep(0, 3), colnames(dat$q_modmat)))
  par$logit_q[] <- c(-1, 2, 3)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  report <- obj$report()
  expected <- plogis(drop(dat$q_modmat %*% par$logit_q))
  expect_equal(unname(exp(report$log_q_obs)), unname(expected))
  expect_true(all(expected > 0 & expected < 1))
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
  rows <- cbind(match(dat$obs$index$year, dat$years), match(dat$obs$index$age, dat$ages))
  expect_equal(unname(report$log_pred[dat$obs_map$type == "index"]),
    unname(log(expected) + log(report$N)[rows] -
             report$Z[rows] * dat$obs$index$samp_time))
})

test_that("logit log predictions and gradients stay finite at extreme predictors", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~1, sd_form = ~1, q_link = "logit", fill_missing = TRUE))
  for (eta in c(-1000, 1000)) {
    par <- make_par(dat)
    par$logit_q[] <- eta
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
    expect_true(is.finite(obj$fn(obj$par)))
    expect_true(all(is.finite(obj$gr(obj$par))))
    expect_equal(unname(obj$report()$log_q_obs),
                 rep(plogis(eta, log.p = TRUE), nrow(dat$obs$index)))
  }
})

test_that("mono retains non-decreasing q and exact plateaus with the logit link", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    index_settings = list(q_form = ~mono(q_block), sd_form = ~1, q_link = "logit", fill_missing = TRUE))
  par <- make_par(dat)
  par$logit_q[] <- -2
  par$dq[] <- rep_len(c(0, .5), length(par$dq))
  eta <- drop(dat$q_modmat %*% par$logit_q + dat$q_mono_modmat %*% par$dq)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  q <- exp(obj$report()$log_q_obs)
  expect_equal(unname(q), unname(plogis(eta)))
  by_age <- tapply(q, dat$obs$index$age, mean)
  expect_true(all(diff(by_age) >= 0))
  expect_true(any(diff(by_age) == 0))
  expect_true(all(q > 0 & q < 1))
})

test_that("simulation and tidied observations use the same logit catchability", {
  dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    F_settings = list(process = "iid", mu_form = ~1),
    index_settings = list(q_form = ~q_block, sd_form = ~1, q_link = "logit", fill_missing = TRUE))
  par <- make_par(dat)
  par$logit_q[] <- seq_along(par$logit_q) / 10
  calls <- list()
  testthat::local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    calls[[length(calls) + 1L]] <<- mean
    rep_len(mean, n)
  }, .package = "stats")
  sim <- nll_fun(par, dat, simulate = TRUE)
  par[intersect(names(par), names(sim))] <- sim[intersect(names(par), names(sim))]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  report <- obj$report()
  expect_equal(tail(calls, 2)[[1]], report$log_pred[dat$is_observed])
  fit <- structure(list(dat = dat, rep = report), class = c("tam_fit", "list"))
  expect_equal(tidy_obs_pred(fit)$index$q,
               unname(plogis(drop(dat$q_modmat %*% par$logit_q))))
})

test_that("logit q coefficients and SEs are reported on their fitted scale", {
  par <- list(logit_q = c(baseline = -1, slope = 2))
  obj <- RTMB::MakeADFun(function(p) sum((p$logit_q - par$logit_q)^2) / .02,
                         par, silent = TRUE)
  fit <- structure(list(obj = obj, sdrep = RTMB::sdreport(obj)),
                   class = c("tam_fit", "list"))
  tab <- tidy_par(fit)$fixed
  expect_identical(tab$par, rep("logit_q", 2))
  expect_equal(tab$est, unname(par$logit_q))
  expect_equal(tab$se, rep(.1, 2))
  expect_identical(tab$se_scale, rep("logit", 2))
  expect_equal(tab$lwr, tab$est - qnorm(.975) * tab$se)
  expect_equal(tab$upr, tab$est + qnorm(.975) * tab$se)
})

test_that("fits, updates and parameter draws retain the selected q link", {
  fit <- fit_tam(cod_obs, years = 2000:2005, ages = 2:6, silent = TRUE,
    index_settings = list(q_form = ~1, sd_form = ~1, q_link = "logit", fill_missing = TRUE))
  expect_true(is.finite(fit$opt$objective))
  expect_true(all(fit$obs_pred$index$q > 0 & fit$obs_pred$index$q <= 1))
  expect_true("logit_q" %in% names(as.list(fit$sdrep, "Estimate")))
  refit <- update(fit, start_par = as.list(fit$sdrep, "Estimate"), silent = TRUE)
  expect_identical(refit$dat$index_settings$q_link, "logit")
  expect_equal(refit$obs_pred$index$q, fit$obs_pred$index$q, tolerance = 1e-4)
  draw <- tinyAM:::.draw_none(fit)
  draw$logit_q[] <- 5
  report <- RTMB::MakeADFun(function(p) nll_fun(p, fit$dat), draw, silent = TRUE)$report()
  expect_equal(unname(exp(report$log_q_obs)), rep(plogis(5), nrow(fit$dat$obs$index)))
  testthat::local_mocked_bindings(rnorm = function(n, mean = 0, sd = 1) {
    rep_len(mean, n)
  }, .package = "stats")
  sims <- sim_tam(fit, n = 1, par_uncertainty = "none", redraw_random = FALSE,
                  progress = FALSE, seed = 812)
  expect_equal(sims$index$obs, unname(fit$obs_pred$index$pred), tolerance = 1e-6)
})
