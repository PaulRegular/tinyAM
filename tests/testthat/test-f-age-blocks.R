test_that("explicit F boundaries share states and exclude unfished ages", {
  d <- make_test_dat(ages = 2:6, F_settings = list(age_breaks = c(2, 3, 4, 5, 7)))
  expect_equal(as.character(d$F_settings$age_blocks), c("2", "3", "4", "5-6", "5-6"))
  expect_equal(ncol(make_par(d)$log_f), 4L)
  expect_equal(names(dimnames(make_par(d)$log_f)), c("year", "age_block"))
  for (bad in list(c(2, 7, 7), c(1, 7), c(2, 6), c(2.5, 7), c(7), c(2, NA, 7))) {
    expect_error(make_test_dat(ages = 2:6, F_settings = list(age_breaks = bad)), "age_breaks")
  }
  expect_error(make_test_dat(ages = 2:6, F_settings = list(age_breaks = c(3:5, 7))), "Positive catch")
  obs <- cod_obs; obs$catch$obs[obs$catch$age == 2] <- 0
  d <- make_test_dat(data = obs, ages = 2:6, F_settings = list(age_breaks = c(3:5, 7)))
  p <- make_par(d); p$log_f[] <- log(.2)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  expect_true(is.finite(obj$fn()))
  expect_true(all(is.finite(obj$gr())))
  r <- obj$report()
  expect_equal(unname(r$F[, 1L]), rep(0, length(d$years)))
  expect_equal(r$F[, "5"], r$F[, "6"])
  expect_true(all(exp(r$log_pred[d$obs_map$type == "catch" & d$obs_map$age == 2]) == 0))
  sim <- nll_fun(p, d, simulate = TRUE)
  expect_equal(ncol(sim$log_f), 3L)
})

test_that("F block means and SDs describe unique states", {
  expect_error(make_test_dat(ages = 2:6, F_settings = list(process = "iid", mu_form = ~age,
    age_breaks = c(2, 4, 7))), "F age block.*mu_form")
  expect_error(make_test_dat(ages = 2:6, F_settings = list(sd_form = ~age,
    age_breaks = c(2, 4, 7))), "F age block.*ages")
  for (process in c("iid", "rw", "ar1", "cor_rw", "off")) {
    d <- make_test_dat(ages = 2:6, F_settings = list(process = process,
      mu_form = if (process %in% c("rw", "cor_rw")) NULL else ~1,
      age_breaks = c(2, 4, 7)))
    p <- make_par(d)
    expect_equal(if (process == "off") length(p$log_f) else ncol(p$log_f), if (process == "off") 0L else 2L)
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    expect_true(all(is.finite(obj$gr())))
    if (process != "off") expect_equal(nrow(d$process_sd$F$matrix), 2L)
  }
})

test_that("explicit M boundaries retain legacy blocks and cutting helpers", {
  legacy <- make_test_dat(ages = 2:6, M_settings = list(process = "iid", age_breaks = c(3, 5, 6)))
  expect_identical(legacy$M_settings$age_blocks, cut_ages(3:6, c(3, 5, 6)))
  equivalent <- make_test_dat(ages = 2:6, M_settings = list(process = "iid", age_breaks = c(3, 5, 6, 7)))
  expect_identical(legacy$M_settings$age_blocks, equivalent$M_settings$age_blocks)
  expect_equal(nll_fun(make_par(legacy), legacy), nll_fun(make_par(equivalent), equivalent))
  explicit <- make_test_dat(ages = 2:6, M_settings = list(process = "iid", age_breaks = c(3, 5, 7)))
  expect_equal(as.character(explicit$M_settings$age_blocks), c("3-4", "3-4", "5-6", "5-6"))
  expect_error(cut_ages(2:6, c(2, 7)), "last break")
})

test_that("a single F block has no unused age correlation parameter", {
  for (process in c("cor_rw", "ar1")) {
    d <- make_test_dat(years = 2000:2005, ages = 2:6,
      F_settings = list(process = process, mu_form = if (process == "ar1") ~1 else NULL,
        age_breaks = c(2, 7)))
    p <- make_par(d)
    expect_null(p$atanh_rho_f)
    if (process == "ar1") expect_identical(names(p$logit_phi_f), "year")
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    expect_true(all(is.finite(obj$gr())))
    if (process == "ar1") {
      i <- which(names(obj$par) == "logit_phi_f")
      a <- b <- obj$par; a[i] <- a[i] + 1e-5; b[i] <- b[i] - 1e-5
      expect_equal(as.numeric(obj$gr()[i]), (obj$fn(a) - obj$fn(b)) / 2e-5, tolerance = 1e-6)
    }
  }
})

test_that("F default and explicit singleton densities agree", {
  d <- make_test_dat(ages = 2:6, F_settings = list(process = "cor_rw"))
  b <- make_test_dat(ages = 2:6, F_settings = list(process = "cor_rw", age_breaks = 2:7))
  p <- make_par(d); q <- make_par(b)
  set.seed(91); p$log_f[] <- rnorm(length(p$log_f)); q$log_f[] <- p$log_f
  expect_equal(nll_fun(p, d), nll_fun(q, b))
  a <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  z <- RTMB::MakeADFun(function(p) nll_fun(p, b), q, silent = TRUE)
  expect_equal(a$gr(), z$gr())
})

test_that("unique block densities agree with independent multivariate calculations", {
  d <- make_test_dat(years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off"),
    F_settings = list(process = "cor_rw", age_breaks = c(2, 4, 7)))
  p <- make_par(d); p$log_sd_f <- log(.3); p$atanh_rho_f <- atanh(.4)
  set.seed(32); p$log_f[] <- rnorm(length(p$log_f), sd = .2)
  prior <- d; prior$.likelihood_component <- "process"
  zero <- p; zero$log_f[] <- 0
  delta <- p$log_f[-1L, ] - p$log_f[-nrow(p$log_f), ]
  covariance <- .3^2 * matrix(c(1, .4, .4, 1), 2L)
  expected <- .5 * sum((delta %*% solve(covariance)) * delta)
  expect_equal(nll_fun(p, prior) - nll_fun(zero, prior), expected)
})

test_that("blocked fits retain unique states in downstream workflows", {
  obs <- cod_obs; obs$catch$obs[obs$catch$age == 2] <- 0
  f <- fit_tam(obs, years = 1983:2024, ages = 2:14,
    N_settings = list(process = "iid", init = "exp"), rec_settings = list(init = "free"),
    F_settings = list(process = "cor_rw", age_breaks = c(3:12, 15)),
    M_settings = list(process = "off", mu_supplied = ~I(.3)), silent = TRUE)
  expect_true(all(f$pop$F$est[f$pop$F$age == 2] == 0))
  expect_true("age_block" %in% names(f$random_par$log_f))
  expect_equal(ncol(f$obj$env$parList(par = f$parameter_values)$log_f), 10L)
  d <- tinyAM:::.tidy_process_deviations(f)
  expect_true("age_block" %in% names(d$F))
  expect_equal(nrow(d$F), 41L * 10L)
  expect_true(all(is.finite(d$F$est)))
  tab <- tidy_tam(model_list = list(blocked = f), process_residuals = FALSE)
  expect_s3_class(plot_process_deviations(tab), "plotly")
  projected <- update(f, proj_settings = list(n_proj = 1, catch = 100), start_par = f, silent = TRUE)
  expect_true(all(projected$pop$F$est[projected$pop$F$age == 2] == 0))
  expect_equal(tail(projected$pop$total_yield_pred$est, 1), 100, tolerance = 1e-6)
  starts <- tinyAM:::.start_parameters(f, projected$dat)
  expect_equal(starts$log_f, f$obj$env$parList(par = f$parameter_values)$log_f)
  changed <- make_test_dat(data = obs, F_settings = list(age_breaks = c(3:11, 15)))
  starts <- tinyAM:::.start_parameters(f, changed)
  expect_false("12-14" %in% colnames(starts$log_f))
  expect_warning(AIC(f), "improper")
  ca <- suppressWarnings(cAIC(f, details = TRUE))
  expect_equal(ca$edf_random, ca$edf_rec + ca$edf_N + ca$edf_F + ca$edf_M + ca$edf_other)
  path <- tempfile(fileext = ".html")
  vis_tam(model_list = list(blocked = f, ordinary = default_fit), output_file = path,
    open_file = FALSE, render_args = list(quiet = TRUE))
  expect_true(file.exists(path))
})
