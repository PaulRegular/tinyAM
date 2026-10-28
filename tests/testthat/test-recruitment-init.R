test_that("recruitment boundary integration is opt-in and adds no density", {
  for (bad in list(NULL, "random", NA_character_, c("free", "latent"))) {
    expect_error(make_test_dat(rec_settings = list(init = bad)), "init")
  }
  for (process in c("rw", "iid", "ar1", "off")) {
    d <- make_test_dat(years = 2000:2006, rec_settings = list(
      process = process, init = "latent", mu_form = if (process == "rw") NULL else ~1))
    p <- make_par(d)
    fixed <- d; fixed$rec_settings$init <- "free"
    expect_identical(make_par(fixed), p)
    expect_equal(nll_fun(p, d), nll_fun(p, fixed))
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    expect_true(all(is.finite(obj$gr())))
    sim <- nll_fun(p, d, simulate = TRUE)
    expect_identical(sim$log_r0, p$log_r0)
    expect_false(any(grepl("sd.*init|sd.*r0", names(p))))
  }
  d <- make_test_dat(years = 2000:2006, N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "off", init = "latent", mu_form = ~1))
  d$.likelihood_component <- "process"
  p <- make_par(d); before <- nll_fun(p, d)
  p$log_r0 <- 10
  expect_equal(nll_fun(p, d), before)
})

test_that("stock-recruit boundary years retain their values without extra penalties", {
  d <- make_test_dat(years = 2000:2010, rec_settings = list(
    process = "iid", init = "latent", mu_form = ~bh(ssb, lag = 4)))
  p <- make_par(d)
  expect_equal(names(p$log_r_init), as.character(2001:2003))
  p$log_r_init[] <- c(1, 2, 3)
  fixed <- d; fixed$rec_settings$init <- "free"
  expect_equal(nll_fun(p, d), nll_fun(p, fixed))
  sim <- nll_fun(p, d, simulate = TRUE)
  expect_identical(sim$log_r_init, p$log_r_init)
})

test_that("fits classify initial recruitment correctly and preserve warm starts", {
  legacy <- default_fit
  legacy$refit_args$rec_settings$init <- NULL
  call <- update(legacy, rec_settings = list(process = "iid", mu_form = ~1), evaluate = FALSE)
  expect_identical(call$rec_settings$init, "free")
  f <- update(default_fit, rec_settings = list(init = "latent"),
    start_par = default_fit, silent = TRUE)
  expect_true("log_r0" %in% f$obj$env$.random)
  expect_false("r0" %in% f$fixed_par$par)
  expect_equal(f$random_par$log_r0$year, f$dat$years[1L])
  expect_equal(f$random_par$log_r0$est, f$pop$recruitment$est[1L])
  expect_equal(length(f$opt$par), length(default_fit$opt$par) - 1L)
  starts <- tinyAM:::.start_parameters(f, default_fit$dat)
  expect_equal(starts$log_r0, as.list(f$sdrep, "Estimate")$log_r0)
  expect_warning(tinyAM:::.aic_warnings(f), "improper")
  par <- f$obj$env$parList(par = f$parameter_values)
  sim <- tinyAM:::.sim_obs(f, function(fit) par, redraw_random = TRUE)
  expect_equal(subset(sim$N, year == min(f$dat$years) & age == min(f$dat$ages))$est,
    exp(par$log_r0))
  explicit <- update(default_fit, rec_settings = list(init = "free"),
    start_par = default_fit, silent = TRUE)
  expect_equal(explicit$opt$objective, default_fit$opt$objective, tolerance = 1e-5)
  projected <- update(f, proj_settings = list(n_proj = 1, F_mult = 1), start_par = f, silent = TRUE)
  expect_true("log_r0" %in% projected$obj$env$.random)
  expect_true(any(projected$pop$recruitment$is_proj))
})
