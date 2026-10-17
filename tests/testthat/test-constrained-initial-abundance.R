test_that("mono and capped N0 follow their equations, including ties and the plus group", {
  for (init in c("mono", "capped")) for (process in c("off", "iid", "rw", "ar1")) {
    dat <- make_test_dat(years = 2000:2004, ages = 2:6,
      N_settings = list(init = init, process = process))
    par <- make_par(dat)
    par$log_r0 <- log(200)
    par$d_n0[] <- c(0, .3, 0, .1)
    expected <- par$log_r0 - if (init == "mono") cumsum(par$d_n0) else par$d_n0
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
    rep <- obj$report()
    expect_equal(unname(log(rep$N[1, -1])), unname(expected))
    expect_equal(unname(rep$initial_log_N), unname(expected))
    expect_equal(rep$N[1, 2], 200)
    expect_true(all(rep$N[1, -1] <= 200))
    if (init == "mono") expect_true(all(diff(rep$N[1, ]) <= 0))
    expect_true(all(is.finite(obj$gr())))
    set.seed(6)
    sim <- nll_fun(par, dat, simulate = TRUE)
    par[intersect(names(sim), names(par))] <- sim[intersect(names(sim), names(par))]
    rebuilt <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
    expect_equal(rebuilt$N[1, ], rep$N[1, ])
  }
})

test_that("constrained starts align by age and bounds reject invalid starts", {
  dat <- make_test_dat(years = 2000:2004, ages = 2:6, N_settings = list(init = "mono"))
  par <- make_par(dat)
  start <- list(d_n0 = c(`4` = .7, `99` = .8))
  merged <- tinyAM:::.merge_start_par(par, start)
  expect_equal(merged$d_n0["4"], c(`4` = .7))
  expect_identical(names(merged$d_n0), as.character(3:6))
  par$d_n0[1] <- -.1
  expect_error(fit_tam(cod_obs, years = 2000:2004, ages = 2:6,
    N_settings = list(init = "mono"), start_par = par, silent = TRUE), "non-negative")
})

test_that("constrained initial abundances have derived uncertainty and fixed decrements", {
  for (init in c("mono", "capped")) {
    fit <- suppressWarnings(update(default_fit, N_settings = list(init = init), silent = TRUE))
    expect_false(any(c("log_r0", "d_n0", "log_n0") %in% fit$obj$env$.random))
    expect_true(all(fit$bounds$lower[names(fit$obj$par) == "d_n0"] == 0))
    tab <- subset(fit$fixed_par, par == "n0")
    expect_equal(tab$age, fit$dat$ages[-1])
    expect_equal(tab$est, unname(fit$rep$N[1, -1]))
    expect_equal(tab$se_scale, rep("log", nrow(tab)))
    expect_true(all(is.finite(tab$se)))
    drops <- subset(fit$fixed_par, par == "d_n0")
    groups <- tinyAM:::.fixed_parameter_groups(fit$fixed_par)
    expect_equal(unique(groups$`Initial steps`$par), "d_n0")
    expect_true(all(drops$est >= 0))
    expect_equal(drops$se_scale, rep("reported", nrow(drops)))
  }
})
