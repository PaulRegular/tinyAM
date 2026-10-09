test_that("stationary F processes require an explicit estimated mean", {
  for (process in c("iid", "ar1")) {
    for (settings in list(list(process = process),
                          list(process = process, mu_form = NULL))) {
      expect_error(make_test_dat(F_settings = settings), "mu_form is required")
    }
    dat <- make_test_dat(F_settings = list(process = process, mu_form = ~ 1))
    par <- make_par(dat)
    expect_identical(names(par$log_mu_f), "(Intercept)")
    expect_true(is.finite(nll_fun(par, dat)))
  }
  expect_null(make_test_dat(F_settings = list(process = "rw", mu_form = NULL))$F_settings$mu_form)
})

test_that("F means must have complete, full-rank historical designs", {
  expect_error(make_test_dat(F_settings = list(process = "iid", mu_form = ~ 0)),
               "finite, non-empty")
  expect_error(make_test_dat(F_settings = list(process = "ar1", mu_form = ~ age + I(2 * age))),
               "rank-deficient")
  expect_error(make_test_dat(years = 2000:2005,
    F_settings = list(process = "iid", mu_form = ~ I(year > 2005)),
    proj_settings = list(n_proj = 2, n_mean = 1, F_mult = 1)), "rank-deficient")
  expect_error(make_test_dat(
    F_settings = list(process = "iid", mu_form = ~ I(log(age - 2)))),
    "every modeled year and age")
})

test_that("time-invariant F mean columns cancel from RW densities", {
  plain <- make_test_dat(years = 2000:2005, ages = 2:6)
  expect_warning(blocked <- make_test_dat(years = 2000:2005, ages = 2:6,
    F_settings = list(process = "rw", mu_form = ~ factor(age))), "time-invariant")
  expect_null(blocked$F_settings$mu_form)
  expect_identical(make_par(plain), make_par(blocked))
  par <- make_par(plain)
  par$log_f[] <- seq_along(par$log_f) / 30
  expect_identical(nll_fun(par, plain), nll_fun(par, blocked))
  set.seed(531); a <- nll_fun(par, plain, simulate = TRUE)
  set.seed(531); b <- nll_fun(par, blocked, simulate = TRUE)
  expect_identical(a, b)
})

test_that("RW means retain identifiable drift and reject temporal redundancy", {
  expect_warning(dat <- make_test_dat(years = 2000:2005, ages = 2:6,
    F_settings = list(process = "rw", mu_form = ~ factor(age) + I(year - 2000))),
    "time-invariant")
  expect_identical(colnames(dat$F_modmat), "I(year - 2000)")
  par <- make_par(dat)
  expect_length(par$log_mu_f, 1L)
  par$log_mu_f[] <- .2
  par$log_f[] <- -2
  report <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
  expected <- matrix(rep(.2 * (0:5), 5), 6, 5,
                     dimnames = dimnames(report$mu_F))
  expect_equal(log(report$mu_F), expected)
  density <- dprocess_rw
  residual <- NULL
  local_mocked_bindings(dprocess_rw = function(x, sd = 1) {
    residual <<- x
    density(x, sd)
  })
  nll_fun(par, dat)
  expect_equal(residual, par$log_f - expected)
  expect_error(make_test_dat(F_settings = list(process = "rw",
    mu_form = ~ 0 + year + I(year + age))), "not identifiable from RW increments")
})
