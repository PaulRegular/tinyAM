test_that("process SD formulas align to ages and unique M states", {
  obs <- cod_obs
  obs$weight$age_group <- factor(ifelse(obs$weight$age < 5, "young", "old"), levels = c("young", "old"))
  dat <- make_test_dat(data = obs, years = 1983:1988, ages = 2:6,
    N_settings = list(process = "iid", sd_form = ~ age_group),
    F_settings = list(process = "cor_rw", sd_form = ~ factor(age)),
    M_settings = list(process = "iid", mu_supplied = ~ I(.3), age_breaks = c(3, 5, 6), sd_form = ~ age_group))
  expect_equal(rownames(dat$process_sd$N$matrix), as.character(3:6))
  expect_equal(rownames(dat$process_sd$F$matrix), as.character(2:6))
  expect_equal(nrow(dat$process_sd$M$matrix), length(dat$M_settings$age_block_start))
  par <- make_par(dat)
  expect_false(any(c("log_sd_n", "log_sd_f", "log_sd_m") %in% names(par)))
  par$sd_beta_n <- c(log(.2), log(2))
  expect_equal(unname(tinyAM:::.process_sd(par, dat, "N")), c(.2, .2, .4, .4))
  expect_true(is.finite(nll_fun(par, dat)))
  set.seed(928)
  expect_true(all(is.finite(nll_fun(par, dat, simulate = TRUE)$log_f)))
})

test_that("scaled densities include the change-of-scale Jacobian", {
  x <- matrix(c(.1, .3, -.2, .4, .2, .8), 3, 2)
  sd <- c(.2, .7)
  for (process in c("iid", "rw", "ar1")) {
    z <- if (process == "rw") x[-1, ] - x[-3, ] else x
    phi <- if (process == "ar1") c(.3, .5) else c(0, 0)
    covariance <- kronecker(outer(sd, sd) * phi[1]^abs(outer(1:2, 1:2, "-")),
      phi[2]^abs(outer(seq_len(nrow(z)), seq_len(nrow(z)), "-"))) / prod(1 - phi^2)
    expected <- -.5 * (length(z) * log(2 * pi) + as.numeric(determinant(covariance)$modulus) + sum(z * solve(covariance, c(z))))
    expect_equal(tinyAM:::.dprocess_scaled(x, sd, process, phi), expected, tolerance = 1e-10)
    if (process == "ar1") expect_equal(tinyAM:::.dprocess_scaled(x, rep(.4, 2), process, phi), dprocess_ar1(x, phi, .4))
    if (process == "rw") expect_equal(tinyAM:::.dprocess_scaled(x, rep(.4, 2), process), dprocess_rw(x, .4))
  }
})

test_that("process SD guards reject unsupported or conflicting designs", {
  for (form in list(~ rw(age), ~ iid(age), ~ year, ~ 0, ~ age + I(2 * age))) {
    expect_error(make_test_dat(F_settings = list(process = "rw", sd_form = form)), "sd_form")
  }
  expect_error(make_test_dat(M_settings = list(process = "iid", mu_supplied = ~ I(.3),
    sd_form = ~ factor(age), age_breaks = c(3, 14))), "within an M age block")
  obs <- cod_obs
  obs$weight$varying <- obs$weight$year
  expect_error(make_test_dat(data = obs, F_settings = list(process = "rw", sd_form = ~ varying)), "varies across years")
  obs$weight$varying[] <- NA
  expect_error(make_test_dat(data = obs, F_settings = list(process = "rw", sd_form = ~ varying)), "complete covariates")
  expect_error(make_test_dat(N_settings = list(process = "off", sd_form = ~ age)), "no process")
})

test_that("default SD formulas retain parameter layouts and simulation draws", {
  dat <- make_test_dat()
  explicit <- make_test_dat(N_settings = list(process = "iid", sd_form = ~ 1), F_settings = list(process = "rw", sd_form = ~ 1))
  expect_identical(make_par(dat), make_par(explicit))
  expect_equal(nll_fun(make_par(dat), dat), nll_fun(make_par(explicit), explicit))
  set.seed(826)
  a <- nll_fun(make_par(dat), dat, simulate = TRUE)
  set.seed(826)
  b <- nll_fun(make_par(explicit), explicit, simulate = TRUE)
  expect_identical(a, b)
})

test_that("SD coefficients stay on log SD scale and profiles have RTMB intervals", {
  fit <- suppressWarnings(update(default_fit, F_settings = list(process = "rw", sd_form = ~ age), silent = TRUE))
  coefficients <- fit$fixed_par[fit$fixed_par$par == "sd_beta_f", ]
  expect_equal(coefficients$est, unname(as.list(fit$sdrep, "Estimate")$sd_beta_f))
  expect_true(all(coefficients$se_scale == "log SD coefficient"))
  expect_equal(sort(unique(fit$pop$sd_F$age)), AGES)
  expect_true(all(fit$pop$sd_F$est > 0))
  expect_true(all(c("se", "lwr", "upr") %in% names(fit$pop$sd_F)))
})
