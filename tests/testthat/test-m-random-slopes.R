m_slope_dat <- function(formula = ~ 0 + capelin + iid(M_group, by = capelin),
                        years = 1983:1994, projection = NULL, data = NULL,
                        breaks = c(2, 4, 6, 8, 10)) {
  if (is.null(data)) {
    data <- cod_obs
    d <- data$pop_covariates
    d$M_group <- factor(pmin(3, pmax(0, floor((d$age - 2) / 2))))
    d$capelin <- sin((d$year - 1983) / 2)
    data$pop_covariates <- d
  }
  prepare_tam(data, years = years, ages = 2:10,
    N_settings = list(process = "off"),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "ar1", mu_supplied = ~ I(.3),
      mu_form = formula, age_breaks = breaks, first_dev_year = min(years)),
    proj_settings = projection)
}

test_that("AR1 M accepts only one age-defined numeric-by IID slope", {
  dat <- m_slope_dat()
  p <- make_par(dat)
  expect_identical(dim(p$log_m), c(12L, 4L))
  term <- dat$M_terms[[1L]]
  expect_length(p[[term$parameter]], 4L)
  expect_true(term$parameter %in% tinyAM:::.formula_random_parameters(dat))
  expect_length(make_par(m_slope_dat(~ 0 + iid(M_group, by = capelin)))$mu_m, 0L)
  expect_true(is.finite(nll_fun(p, dat)))
  expect_identical(dim(make_par(m_slope_dat(~ 0 + iid(age, by = capelin),
    breaks = 2:10))$log_m), c(12L, 9L))
  coarse <- dat$obs
  coarse$pop_covariates$M_group <- factor(as.integer(coarse$pop_covariates$age >= 6))
  expect_length(make_par(m_slope_dat(data = coarse))[[term$parameter]], 2L)
  for (formula in list(~ 0 + iid(year, by = capelin), ~ 0 + iid(M_group),
    ~ 0 + rw(year, by = capelin), ~ 0 + ar1(year, by = capelin),
    ~ 0 + iid(M_group, by = capelin) + iid(age, by = capelin))) {
    expect_error(m_slope_dat(formula), "one iid\\(age_group")
  }
  bad <- dat$obs
  bad$pop_covariates$M_group <- factor(paste(bad$pop_covariates$age, bad$pop_covariates$year))
  expect_error(m_slope_dat(data = bad), "constant across years")
  expect_error(prepare_tam(dat$obs, years = dat$years, ages = dat$ages,
    F_settings = list(process = "ar1", mu_form = ~ 0 + iid(M_group, by = capelin))),
    "IID residuals")
})

test_that("M block errors identify fixed, supplied and random mean conflicts", {
  dat <- m_slope_dat(~ 0 + capelin:factor(age))
  expect_error(make_par(dat), "mu_form.*differs between ages 2 and 3 in year 1984")
  expect_error(make_par(m_slope_dat(~ 0 + iid(age, by = capelin))),
    "mu_form M_iid_age_capelin differs")
  dat <- m_slope_dat()
  dat$log_mu_supplied_m <- log(.3) + .01 * dat$pop_data$age
  expect_error(make_par(dat), "mu_supplied differs between ages 2 and 3")
  expect_error(make_par(dat), "Align mean-effect groups.*finer age_breaks")
})

test_that("M slope cautions count independent groups rather than duplicated ages", {
  dat <- m_slope_dat()
  notes <- tinyAM:::.mean_process_advisories(dat)
  expect_true(all(c("M_random_slopes_AR1", "M_random_slope_groups") %in% notes$issue))
  expect_match(notes$detail[notes$issue == "M_random_slope_groups"], "Only 4 independent")
  dat$M_terms[[1L]]$sd_parameter <- NULL
  expect_identical(tinyAM:::.mean_process_advisories(dat)$issue, "M_random_slopes_AR1")
  dat$M_settings$process <- "off"
  expect_equal(nrow(tinyAM:::.mean_process_advisories(dat)), 0L)
})

test_that("random slope integration matches the unique-state AR1 covariance", {
  dat <- m_slope_dat(years = 1983:1989)
  term <- dat$M_terms[[1L]]
  rows <- dat$pop_data$age %in% dat$M_settings$age_block_start
  z <- tinyAM:::.q_term_design(term, nrow(dat$pop_data))[rows, , drop = FALSE]
  x <- dat$pop_data$capelin[rows]
  y <- log(.3) + .1 * x + seq(-.15, .15, length.out = length(x))
  p <- make_par(dat)[c("mu_m", term$parameter, term$sd_parameter, "log_sd_m", "logit_phi_m")]
  p$mu_m[] <- .1
  p[[term$sd_parameter]][] <- log(.15)
  p$log_sd_m[] <- log(.12)
  p$logit_phi_m[] <- qlogis(c(.25, .6))
  ny <- length(dat$years)
  nb <- ncol(z)
  obj <- RTMB::MakeADFun(function(p) {
    effect <- tinyAM:::.mean_effects(p, dat, "M")
    residual <- matrix(y - log(.3) - x * p$mu_m - effect$contribution[rows], ny, nb)
    effect$nll - dprocess_ar1(residual, phi = plogis(p$logit_phi_m), sd = exp(p$log_sd_m))
  }, p, random = term$parameter, silent = TRUE)
  covariance <- .15^2 * tcrossprod(z) + .12^2 / ((1 - .25^2) * (1 - .6^2)) *
    kronecker(.25^abs(outer(seq_len(nb), seq_len(nb), `-`)),
      .6^abs(outer(seq_len(ny), seq_len(ny), `-`)))
  residual <- y - log(.3) - .1 * x
  ch <- chol(covariance)
  exact <- length(y) / 2 * log(2 * pi) + sum(log(diag(ch))) +
    sum(backsolve(ch, residual, transpose = TRUE)^2) / 2
  expect_equal(as.numeric(obj$fn(obj$par)), exact, tolerance = 1e-8)
  derivative <- vapply(seq_along(obj$par), function(j) {
    lo <- hi <- obj$par
    lo[j] <- lo[j] - 1e-5
    hi[j] <- hi[j] + 1e-5
    (obj$fn(hi) - obj$fn(lo)) / 2e-5
  }, numeric(1))
  expect_equal(as.numeric(obj$gr(obj$par)), derivative, tolerance = 1e-5)
})

test_that("simulation reuses slopes and preserves absolute block states", {
  dat <- m_slope_dat(projection = list(n_proj = 2, n_mean = 2, F_mult = 1))
  p <- make_par(dat)
  term <- dat$M_terms[[1L]]
  set.seed(9302)
  sim <- nll_fun(p, dat, simulate = TRUE)
  expect_length(sim[[term$parameter]], 4L)
  p[intersect(names(p), names(sim))] <- sim[intersect(names(p), names(sim))]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  report <- obj$report()
  mu <- matrix(log(.3) + drop(dat$M_modmat %*% p$mu_m) +
    tinyAM:::.mean_effects(p, dat, "M")$contribution, length(dat$years), length(dat$ages))
  expect_equal(log(report$mu_M), mu, ignore_attr = TRUE)
  expect_equal(report$M, exp(p$log_m[, dat$M_settings$age_blocks]), ignore_attr = TRUE)
  expect_equal(report$M[, 1], report$M[, 2], ignore_attr = TRUE)
  expect_equal(dat$pop_data$capelin[dat$pop_data$year == max(dat$years)],
    dat$pop_data$capelin[dat$pop_data$year == max(dat$years) - 2L])
  baseline <- dat
  baseline$M_terms <- list()
  p[[term$parameter]][] <- 0
  expected <- tinyAM:::.mean_effects(p, dat, "M")$nll
  expect_equal(as.numeric(nll_fun(p, dat) - nll_fun(p, baseline)), expected, tolerance = 1e-9)
  short <- m_slope_dat(years = 1983:1989)
  aligned <- tinyAM:::.merge_start_par(make_par(short), p)
  expect_equal(aligned[[term$parameter]], p[[term$parameter]])
  expect_equal(rownames(aligned$log_m), as.character(1983:1989))
})

test_that("M random-slope fitting warns before optimization", {
  dat <- m_slope_dat()
  local_mocked_bindings(nlminb = function(...) stop("stop before optimizer"), .package = "stats")
  expect_warning(expect_error(fit_tam(dat$obs, years = dat$years, ages = dat$ages,
    F_settings = dat$F_settings[c("process", "mu_form")],
    M_settings = dat$M_settings[c("process", "mu_form", "mu_supplied", "age_breaks", "first_dev_year")],
    silent = TRUE), "stop before optimizer"), "Formula-effect caution")
})
