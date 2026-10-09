mean_test_dat <- function(component, type, ...) {
  formula <- switch(type, iid = ~ iid(year, sd = .2), rw = ~ rw(year, sd = .2),
                    ar1 = ~ ar1(year, sd = .2, phi = .6))
  data <- lapply(cod_obs, function(d) { d$age_group <- factor(d$age); d })
  args <- list(data = data, years = 1983:1989, ages = 2:6,
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3)))
  args[[paste0(component, "_settings")]]$mu_form <- formula
  if (component == "M") args$M_settings$mu_form <- update(formula, ~ 0 + .)
  args <- utils::modifyList(args, list(...))
  do.call(prepare_tam, args)
}

test_that("F/M formulas reuse unique densities and keep absolute latent states", {
  for (component in c("F", "M")) for (type in c("iid", "rw", "ar1")) {
    dat <- mean_test_dat(component, type)
    p <- make_par(dat)
    term <- dat[[paste0(component, "_terms")]][[1]]
    expect_true(term$parameter %in% names(p))
    expect_true(term$parameter %in% tinyAM:::.formula_random_parameters(dat))
    expect_length(p[[term$parameter]], length(dat$years) - as.integer(type == "rw"))
    baseline <- dat
    baseline[[paste0(component, "_terms")]] <- list()
    expected <- tinyAM:::.mean_effects(p, dat, component)$nll
    expect_equal(as.numeric(nll_fun(p, dat) - nll_fun(p, baseline)),
                 as.numeric(expected), tolerance = 1e-9)
    p[[term$parameter]][] <- .1
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
    expect_true(all(is.finite(obj$gr(obj$par))))
    report <- obj$report()
    contribution <- matrix(tinyAM:::.mean_effects(p, dat, component)$contribution,
                           length(dat$years), length(dat$ages))
    expect_equal(as.numeric(log(report[[paste0("mu_", component)]])),
      as.numeric(contribution + if (component == "M") log(.3) else 0))
  }
})

test_that("mean by terms distinguish multiplied and independent group effects", {
  dat <- mean_test_dat("F", "rw", F_settings = list(mu_form =
    ~ factor(age) + rw(year, by = age, sd = .2)))
  p <- make_par(dat)
  term <- dat$F_terms[[1]]
  expect_length(term$groups, 1L)
  p[[term$parameter]][] <- -.1
  effect <- tinyAM:::.mean_effects(p, dat, "F")
  d <- dat$obs$catch
  expect_equal(effect$contribution, ifelse(d$year == min(d$year), 0, -.1 * d$age))
  dat <- mean_test_dat("F", "ar1", F_settings = list(mu_form =
    ~ factor(age) + ar1(year, by = age_group)))
  expect_length(dat$F_terms[[1]]$groups, length(dat$ages))
  p <- make_par(dat)
  expect_length(p$log_sd_mu_F_ar1_year_age_group, 1L)
  expect_length(p$logit_phi_mu_F_ar1_year_age_group, 1L)
})

test_that("mean/residual variance aliases and temporal overlaps are rejected", {
  expect_error(mean_test_dat("F", "rw", F_settings = list(process = "rw")), "IID residuals")
  expect_error(mean_test_dat("M", "ar1", M_settings = list(process = "ar1")), "IID residuals")
  expect_error(mean_test_dat("F", "rw", F_settings = list(
    mu_form = ~ rw(year) + ar1(year))), "one temporal")
  expect_error(mean_test_dat("F", "iid", F_settings = list(
    mu_form = ~ iid(year, by = factor_age))), "existing observation column")
  expect_error(mean_test_dat("F", "iid", F_settings = list(
    mu_form = ~ iid(year, by = age_group))), "IID F residual variance")
  expect_error(mean_test_dat("M", "iid", M_settings = list(process = "iid",
    age_breaks = 2:6, first_dev_year = 1983,
    mu_form = ~ 0 + iid(year, by = age_group))), "Unreplicated IID M")
  expect_error(mean_test_dat("F", "iid", F_settings = list(
    mu_form = ~ factor(age) + iid(age))), "saturated")
  expect_error(mean_test_dat("M", "iid", M_settings = list(mu_form = ~ mono(age))), "catchability curves")
  dat <- mean_test_dat("M", "rw", M_settings = list(process = "iid",
    mu_form = ~ 0 + rw(year, by = age), age_breaks = c(3, 6)))
  expect_error(make_par(dat), "within M age_blocks")
})

test_that("simulation constructs observations from the returned mean and residual states", {
  dat <- mean_test_dat("F", "ar1", M_settings = list(process = "iid",
    mu_form = ~ 0 + rw(year, sd = .1), mu_supplied = ~ I(.3), age_breaks = 3:6))
  p <- make_par(dat)
  p$log_r0 <- log(1e5)
  p$log_sd_r <- log(.1)
  p$log_sd_f <- log(.1)
  p$log_sd_m <- log(.05)
  p$log_mu_f[] <- log(.2)
  p$log_f[] <- log(.2)
  p$log_sd_catch[] <- p$log_sd_index[] <- log(.1)
  set.seed(9001)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  set.seed(9001)
  expect_equal(simulated, nll_fun(p, dat, simulate = TRUE))
  keep <- intersect(names(p), names(simulated))
  p[keep] <- simulated[keep]
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
  report <- obj$report()
  expect_equal(report$F, exp(p$log_f), ignore_attr = TRUE)
  mu <- matrix(tinyAM:::.mean_effects(p, dat, "F")$contribution,
               length(dat$years), length(dat$ages)) + log(.2)
  expect_equal(report$mu_F, exp(mu), ignore_attr = TRUE)
  standardized <- (simulated$log_obs[dat$is_observed] - report$log_pred[dat$is_observed]) /
    report$sd_obs[dat$is_observed]
  expect_lt(abs(mean(standardized)) * sqrt(length(standardized)), 3)
  expect_equal(sd(standardized), 1, tolerance = .2)
  expect_true(all(c(tinyAM:::.formula_random_parameters(dat)) %in% names(simulated)))
})

test_that("mean states align by names and forecast densities integrate out", {
  dat <- mean_test_dat("F", "rw", proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1))
  p <- make_par(dat)
  term <- dat$F_terms[[1]]
  p[[term$parameter]][] <- seq_along(p[[term$parameter]]) / 10
  short <- mean_test_dat("F", "rw", years = 1983:1987)
  merged <- tinyAM:::.merge_start_par(make_par(short), p)
  expect_equal(merged[[term$parameter]], p[[term$parameter]][names(merged[[term$parameter]])])
  n <- length(dat$years)
  objective <- function(dat) {
    p <- make_par(dat)
    RTMB::MakeADFun(function(p) {
      e <- tinyAM:::.mean_effects(p, dat, "F")
      e$nll - sum(RTMB::dnorm(rep(.1, sum(!dat$is_proj)),
        e$contribution[which(!dat$is_proj)], .1, log = TRUE))
    }, p[term$parameter], random = term$parameter, silent = TRUE)
  }
  future <- objective(dat)
  past <- objective(mean_test_dat("F", "rw"))
  expect_equal(future$fn(future$par), past$fn(past$par), tolerance = 1e-9)
  draw <- future$env$parList()
  level <- tinyAM:::.mean_effects(draw, dat, "F")$contribution[seq_len(n)]
  expect_equal(level[(n - 1):n], rep(level[n - 2], 2), tolerance = 1e-7)
})

test_that("integrated mean likelihood equals the exact Gaussian covariance calculation", {
  for (type in c("iid", "rw", "ar1")) {
    dat <- mean_test_dat("F", type)
    term <- dat$F_terms[[1]]
    p <- make_par(dat)[term$parameter]
    year <- match(dat$obs$catch$year, dat$years) - 1L
    covariance <- switch(type,
      iid = .2^2 * outer(year, year, `==`),
      rw = .2^2 * outer(year, year, pmin),
      ar1 = .2^2 / (1 - .6^2) * .6^abs(outer(year, year, `-`)))
    diag(covariance) <- diag(covariance) + .1^2
    y <- seq(-.1, .1, length.out = length(year))
    L <- chol(covariance)
    expected <- length(y) / 2 * log(2 * pi) + sum(log(diag(L))) +
      sum(forwardsolve(t(L), y)^2) / 2
    obj <- RTMB::MakeADFun(function(p) {
      e <- tinyAM:::.mean_effects(p, dat, "F")
      e$nll - sum(RTMB::dnorm(y, e$contribution, .1, log = TRUE))
    }, p, random = term$parameter, silent = TRUE)
    expect_equal(as.numeric(obj$fn(obj$par)), expected, tolerance = 1e-9)
  }
})

test_that("F/M mean effects fit, report signed effects and retain empty fixed designs", {
  fit <- update(default_fit, silent = TRUE, years = 1983:2000,
    F_settings = list(process = "iid", mu_form = ~ factor(age) + rw(year, sd = .1)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3),
                      mu_form = ~ 0 + ar1(year, sd = .05, phi = .6)))
  expect_true(is.finite(fit$opt$objective))
  expect_true(all(c("eta_mu_F_rw_year", "eta_mu_M_ar1_year") %in% fit$obj$env$.random))
  expect_true(all(c("F_rw_year", "M_ar1_year") %in% names(fit$formula_effects$levels)))
  expect_equal(fit$formula_effects$levels$F_rw_year$est[1], 0)
  expect_equal(fit$formula_effects$levels$M_ar1_year$component, rep("M", 18))
  expect_true(all(is.finite(fit$formula_effects$increments$F_rw_year$se)))
  expect_false(any(fit$fixed_par$par == "mu_m"))
  expect_false(any(grepl("increments", names(fit$pop))))
  sim <- tinyAM:::.sim_obs(fit, tinyAM:::.draw_none, redraw_random = TRUE)
  expect_true(all(is.finite(sim$F$est)))
  tabs <- tidy_tam(model_list = list(mean = fit), interval = .9)
  expect_true(all(c("F", "M") %in% unique(do.call(rbind, tabs$formula_effects$levels)$component)))
})

test_that("mortality design cautions apply only to the estimated components", {
  dat <- mean_test_dat("M", "rw", M_settings = list(process = "iid",
    age_breaks = 2:6, mu_form = ~ 0 + rw(year)))
  caution <- tinyAM:::.mean_process_advisories(dat)
  expect_identical(caution$issue, "M_variance_separation")
  expect_match(caution$detail, "process = 'off'")
  expect_match(caution$detail, "supply sd")
  dat$M_settings$process <- "off"
  expect_equal(nrow(tinyAM:::.mean_process_advisories(dat)), 0)
  dat$M_settings$process <- "iid"
  dat$M_terms[[1]]$sd_parameter <- NULL
  expect_equal(nrow(tinyAM:::.mean_process_advisories(dat)), 0)
  dat <- mean_test_dat("F", "rw", F_settings = list(mu_form = ~ factor(age) + rw(year)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3), mu_form = ~ 0 + ar1(year)))
  caution <- tinyAM:::.mean_process_advisories(dat)
  expect_identical(caution$issue, "joint_mortality_means")
  expect_match(caution$detail, "one mortality mean at a time")
  dat$F_terms[[1]]$sd_parameter <- NULL
  expect_equal(nrow(tinyAM:::.mean_process_advisories(dat)), 0)
})

test_that("fitting warns about M variance separation before optimizing", {
  local_mocked_bindings(nlminb = function(...) stop("stop before optimizer"), .package = "stats")
  expect_warning(expect_error(fit_tam(data = cod_obs, years = 1983:1995, ages = 2:6,
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "iid", mu_form = ~ 0 + rw(year),
      mu_supplied = ~ I(.3), age_breaks = 2:6), silent = TRUE), "stop before optimizer"),
    "Mortality mean-process caution")
})
