audit_report <- function(par, dat) {
  RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
}

test_that("mean mortality at a single selected age equals that age's rate", {
  dat <- make_test_dat(years = 2000:2003, ages = 2:6,
    F_settings = list(process = "iid", mean_ages = 3),
    M_settings = list(process = "off", mu_supplied = ~ I(.2), mean_ages = 4))
  rep <- audit_report(make_par(dat), dat)
  expect_equal(unname(rep$F_bar), unname(rep$F[, "3"]))
  expect_equal(unname(rep$M_bar), unname(rep$M[, "4"]))
})

test_that("Baranov catches remain finite at small positive mortality", {
  dat <- make_test_dat(years = 2000:2003, ages = 2:6,
    M_settings = list(process = "off", mu_supplied = ~ I(exp(-50))))
  par <- make_par(dat)
  par$log_f[] <- -50
  rep <- audit_report(par, dat)
  ic <- dat$obs_map$type == "catch"
  # As Z approaches zero, C approaches N * F.
  expect_equal(unname(rep$log_pred[ic]), c(log(rep$N) - 50), tolerance = 1e-12)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
})

test_that("projections retain all catch rows despite recent missing observations", {
  obs <- cod_obs
  obs$catch$obs[obs$catch$year %in% 2003:2005 & obs$catch$age == 2] <- NA
  # A survey may end before the historical model does.
  obs$index <- subset(obs$index, year <= 2004)
  dat <- prepare_tam(data = obs, years = 2000:2005, ages = 2:6,
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 1))
  expect_equal(as.integer(table(dat$obs$catch$year)), rep(5L, 8))
  expect_equal(sort(unique(dat$obs$index$year[dat$obs$index$is_proj])), 2006:2007)
  rep <- NULL
  expect_no_warning(rep <- audit_report(make_par(dat), dat))
  expect_true(all(is.finite(rep$total_catch_pred)))
})

test_that("plus groups preserve separate surveys and available older-age data", {
  obs <- cod_obs
  index <- expand.grid(year = 2000:2001, age = 4:6, survey = c("A", "B"))
  index$obs <- ifelse(index$survey == "A", 1, 10) * index$age
  index$samp_time <- .5
  obs$index <- index
  plus <- tinyAM:::.plus_fun(obs, 4)$index
  expect_equal(plus$obs[plus$survey == "A"], rep(15, 2))
  expect_equal(plus$obs[plus$survey == "B"], rep(150, 2))
  # No original row at age 4 for B: retain the sum of ages 5 and 6.
  obs$index <- subset(index, !(survey == "B" & age == 4))
  plus <- tinyAM:::.plus_fun(obs, 4)$index
  expect_equal(plus$obs[plus$survey == "B"], rep(110, 2))
  expect_true(all(plus$age == 4))
  obs$index$obs[] <- NA_real_
  expect_true(all(is.na(tinyAM:::.plus_fun(obs, 4)$index$obs)))
})

test_that("hindcast scores match observation series and do not require every fold", {
  d <- data.frame(year = 2001, age = 2,
    type = rep(c("catch", "index", "index"), 2),
    survey = rep(c(NA, "A", "B"), 2),
    obs = c(10, 100, 1000, NA, NA, NA),
    pred = c(NA, NA, NA, 10, 100, 1000),
    fold = rep(c(2002, 2000), each = 3), is_proj = rep(c(FALSE, TRUE), each = 3))
  expect_equal(compute_hindcast_rmse(d, log = FALSE), 0)
  expect_equal(compute_hindcast_rmse(d), 0)
  # Repeated observations from retained fits do not multiply their weight.
  repeated <- d[1:3, ]
  repeated$fold <- 2003
  expect_equal(compute_hindcast_rmse(rbind(d, repeated)), 0)
})

test_that("M process blocks require constant mean components, including default blocks", {
  dat <- make_test_dat(years = 2000:2003, ages = 2:6,
    M_settings = list(process = "iid", mu_supplied = ~ I(.1 * age)))
  expect_error(make_par(dat), "M mean structure varies")
  dat$M_settings$process <- "off"
  expect_no_error(make_par(dat))
  # A single arbitrary coefficient combination would hide these opposing slopes.
  dat <- make_test_dat(years = 2000:2003, ages = 2:6,
    M_settings = list(process = "iid", mu_supplied = NULL,
                     mu_form = ~ 0 + age + I(-age / 10), age_breaks = c(3, 6)))
  expect_error(make_par(dat), "M mean structure varies")
})

test_that("invalid grid boundaries cannot enter cohort and mortality matrices", {
  expect_error(make_test_dat(years = 2000), "two historical years")
  for (start in list(2006, 2001.5, c(2000, 2001))) {
    expect_error(make_test_dat(years = 2000:2005,
      M_settings = list(process = "iid", mu_supplied = ~ I(.2), first_dev_year = start)),
      "single historical modeled year")
  }
  for (nm in c("catch", "weight", "maturity")) {
    obs <- cod_obs
    obs[[nm]] <- rbind(obs[[nm]], obs[[nm]][1, ])
    expect_error(check_obs(obs), "exactly one row per year and age")
  }
})

test_that("singleton AR1 axes are not estimated as redundant variance multipliers", {
  captured <- NULL
  local_mocked_bindings(MakeADFun = function(func, parameters, ..., map) {
    captured <<- list(par = parameters, map = map)
    stop("captured parameter map")
  }, .package = "RTMB")
  expect_error(fit_tam(cod_obs, years = 2000:2001, ages = 2:3,
    N_settings = list(process = "ar1", init = "exp"),
    M_settings = list(process = "ar1", mu_supplied = ~ I(.2))), "captured parameter map")
  expect_true(all(is.na(captured$map$logit_phi_n)))
  expect_true(all(is.na(captured$map$logit_phi_m)))
  expect_equal(unname(plogis(captured$par$logit_phi_n)), c(0, 0))
  expect_equal(unname(plogis(captured$par$logit_phi_m)), c(0, 0))
})

test_that("parameter plots allow negative effects and large estimates", {
  d <- data.frame(par = c("mu_m", "r0"), coef = NA_character_,
                  est = c(-.5, 100), lwr = c(-1, 80), upr = c(.2, 120))
  p <- plotly::plotly_build(plot_par(d))
  limits <- p$x$layout$xaxis$range
  expect_true(is.null(limits) || (min(limits) <= -1 && max(limits) >= 120))
})

test_that("missing formula rows and invalid supplied scales cannot misalign likelihoods", {
  obs <- cod_obs
  obs$index$x <- 1
  obs$index$x[which(obs$index$year == 2001)[1]] <- NA_real_
  expect_error(prepare_tam(obs, years = 2000:2003,
    index_settings = list(q_form = ~ x, sd_form = ~ 1, fill_missing = TRUE)), "one finite row")
  obs$weight$x <- .2
  obs$weight$x[which(obs$weight$year == 2001)[1]] <- NA_real_
  expect_error(prepare_tam(obs, years = 2000:2003,
    M_settings = list(process = "off", mu_supplied = ~ x)), "positive finite value")
  expect_error(make_test_dat(M_settings = list(mu_supplied = ~ I(0))), "positive finite value")
  expect_error(make_test_dat(years = c(2000, 2001.5)), "integer values")
  expect_error(make_test_dat(proj_settings = list(n_proj = 1, n_mean = 1, F_mult = -1)),
    "non-negative multipliers")
  for (nm in names(obs)) {
    bad <- cod_obs
    bad[[nm]]$obs[1] <- -1
    expect_error(check_obs(bad), "non-negative values")
  }
})
