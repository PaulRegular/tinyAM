test_that("SAM observations retain inputs, keys, timing and missing values", {
  f <- sam_fixture()
  before <- unserialize(serialize(f, NULL))
  obs <- sam_to_tam_obs(f)
  expect_true(check_obs(obs))
  expect_equal(nrow(obs$catch), 9)
  expect_true(all(is.na(obs$catch$obs[obs$catch$year == 2002])))
  expect_equal(obs$catch$obs[obs$catch$year == 2001 & obs$catch$age == 1], 0)
  expect_true(all(obs$index$samp_time[obs$index$fleet_id == 2] == .125))
  expect_true(all(is.na(obs$index$effort)))
  expect_equal(obs$weight$M_assumption[obs$weight$age == 2], rep(.3, 3))
  expect_equal(obs$index$q_key[obs$index$age == 1], rep(0, 5))
  expect_identical(f, before)
})

test_that("SAM settings preserve sharing and support explicit nested overrides", {
  f <- sam_fixture()
  s <- sam_to_tam_settings(f)
  d <- do.call(make_dat, c(list(obs = sam_to_tam_obs(f)), s))
  expect_equal(s$N_settings, list(process = "iid", init = "free"))
  expect_null(s$proj_settings)
  expect_equal(ncol(d$q_modmat), 2)
  expect_equal(ncol(d$sd_catch_modmat), 2)
  expect_equal(ncol(d$sd_index_modmat), 3)
  expect_true(is.finite(nll_fun(make_par(d), d)))
  o <- sam_to_tam_settings(f, list(N_settings = list(init = "exp"), proj_settings = NULL))
  expect_equal(o$N_settings, list(process = "iid", init = "exp"))
  expect_null(o$proj_settings)
  expect_error(sam_to_tam_settings(f, list(bogus = 1)), "Unknown")
  expect_error(sam_to_tam_settings(f, list(obs = list())), "Unknown")
})

test_that("audit distinguishes exact mappings, approximations and unresolved inputs", {
  f <- sam_fixture()
  status <- function(a, field) a$tam_status[match(field, a$sam_setting)]
  a <- sam_to_tam_audit(f)
  expect_equal(status(a, "keyLogFpar"), "supported")
  expect_equal(status(a, "keyVarObs"), "supported")
  expect_equal(status(a, "initState"), "partially_supported")
  expect_equal(status(a, "fbarRange"), "partially_supported")
  f$conf$corFlag <- 2
  expect_equal(status(sam_to_tam_audit(f), "corFlag"), "unsupported")
  s <- sam_to_tam_settings(f, list(index_settings = list(q_form = ~1), N_settings = list(process = "off")))
  a <- sam_to_tam_audit(f, s)
  expect_equal(status(a, "keyLogFpar"), "unsupported")
  expect_equal(status(a, "keyVarLogN"), "unsupported")
  f$conf$initState <- NULL
  f$conf$newFeature <- 1
  a <- sam_to_tam_audit(f)
  expect_equal(status(a, "initState"), "not_checked")
  expect_equal(status(a, "newFeature"), "not_checked")
  f$data$propMat[1, ] <- NA
  expect_true(all(is.na(sam_to_tam_obs(f)$maturity$obs[sam_to_tam_obs(f)$maturity$year == 2000])))
  expect_equal(status(sam_to_tam_audit(f), "selected years"), "unsupported")
  s <- sam_to_tam_settings(f, list(years = 2001:2002))
  expect_equal(status(sam_to_tam_audit(f, s), "selected years"), "partially_supported")
})

test_that("fixed q cells and single SD groups use valid designs", {
  f <- sam_fixture()
  f$conf$keyLogFpar[2, 1] <- -1
  f$conf$keyVarObs[1, ] <- 0
  s <- sam_to_tam_settings(f)
  obs <- sam_to_tam_obs(f)
  d <- do.call(make_dat, c(list(obs = obs), s))
  expect_equal(ncol(d$sd_catch_modmat), 1)
  expect_true(all(d$q_modmat[d$obs$index$q_key == -1, ] == 0))
  f$conf$keyLogFpar[,] <- -1
  s <- sam_to_tam_settings(f)
  d <- do.call(make_dat, c(list(obs = sam_to_tam_obs(f)), s))
  expect_equal(ncol(d$q_modmat), 0)
  expect_true(is.finite(nll_fun(make_par(d), d)))
})

test_that("untranslatable observation structures are rejected", {
  expect_error(sam_to_tam_obs(list()), "fitted SAM")
  f <- sam_fixture(); f$opt$objective <- NA_real_
  expect_error(sam_to_tam_obs(f), "finite optimizer")
  f <- sam_fixture(); f$data$fleetTypes[2] <- 0
  expect_error(sam_to_tam_obs(f), "catch fleet")
  f <- sam_fixture(); f$data$fleetTypes[2] <- 3
  expect_error(sam_to_tam_obs(f), "fleet type|type-2|Type")
  f <- sam_fixture(); f$data$aux[2, ] <- f$data$aux[1, ]
  expect_error(sam_to_tam_obs(f), "Duplicate")
  f <- sam_fixture(); f$conf$maxAgePlusGroup[2] <- 1
  expect_error(sam_to_tam_obs(f), "Survey plus groups")
  f <- sam_fixture(); f$conf$maxAgePlusGroup[1] <- 0
  expect_error(sam_to_tam_obs(f), "terminal catch plus group")
})

test_that("SAM comparison reports native states and honest uncertainty", {
  skip_if_not_installed("stockassessment")
  f <- sam_fixture()
  before <- unserialize(serialize(f, NULL))
  x <- sam_to_tam_comparison(f)
  expect_s3_class(x, "tam_comparison")
  expect_false(inherits(x, "tam_fit"))
  expect_null(x$obj)
  expect_equal(x$pop$N$est, as.vector(stockassessment::ntable(f)))
  expect_equal(x$pop$F$est, as.vector(stockassessment::faytable(f)))
  expect_equal(x$pop$ssb$est, 1:3)
  expect_equal(x$pop$ssb$se, rep(.1, 3))
  expect_equal(x$pop$ssb$se_scale, rep("log", 3))
  expect_equal(x$pop$ssb$upr, (1:3) * exp(qnorm(.975) * .1))
  expect_true(all(is.na(x$pop$abundance$se)))
  expect_equal(x$obs_pred$index$sd[x$obs_pred$index$age == 2], rep(.4, 3))
  expect_equal(x$obs_pred$index$q[x$obs_pred$index$age == 2], rep(.6, 3))
  tabs <- tidy_tam(x, interval = .8)
  expect_equal(tabs$pop$ssb$upr, (1:3) * exp(qnorm(.9) * .1))
  expect_error(fit_tam(x), "Missing required tables")
  expect_error(sim_tam(x), "tam_fit")
  expect_error(fit_retro(x), "tam_fit")
  expect_error(update(x), "reporting")
  x$dat$years <- 2001:2002
  expect_equal(sort(unique(tidy_tam(x)$pop$N$year)), 2001:2002)
  expect_identical(f, before)
  f$sdrep <- NULL
  expect_null(sam_to_tam_comparison(f)$pop$ssb)
})
