catch_only_data <- function(fleets = FALSE) {
  x <- cod_obs
  x$index <- data.frame(year = integer(), age = integer(), obs = numeric(),
    survey = character(), samp_time = numeric())
  if (fleets) x$catch <- rbind(transform(x$catch, fleet = "Trawl", obs = .6 * obs),
    transform(x$catch, fleet = "Gillnet", obs = .4 * obs))
  x
}

catch_only_fixture <- function(proj_settings = NULL, fleets = FALSE, index_settings = list()) {
  d <- prepare_tam(catch_only_data(fleets), years = 1983:1990, ages = 2:5,
    proj_settings = proj_settings, index_settings = index_settings)
  p <- make_par(d)
  o <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  value <- as.numeric(o$fn(o$par))
  structure(list(dat = d, obj = o, rep = o$report(o$par), parameter_values = o$par,
    parameter_map = list(), opt = list(par = o$par, objective = value), sdrep = NULL), class = "tam_fit")
}

test_that("catch-only preparation creates no unused survey parameters", {
  x <- catch_only_data(); before <- x
  f <- catch_only_fixture()
  d <- f$dat
  expect_identical(nrow(d$obs$index), 0L)
  expect_equal(dim(d$q_modmat), c(0L, 0L))
  expect_equal(dim(d$sd_index_modmat), c(0L, 0L))
  expect_length(d$q_terms, 0L)
  expect_length(d$sd_index_terms, 0L)
  p <- make_par(d)
  expect_length(p$log_q, 0L)
  expect_length(p$log_sd_index, 0L)
  expect_false(any(grepl("(^log_q|^log_sd_index|^eta_q|^eta_sd_index)", names(f$obj$par))))
  expect_true(all(d$obs_map$type == "catch"))
  expect_equal(length(d$observed), sum(is.finite(d$obs$catch$obs) & d$obs$catch$obs > 0))
  expect_true(is.finite(f$opt$objective))
  expect_true(all(is.finite(f$obj$gr(f$obj$par))))
  explicit <- catch_only_fixture(index_settings = list(q_form = ~0, sd_form = ~0))
  expect_identical(f$obj$par, explicit$obj$par)
  expect_equal(f$opt$objective, explicit$opt$objective, tolerance = 1e-12)
  expect_equal(f$obj$gr(f$obj$par), explicit$obj$gr(explicit$obj$par), tolerance = 1e-12)
  # Unused survey formulas do not create coefficients, effects or correlations.
  scoped <- catch_only_fixture(index_settings = list(q_form = ~0 + for_survey("RV", ~s(age,k=4)),
    sd_form = ~iid(age), sd_supplied = ~I(.2)))
  expect_identical(f$obj$par, scoped$obj$par)
  expect_equal(f$opt$objective, scoped$opt$objective)
  expect_identical(x, before)
})

test_that("catch-only projections and fleet preparation retain an empty index", {
  for (fleets in c(FALSE, TRUE)) {
    f <- catch_only_fixture(list(n_proj = 2, n_mean = 2, F_mult = 1), fleets = fleets)
    expect_equal(tail(f$dat$years, 2L), c(1991, 1992))
    expect_identical(f$dat$obs$index$is_proj, logical())
    expect_true(is.finite(f$opt$objective))
    expect_true(all(is.finite(f$obj$gr(f$obj$par))))
    tabs <- tidy_obs_pred(f)
    expect_identical(nrow(tabs$index), 0L)
    expect_true(all(c("year", "age", "obs", "survey", "samp_time", "pred", "sd", "q",
      "q_lwr", "q_upr", "q_se", "q_se_scale", "q_unit", "std_res") %in% names(tabs$index)))
    expect_type(tabs$index$pred, "double")
    expect_type(tabs$index$std_res, "double")
    expect_type(tabs$index$q_se_scale, "character")
    f$obs_pred <- tabs
    out <- suppressMessages(tidy_tam(model_list = list(A = f, B = f), process_residuals = FALSE))
    expect_identical(nrow(out$obs_pred$index), 0L)
    expect_true(all(c("model", "pred", "obs") %in% names(out$obs_pred$index)))
    expect_true(nrow(out$obs_pred$catch) > 0L)
    expect_true(nrow(out$pop$ssb) > 0L)
    draws <- sim_tam(f, n = 2, par_uncertainty = "none", progress = FALSE, seed = 321)
    expect_identical(nrow(draws$index), 0L)
    expect_true(all(c("sim", "year", "age", "survey", "obs") %in% names(draws$index)))
    expect_type(draws$index$obs, "double")
  }
})

test_that("catch-only OSA tables retain numeric empty survey columns", {
  f <- catch_only_fixture(); f$obs_pred <- tidy_obs_pred(f)
  original <- f$parameter_values
  local_mocked_bindings(oneStepPredict = function(obj, method, subset = NULL, ...) {
    selected <- if (is.null(subset)) seq_along(f$dat$observed) else subset
    data.frame(residual = rep(.2, length(selected)))
  }, .package = "RTMB")
  osa <- add_osa_residuals(f, subset = 1:2, trace = FALSE)
  expect_identical(nrow(osa$obs_pred$index), 0L)
  expect_type(osa$obs_pred$index$osa_res, "double")
  expect_type(osa$obs_pred$index$osa_method, "character")
  expect_equal(sum(is.finite(osa$obs_pred$catch$osa_res)), 2L)
  extracted <- residuals(osa, type = "osa")
  expect_true(all(extracted$type == "catch"))
  expect_equal(sum(is.finite(extracted$residual)), 2L)
  expect_identical(f$parameter_values, original)
})

test_that("all-empty nested tables preserve union schema and identifier types", {
  x <- list(Official = list(index = data.frame(year = integer(), survey = character())),
    tinyAM = list(index = data.frame(year = integer(), survey = character(), pred = numeric(),
      osa_method = character())), Missing = list())
  before <- x
  for (type in c("auto", "numeric", "character", "factor")) {
    out <- stack_nested(x, label = "model", label_type = type)$index
    expect_identical(nrow(out), 0L)
    expect_identical(names(out), c("model", "year", "survey", "pred", "osa_method"))
    expect_type(out$year, "integer")
    expect_type(out$pred, "double")
    expect_type(out$osa_method, "character")
    if (type == "factor") expect_identical(levels(out$model), c("Official", "tinyAM"))
    if (type == "numeric") expect_type(out$model, "double")
    if (type == "character") expect_type(out$model, "character")
  }
  expect_false("model" %in% names(stack_nested(x, label = NULL)$index))
  expect_identical(x, before)
  partial <- x
  partial$tinyAM$index <- data.frame(year = 2000L, survey = "RV", pred = 2, osa_method = "oneStepGeneric")
  stacked <- stack_nested(partial, label_type = "factor")$index
  expect_identical(nrow(stacked), 1L)
  expect_equal(stacked$pred, 2)
})
