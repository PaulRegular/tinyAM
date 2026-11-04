osa_fixture <- function() {
  d<-prepare_tam(cod_obs,years=1983:1987,ages=2:4,
    N_settings=list(process="off"),rec_settings=list(process="iid",mu_form=~1),
    F_settings=list(process="iid",mu_form=~1),index_settings=list(q_form=~1))
  p<-make_par(d);p$log_r0<-log(100);p$rec_beta[]<-log(100)
  p$log_r[]<-log(100);p$log_mu_f[]<-log(.2);p$log_f[]<-log(.2)
  o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,random=c("log_r","log_f"),silent=TRUE)
  o$fn(o$par)
  f<-structure(list(dat=d,obj=o,rep=o$report(),parameter_values=o$env$last.par,
    parameter_map=list(),sdrep=NULL),class="tam_fit")
  f$obs_pred<-tidy_obs_pred(f)
  f
}

test_that("OSA requests reuse caches, align rows and never refit", {
  f<-osa_fixture();before<-f$parameter_values
  mode<-f$obj$env$last.par;best<-f$obj$env$last.par.best
  calls<-0L
  local_mocked_bindings(oneStepPredict=function(obj,method,trace,...) {
    calls<<-calls+1L
    expect_false(identical(obj,f$obj))
    expect_equal(obj$env$last.par,f$parameter_values)
    expect_identical(method,"oneStepGaussianOffMode")
    expect_false(trace)
    data.frame(residual=seq_len(sum(f$dat$is_observed))/10)
  },.package="RTMB")
  local_mocked_bindings(fit_tam=function(...)stop("Must not refit"))
  f$obs_pred$catch$osa_res<-.7
  f$obs_pred$index<-f$obs_pred$index[rev(seq_len(nrow(f$obs_pred$index))),]
  out<-residuals(f,type="osa",trace=FALSE)
  expect_equal(calls,1L)
  catch_reference <- split(seq_len(sum(f$dat$is_observed))/10,
    f$dat$obs_map$type[f$dat$is_observed])$catch
  expect_equal(out$residual[out$type=="catch" & is.finite(out$residual)], catch_reference)
  reference<-split(seq_len(sum(f$dat$is_observed))/10,f$dat$obs_map$type[f$dat$is_observed])$index
  rows<-f$dat$obs$index[f$dat$obs$index$obs>0 & !is.na(f$dat$obs$index$obs),]
  expected<-reference[match(paste(out$year,out$age,out$survey),paste(rows$year,rows$age,rows$survey))]
  valid<-out$type=="index" & is.finite(out$residual)
  expect_equal(out$residual[valid],expected[valid])
  expect_identical(f$parameter_values,before);expect_identical(f$obj$env$last.par,mode)
  expect_identical(f$obj$env$last.par.best,best)
  expect_false("osa_res" %in% names(f$obs_pred$index))
  f$obs_pred$index$osa_res<-.5
  residuals(f,type="osa")
  expect_equal(calls,1L)
  expect_error(residuals(f,trace=FALSE),"only used")
})

test_that("actual OSA computation preserves objective values and missing masks", {
  f<-osa_fixture();mode<-f$obj$env$last.par;before<-f$parameter_values
  outer<-f$obj$par;objective<-f$obj$fn(outer)
  stored<-f$obj$env$last.par
  x<-residuals(f,type="osa",trace=FALSE,parallel=FALSE)
  expect_equal(nrow(x),sum(vapply(f$obs_pred,nrow,integer(1))))
  valid<-is.finite(x$obs) & x$obs>0 & !x$is_proj
  expect_true(all(is.finite(x$residual[valid])))
  expect_true(all(is.na(x$residual[!valid])))
  expect_identical(f$parameter_values,before)
  expect_identical(f$obj$env$last.par,stored)
  expect_equal(f$obj$fn(outer),objective,tolerance=1e-10)
  expect_false("osa_res" %in% names(f$obs_pred$catch))
  local_mocked_bindings(oneStepPredict=function(...)stop("OSA probe failed"),.package="RTMB")
  expect_error(residuals(f,type="osa"),"Could not calculate OSA residuals.*OSA probe failed")
  expect_identical(f$parameter_values,before)
})

test_that("partial OSA results retain requested order and original row identities", {
  f <- osa_fixture()
  n <- length(f$dat$observed)
  selected <- c(n, 1L)
  local_mocked_bindings(oneStepPredict = function(obj, subset, ...) {
    expect_identical(subset, selected)
    data.frame(residual = c(.4, -.7))
  }, .package = "RTMB")
  tables <- tidy_obs_pred(f, add_osa_res = TRUE, subset = selected, trace = FALSE)
  actual <- c(tables$catch$osa_res, tables$index$osa_res)[f$dat$is_observed]
  expected <- rep(NA_real_, n); expected[selected] <- c(.4, -.7)
  expect_equal(actual, expected)
  local_mocked_bindings(oneStepPredict = function(obj, conditional, ...) {
    data.frame(residual = rep(.3, n - length(conditional)))
  }, .package = "RTMB")
  tables <- tidy_obs_pred(f, add_osa_res = TRUE, conditional = seq_len(n - 1L), trace = FALSE)
  actual <- c(tables$catch$osa_res, tables$index$osa_res)[f$dat$is_observed]
  expect_true(all(is.na(actual[-n])))
  expect_equal(actual[n], .3)
})

test_that("OSA methods and numerical controls are forwarded once", {
  f <- osa_fixture()
  calls <- list()
  local_mocked_bindings(oneStepPredict = function(obj, method, ...) {
    controls <- list(...)
    expect_false(controls$discrete)
    expect_false("reverse" %in% names(controls))
    calls[[length(calls) + 1L]] <<- list(method = method, controls = controls)
    data.frame(residual = rep(.2, length(f$dat$observed)))
  }, .package = "RTMB")
  methods <- c("oneStepGaussianOffMode", "oneStepGaussian", "oneStepGeneric", "fullGaussian", "cdf")
  for (method in methods) {
    out <- residuals(f, type = "osa", method = method, trace = FALSE)
    expect_identical(unique(out$osa_method), method)
    expect_identical(attr(out, "osa")$method, method)
    expect_identical(attr(out, "osa")$settings$method, method)
  }
  expect_identical(vapply(calls, `[[`, character(1), "method"), methods)
  out <- tidy_obs_pred(f, add_osa_res = TRUE, method = "oneStepGeneric", trace = FALSE)
  expect_identical(unique(out$catch$osa_method), "oneStepGeneric")
  expect_identical(attr(out, "osa")$method, "oneStepGeneric")
  expect_error(residuals(f, type = "osa", methods = "oneStepGeneric"), "Use.*method.*not.*methods")
  expect_error(tidy_obs_pred(f, add_osa_res = TRUE, methods = "oneStepGeneric"), "Use.*method.*not.*methods")
  expect_error(add_osa_residuals(f, methods = "oneStepGeneric"), "Use.*method.*not.*methods")
  expect_error(residuals(f, type = "osa", method = "OneStepGeneric"), "method.*must be one of")
  expect_error(residuals(f, method = "cdf"), "only used with")
})

test_that("cached OSA extraction is cheap but explicit requests refresh both tables", {
  f <- osa_fixture()
  f$obs_pred$catch$osa_res <- .7
  f$obs_pred$index$osa_res <- .5
  calls <- 0L
  local_mocked_bindings(oneStepPredict = function(obj, method, ...) {
    calls <<- calls + 1L
    data.frame(residual = rep(-.2, length(f$dat$observed)))
  }, .package = "RTMB")
  cached <- residuals(f, type = "osa")
  expect_identical(calls, 0L)
  expect_true(all(is.na(cached$osa_method)))
  expect_true(all(cached$residual[cached$type == "catch" & is.finite(cached$residual)] == .7))
  refreshed <- residuals(f, type = "osa", method = "oneStepGaussianOffMode")
  expect_identical(calls, 1L)
  expect_true(all(refreshed$residual[is.finite(refreshed$residual)] == -.2))
  refreshed <- residuals(f, type = "osa", trace = FALSE)
  expect_identical(calls, 2L)
  expect_true(all(refreshed$residual[is.finite(refreshed$residual)] == -.2))
  expect_true(all(f$obs_pred$catch$osa_res == .7))
  expect_true(all(f$obs_pred$index$osa_res == .5))
})

test_that("explicit upstream numerical controls are retained", {
  f <- osa_fixture()
  local_mocked_bindings(oneStepPredict = function(obj, method, discrete, reverse, range, splineApprox, ...) {
    expect_identical(method, "oneStepGeneric")
    expect_true(discrete)
    expect_true(reverse)
    expect_identical(range, c(-10, 10))
    expect_false(splineApprox)
    data.frame(residual = rep(.2, length(f$dat$observed)))
  }, .package = "RTMB")
  out <- add_osa_residuals(f, method = "oneStepGeneric", discrete = TRUE,
    reverse = TRUE, range = c(-10, 10), splineApprox = FALSE)
  expect_identical(out$osa$requested_settings$reverse, TRUE)
  expect_identical(out$osa$settings$discrete, TRUE)
  expect_identical(out$osa$settings$range, c(-10, 10))
})

test_that("retaining OSA results preserves custom tables, intervals and fitted state", {
  f <- osa_fixture()
  f$obs_pred$index <- f$obs_pred$index[rev(seq_len(nrow(f$obs_pred$index))), ]
  f$obs_pred$index$custom_note <- paste0("row", seq_len(nrow(f$obs_pred$index)))
  f$obs_pred$index$q_lwr <- seq_len(nrow(f$obs_pred$index)) / 100
  f$obs_pred$index$q_upr <- f$obs_pred$index$q_lwr + 1
  attr(f$obs_pred, "interval") <- .8
  attr(f$obs_pred$index, "custom_metadata") <- list(source = "retained")
  f$sdrep <- list(untouched = "uncertainty")
  original <- f$obs_pred
  parameters <- f$parameter_values
  mode <- f$obj$env$last.par
  best <- f$obj$env$last.par.best
  local_mocked_bindings(oneStepPredict = function(obj, ...) {
    data.frame(residual = seq_along(f$dat$observed) / 10)
  }, .package = "RTMB")
  local_mocked_bindings(fit_tam = function(...) stop("Must not refit"))
  out <- add_osa_residuals(f, method = "oneStepGeneric", trace = FALSE)
  expect_s3_class(out, "tam_fit")
  for (name in names(original)) {
    for (column in names(original[[name]])) expect_identical(out$obs_pred[[name]][[column]], original[[name]][[column]])
  }
  expect_identical(attr(out$obs_pred, "interval"), .8)
  expect_identical(attr(out$obs_pred$index, "custom_metadata"), list(source = "retained"))
  expect_identical(f$obs_pred, original)
  expect_identical(f$parameter_values, parameters)
  expect_identical(f$obj$env$last.par, mode)
  expect_identical(f$obj$env$last.par.best, best)
  expect_identical(out$sdrep, f$sdrep)
  expect_identical(out$obj, f$obj)
  expect_identical(out$osa$method, "oneStepGeneric")
  expect_identical(attr(out$obs_pred, "osa"), out$osa)
  expect_identical(attr(out$obs_pred$index, "osa"), out$osa)
  expect_identical(out$osa$n_requested, length(f$dat$observed))
  expect_identical(out$osa$n_failed, 0L)
  raw <- out$osa$coverage
  rows <- match(paste(out$obs_pred$index$year, out$obs_pred$index$age, out$obs_pred$index$survey),
    paste(raw$year[raw$type == "index"], raw$age[raw$type == "index"], raw$survey[raw$type == "index"]))
  expect_equal(out$obs_pred$index$osa_res,
    raw$observation[raw$type == "index"][rows] / 10)
  local_mocked_bindings(oneStepPredict = function(...) stop("OSA failure"), .package = "RTMB")
  expect_error(add_osa_residuals(f, method = "cdf"), "Could not calculate OSA residuals.*OSA failure")
  expect_identical(f$obs_pred, original)
  expect_false("osa" %in% names(f))
})

test_that("OSA coverage distinguishes conditioning, unrequested and failed observations", {
  f <- osa_fixture()
  selected <- c(3L, 1L, 4L)
  local_mocked_bindings(oneStepPredict = function(obj, subset, conditional, ...) {
    expect_identical(subset, selected)
    expect_identical(conditional, 2L)
    data.frame(residual = c(.1, NA_real_, Inf))
  }, .package = "RTMB")
  expect_warning(out <- add_osa_residuals(f, subset = selected, conditional = 2L, trace = FALSE),
    "unavailable residuals for 2 requested observations")
  expect_identical(out$osa$n_requested, 3L)
  expect_identical(out$osa$n_available, 1L)
  expect_identical(out$osa$n_failed, 2L)
  expect_identical(out$osa$coverage$status[1:5], c("failed", "conditional", "available", "failed", "not_requested"))
  extracted <- residuals(out, type = "osa")
  expect_identical(sum(is.finite(extracted$residual)), 1L)
  expect_equal(extracted$residual[is.finite(extracted$residual)], .1)
  expect_error(add_osa_residuals(f, subset = c(1, 1)), "distinct observation indices")
  expect_error(add_osa_residuals(f, subset = -1), "distinct observation indices")
  expect_error(add_osa_residuals(f, subset = 1e30), "distinct observation indices")
  expect_error(add_osa_residuals(f, subset = 1, conditional = 1), "must not overlap")
})

test_that("keyed OSA attachment distinguishes surveys and fleets", {
  source <- list(
    catch = data.frame(year = 2000, age = 3, fleet = c("A", "B"), osa_res = c(.1, .2), osa_method = "cdf"),
    index = data.frame(year = 2000, age = 3, survey = c("A", "B"), osa_res = c(.3, .4), osa_method = "cdf")
  )
  target <- lapply(source, function(d) d[2:1, setdiff(names(d), c("osa_res", "osa_method")), drop = FALSE])
  target$catch$unit <- "thousand fish"
  out <- .attach_osa_tables(target, source, list(method = "cdf"))
  expect_equal(out$catch$osa_res, c(.2, .1))
  expect_equal(out$index$osa_res, c(.4, .3))
  expect_identical(out$catch$unit, target$catch$unit)
  target$catch$fleet <- NULL
  expect_error(.attach_osa_tables(target, source), "must retain.*keys")
})

test_that("OSA retention rebuilds missing tables and handles empty requests", {
  f <- osa_fixture()
  f$obs_pred <- NULL
  f$pop <- structure(list(), interval = .8)
  local_mocked_bindings(oneStepPredict = function(...) stop("No observations requested"), .package = "RTMB")
  out <- add_osa_residuals(f, subset = integer())
  expect_identical(attr(out$obs_pred, "interval"), .8)
  expect_true(all(is.na(out$obs_pred$catch$osa_res)))
  expect_identical(out$osa$n_requested, 0L)
  expect_true(all(out$osa$coverage$status == "not_requested"))
  expect_null(f$obs_pred)
})

test_that("oneStepGeneric is available without changing fitted parameters", {
  f <- osa_fixture()
  before <- f$parameter_values
  mode <- f$obj$env$last.par
  out <- residuals(f, type = "osa", method = "oneStepGeneric",
    subset = seq_len(2L), trace = FALSE, parallel = FALSE)
  expect_equal(sum(is.finite(out$residual)), 2L)
  expect_identical(unique(out$osa_method), "oneStepGeneric")
  expect_identical(f$parameter_values, before)
  expect_identical(f$obj$env$last.par, mode)
})
