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
  expect_true(all(out$residual[out$type=="catch" & is.finite(out$residual)]==.7))
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
