test_that("F means use historical support and terminal projection rows", {
  for (form in list(~factor(year) + factor(age), ~poly(year, 2) + factor(age), ~year:factor(age))) {
    historical <- prepare_tam(cod_obs, years=1983:1990, ages=2:6,
      F_settings=list(process="off",mu_form=form))
    projected <- prepare_tam(cod_obs, years=1983:1990, ages=2:6,
      F_settings=list(process="off",mu_form=form), proj_settings=list(n_proj=2))
    past <- projected$pop_data$year<=1990
    expect_equal(projected$F_modmat[past,],historical$F_modmat,ignore_attr=TRUE)
    for (age in 2:6) expect_equal(projected$F_modmat[!past & projected$pop_data$age==age,,drop=FALSE],
      historical$F_modmat[rep(which(historical$pop_data$year==1990 & historical$pop_data$age==age),2),,drop=FALSE],ignore_attr=TRUE)
    expect_equal(qr(projected$F_modmat[past,])$rank,ncol(projected$F_modmat))
  }
})

test_that("unused F factor levels are dropped but genuine redundancy is rejected", {
  data <- cod_obs
  data$pop_covariates$period <- factor(ifelse(data$pop_covariates$year<1990,"early","late"), levels=c("late","early","unused"))
  d <- prepare_tam(data,years=1983:1990,ages=2:6,F_settings=list(process="off",mu_form=~period+factor(age)))
  expect_false(any(grepl("unused",colnames(d$F_modmat))))
  expect_equal(d$F_modmat[,"periodearly"],as.numeric(d$pop_data$period=="early"),ignore_attr=TRUE)
  expect_error(prepare_tam(data,years=1983:1990,ages=2:6,
    F_settings=list(process="off",mu_form=~year+I(2*year))),"rank-deficient")
})

test_that("warm mean coefficients align for peels but not changed reference coding", {
  form <- ~factor(year)+factor(age)
  d <- prepare_tam(cod_obs,years=1983:1990,ages=2:6,F_settings=list(process="off",mu_form=form))
  p <- make_par(d); p$log_mu_f[] <- seq_along(p$log_mu_f)/10
  fake <- structure(list(dat=d,sdrep=list()),class="tam_fit")
  fake$sdrep <- structure(list(),class="sdreport")
  fake$parameter_values <- unlist(p)
  fake$obj <- list(env=list(parList=function(par)p))
  peel <- prepare_tam(cod_obs,years=1983:1989,ages=2:6,F_settings=list(process="off",mu_form=form))
  start <- .merge_start_par(make_par(peel),.start_parameters(fake,peel))
  expect_equal(start$log_mu_f,p$log_mu_f[names(start$log_mu_f)])
  changed <- prepare_tam(cod_obs,years=1984:1990,ages=2:6,F_settings=list(process="off",mu_form=form))
  expect_message(starts <- .start_parameters(fake,changed),"basis changed")
  expect_null(starts$log_mu_f)
})

test_that("historical mean designs retain explicit contrast reference levels", {
  data<-cod_obs
  data$pop_covariates$phase<-factor(ifelse(data$pop_covariates$year<1987,"a","b"),levels=c("a","b","unused"))
  contrasts(data$pop_covariates$phase)<-contr.treatment(c("a","b","unused"),base=2)
  d<-prepare_tam(data,years=1983:1990,ages=2:6,F_settings=list(process="off",mu_form=~phase),proj_settings=list(n_proj=1))
  expect_named(as.data.frame(d$F_modmat),c("(Intercept)","phasea"))
  expect_equal(unname(d$F_modmat[,"phasea"]),as.numeric(d$pop_data$year<1987))
})
