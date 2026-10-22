test_that("deterministic recruitment follows its formula with free boundaries", {
  obs <- cod_obs
  obs$pop_covariates$temp <- ifelse(obs$pop_covariates$age==2,sin(obs$pop_covariates$year),NA)
  d <- make_test_dat(data=obs,years=1983:1995,ages=2:8,
    rec_settings=list(process="off",mu_form=~temp),proj_settings=list(n_proj=2))
  p <- make_par(d);p$log_r0<-log(100);p$rec_beta[]<-c(2,.2)
  expect_false(any(c("log_r","log_sd_r","logit_phi_r") %in% names(p)))
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  r <- o$report()
  expect_equal(unname(log(r$recruitment)[1]),p$log_r0)
  expect_equal(unname(log(r$recruitment)[-1]),drop(d$rec$matrix[-1,]%*%p$rec_beta),ignore_attr=TRUE)
  expect_null(r$rec_residual)
  expect_equal(.rec_nll(log(r$recruitment),r$rec_log_mean,p,d),0)
  expect_true(all(is.finite(o$gr(o$par))))
  set.seed(12);sim<-nll_fun(p,d,simulate=TRUE)
  expect_equal(unname(sim$log_r),unname(log(r$recruitment)[-1]))
})

test_that("deterministic stock-recruit models use the same chronological parents", {
  for(curve in c("bh","ricker")) for(lag in c(0,2,4)) {
    obs<-cod_obs;obs$maturity$obs[obs$maturity$age==2]<-0
    d<-make_test_dat(data=obs,years=1983:1998,ages=2:8,
      rec_settings=list(process="off",mu_form=as.formula(paste0("~",curve,"(ssb,lag=",lag,")"))))
    p<-make_par(d)
    p$log_r0<-1
    if(!is.null(p$log_r_init))p$log_r_init[]<-seq_along(p$log_r_init)/10
    p<-.initialize_rec_curve(p,d)
    o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
    r<-o$report();i<-d$rec$eligible
    expect_equal(unname(log(r$recruitment)[d$rec$boundary]),unname(c(p$log_r0,p$log_r_init)))
    expect_equal(unname(log(r$recruitment)[i]),unname(.rec_log_curve(log(r$ssb[i-lag]),p,curve)),tolerance=1e-10)
    expect_true(all(is.finite(o$gr(o$par))))
    expect_equal(unname(.rec_log_curve(log(sum(.population_states(p,d,r$Z)$ssb_mat[i[1]-lag,])),p,curve)),unname(p$log_r0),tolerance=1e-8)
    # Other random processes may change parents; the deterministic relationship still holds.
    sim<-nll_fun(p,d,simulate=TRUE)
    ps<-p
    for(nm in intersect(c("log_n","log_f","log_m"),names(sim)))ps[[nm]]<-sim[[nm]]
    simulated<-RTMB::MakeADFun(function(p)nll_fun(p,d),ps,silent=TRUE)$report()
    expect_equal(unname(sim$log_r),unname(log(simulated$recruitment)[i]),tolerance=1e-10)
  }
})

test_that("deterministic recruitment rejects unused settings and supports saturated fixed means", {
  expect_error(make_test_dat(rec_settings=list(process="off")),"mu_form is required")
  for(nm in c("sd","phi")) {
    settings<-list(process="off",mu_form=~1);settings[[nm]]<-.2
    expect_error(make_test_dat(rec_settings=settings),"no residual process")
  }
  d<-make_test_dat(years=1983:1988,ages=2:6,rec_settings=list(process="off",mu_form=~factor(year)))
  expect_equal(ncol(d$rec$matrix),5)
  expect_error(make_test_dat(rec_settings=list(process="off",mu_form=~year+I(year*2))),"redundant")
  expect_error(make_test_dat(rec_settings=list(process="off",mu_form=~bh(ssb,lag=100))),"historical")
})

test_that("deterministic fits retain uncertainty, updates and evaluation workflows", {
  g<-expand.grid(year=2000:2024,age=1:5)
  data<-list(catch=transform(g,obs=100,unit="thousand fish"),
    index=transform(g,obs=100,survey="S",samp_time=.3,unit="thousand fish"),
    weight=transform(g,obs=age,unit="kg"),maturity=transform(g,obs=as.numeric(age>=3)))
  settings<-list(N_settings=list(process="off"),rec_settings=list(process="off",mu_form=~1),
    F_settings=list(process="off",mu_form=~1),index_settings=list(q_form=~1))
  d<-do.call(prepare_tam,c(list(data=data),settings));p<-make_par(d)
  p$log_r0<-log(800);p$rec_beta[]<-log(1000);p$log_mu_f[]<-log(.4)
  p$log_q[]<-log(.5);p$log_sd_catch[]<-log(.15);p$log_sd_index[]<-log(.2)
  set.seed(144);s<-nll_fun(p,d,simulate=TRUE)
  data<-d$obs
  for(type in c("catch","index"))data[[type]]$obs<-exp(s$log_obs[d$obs_map$type==type])
  f<-do.call(fit_tam,c(list(data=data,start_par=p,silent=TRUE),settings))
  expect_true(f$is_converged);expect_true(f$sdrep$pdHess)
  expect_length(f$obj$env$.random,0)
  expect_null(tidy_recruitment(f)$residual)
  expect_true(all(is.finite(f$pop$recruitment$se)))
  expect_equal(attr(logLik(f),"df"),length(f$opt$par))
  expect_equal(cAIC(f,details=TRUE)$edf_rec,0)
  shorter<-suppressMessages(update(f,years=2000:2023,start_par=f,silent=TRUE))
  expect_true(shorter$is_converged)
  expect_equal(update(f,rec_settings=list(process="iid"),evaluate=FALSE)$rec_settings$process,"iid")
  stochastic<-prepare_tam(data,N_settings=list(process="off"),rec_settings=list(process="iid",mu_form=~1),
    F_settings=list(process="off",mu_form=~1),index_settings=list(q_form=~1))
  start<-.merge_start_par(make_par(stochastic),.start_parameters(f,stochastic))
  expect_equal(unname(start$log_r),unname(log(f$rep$recruitment[-1])))
  future::plan(future::sequential)
  h<-fit_hindcast(f,folds=1,start_from_fit=TRUE,progress=FALSE)
  expect_true(all(h$fold_diagnostics$is_converged))
  expect_true(nrow(h$scores)>0)
  sim<-sim_tam(f,n=1,redraw_random=TRUE,par_uncertainty="none")
  expect_equal(sim$recruitment$est,f$pop$recruitment$est)
  set.seed(68);fixed<-sim_tam(f,n=1,redraw_random=TRUE,par_uncertainty="fixed")
  set.seed(68);joint<-sim_tam(f,n=1,redraw_random=TRUE,par_uncertainty="joint")
  expect_identical(joint,fixed)
  expect_true(all(is.finite(joint$recruitment$est)))
})
