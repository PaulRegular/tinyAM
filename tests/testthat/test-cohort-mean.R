cohort_data <- function() {
 d <- cod_obs
 d$pop_covariates$x <- sin(d$pop_covariates$year / 3) + d$pop_covariates$age / 20
 d$pop_covariates$x[d$pop_covariates$age == 2 | d$pop_covariates$year == 1983] <- NA_real_
 d
}

test_that('cohort covariates multiply survivors after plus-group combination', {
 for (process in c('off','iid','ar1','rw')) {
  d <- prepare_tam(cohort_data(), years=1983:1990,ages=2:6,
   N_settings=list(process=process,mu_form=~0+x),F_settings=list(process='iid',mu_form=~1))
  p<-make_par(d);p$mu_n[]<-.3
  z<-matrix(.5,length(d$years),length(d$ages))
  v<-.population_states(p,d,z)
  expected<-v$log_N[-nrow(z),-ncol(z)]-z[-nrow(z),-ncol(z)]
  expected[,ncol(expected)]<-log(exp(expected[,ncol(expected)])+exp(v$log_N[-nrow(z),ncol(z)]-z[-nrow(z),ncol(z)]))
  x<-matrix(d$N_modmat%*%p$mu_n,nrow(z))[-1,-1]
  expect_equal(v$pred_log_N[-1,-1],expected+x,ignore_attr=TRUE)
  expect_equal(v$log_N[1,],c(0,-.5,-1,-1.5,-2),ignore_attr=TRUE)
  expect_equal(v$log_N[,1],rep(0,8),ignore_attr=TRUE)
  if (process=='off') expect_equal(v$log_N[-1,-1],expected+x,ignore_attr=TRUE)
  o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  expect_true(all(is.finite(o$gr(o$par))))
  expect_equal(.check_tam_structure(d)$status[.check_tam_structure(d)$check=='N_modmat'],'pass')
 }
})

test_that('cohort formulas reject unsupported and redundant designs', {
 for (formula in list(~iid(year),~rw(year),~logistic(age),~(1|age),~offset(year))) {
  expect_error(prepare_tam(cod_obs,years=1983:1990,ages=2:6,N_settings=list(mu_form=formula)), 'ordinary fixed')
 }
 expect_error(prepare_tam(cod_obs,years=1983:1990,ages=2:6,N_settings=list(mu_form=~missing_x)), 'pop_covariates')
 expect_error(prepare_tam(cod_obs,years=1983:1990,ages=2:6,N_settings=list(process='rw',mu_form=~1)), 'cancel')
 expect_error(prepare_tam(cod_obs,years=1983:1990,ages=2:6,N_settings=list(process='rw',mu_form=~0+age)), 'cancel')
 expect_error(prepare_tam(cohort_data(),years=1983:1990,ages=2:6,N_settings=list(mu_form=~x+I(2*x))), 'redundant')
})

test_that('cohort projection designs hold terminal values without changing historical bases', {
 d<-cohort_data()
 a<-prepare_tam(d,years=1983:1990,ages=2:6,N_settings=list(mu_form=~poly(x,2)))
 b<-prepare_tam(d,years=1983:1990,ages=2:6,N_settings=list(mu_form=~poly(x,2)),proj_settings=list(n_proj=2))
 expect_equal(b$N_modmat[b$pop_data$year<=1990,],a$N_modmat)
 for(age in 3:6) expect_equal(b$N_modmat[b$pop_data$age==age & b$pop_data$year>1990,], b$N_modmat[rep(which(b$pop_data$age==age & b$pop_data$year==1990),2),], ignore_attr=TRUE)
 p<-make_par(b);p$mu_n[]<-c(.1,.2,.3)
 set.seed(4); s<-.population_states(p,b,matrix(.5,10,5),simulate=TRUE)
 expect_true(all(is.finite(s$log_N)))
 expect_true(all(is.finite(s$ssb_mat)))
 expect_equal(s$pred_log_N[1,2],-.5)
})

test_that('zero cohort coefficients preserve the existing likelihood and gradients', {
 for (process in c('off','iid','ar1','rw')) {
  a<-prepare_tam(cohort_data(),years=1983:1990,ages=2:6,N_settings=list(process=process),F_settings=list(process='iid',mu_form=~1))
  b<-prepare_tam(cohort_data(),years=1983:1990,ages=2:6,N_settings=list(process=process,mu_form=~0+x),F_settings=list(process='iid',mu_form=~1))
  pa<-make_par(a);pb<-make_par(b)
  oa<-RTMB::MakeADFun(function(p)nll_fun(p,a),pa,silent=TRUE)
  ob<-RTMB::MakeADFun(function(p)nll_fun(p,b),pb,silent=TRUE)
  expect_identical(pb[names(pa)],pa)
  expect_equal(ob$fn(ob$par),oa$fn(oa$par),tolerance=1e-10)
  expect_equal(as.numeric(ob$gr(ob$par)[names(ob$par)!='mu_n']),as.numeric(oa$gr(oa$par)),tolerance=1e-10)
 }
})
