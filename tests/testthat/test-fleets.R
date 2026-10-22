multi_data <- function(years=2000:2006,ages=1:6) {
 g<-expand.grid(year=years,age=ages)
 list(catch=rbind(transform(g,obs=100,fleet='A',unit='thousand fish'),transform(g,obs=50,fleet='B',unit='thousand fish')),
 index=transform(g,obs=100,survey='S',samp_time=.4,q_block=factor(age),unit='thousand fish'),
 weight=transform(g,obs=age^2/10,unit='kg'),maturity=transform(g,obs=.5*plogis(age-3)))
}

test_that('fleet catches use common total mortality and sum correctly', {
 obs<-multi_data();obs$weight$catch_weight<-2*obs$weight$obs
 original<-obs
 d<-prepare_tam(obs,ages=1:4,N_settings=list(process='off'),F_settings=list(process='off',mu_form=~fleet+age),index_settings=list(q_form=~1))
 p<-make_par(d);p$log_mu_f[]<-c(log(.2),log(2),.1)
 o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE);r<-o$report(o$par)
 expect_equal(r$F,r$F_fleet[,,1]+r$F_fleet[,,2],ignore_attr=TRUE)
 expect_equal(r$Z,r$F+.2,ignore_attr=TRUE)
 expect_true(all(is.finite(o$gr(o$par))))
 for (fleet in d$fleets) {
  i<-d$obs_map$type=='catch' & !is.na(d$obs_map$fleet) & d$obs_map$fleet==fleet
  expected<-r$N*r$F_fleet[,,fleet]/r$Z*(-expm1(-r$Z))
  expect_equal(exp(r$log_pred[i]),as.vector(expected),ignore_attr=TRUE)
 }
 expect_equal(r$total_catch_pred,rowSums(r$catch_fleet_pred),ignore_attr=TRUE)
 expect_equal(r$total_yield_pred,rowSums(r$yield_fleet_pred),ignore_attr=TRUE)
 expect_equal(r$W_catch[,4],rowSums(r$N_plus*d$W_catch_plus_input)/rowSums(r$N_plus),ignore_attr=TRUE)
 expect_equal(r$W_catch,2*r$W,ignore_attr=TRUE)
 expect_identical(obs,original)
 for(f in d$fleets)expect_equal(subset(d$obs$catch,fleet==f & age==4)$obs,rep(if(f=='A')300 else 150,7))
})

test_that('F process densities are independent across fleets, including age correlation', {
 d<-prepare_tam(multi_data(),ages=1:4,F_settings=list(process='cor_rw',sd_form=~fleet+age))
 p<-make_par(d);p$sd_beta_f[]<-c(log(.2),log(1.4),.05)
 set.seed(1);p$log_f[]<-rnorm(length(p$log_f));mean<-.F_matrix(0,d)
 for(rho in c(-.5,0,.6)) {
  p$atanh_rho_f<-atanh(rho)
  sd<-.process_sd(p,d,'F'); expected<-0
  for(f in 1:2) {
   cols<-(f-1)*4+1:4
   Sigma<-outer(sd[cols],sd[cols])*rho^abs(outer(1:4,1:4,'-'))
   increments<-diff(p$log_f[,cols,drop=FALSE])
   for(i in seq_len(nrow(increments)))expected<-expected-.5*(4*log(2*pi)+as.numeric(determinant(Sigma,logarithm=TRUE)$modulus)+sum(increments[i,]*solve(Sigma,increments[i,])))
  }
  expect_equal(.F_density(p$log_f,mean,p,d),expected,tolerance=1e-9)
  shifted<-sweep(p$log_f,2,seq_len(ncol(p$log_f)),'+')
  expect_equal(.F_density(shifted,mean,p,d),expected,tolerance=1e-9)
 }
})

test_that('fleet simulations share age correlation but no cross-fleet innovations', {
 T<-4000L; A<-3L
 d<-list(years=seq_len(T),ages=1:A,fleets=c('A','B'),is_proj=rep(FALSE,T),
  F_settings=list(process='cor_rw'),process_sd=list(F=list(default=FALSE,parameter='sd_beta_f',matrix=diag(2*A))))
 p<-list(log_f=matrix(0,T,2*A),sd_beta_f=log(c(.2,.3,.4,.5,.6,.7)),atanh_rho_f=atanh(.6))
 set.seed(29);x<-.F_states(p,d,matrix(0,T,2*A),simulate=TRUE)$states
 expect_equal(x[1,],rep(0,2*A))
 delta<-diff(x)
 expect_equal(diag(cov(delta)),exp(2*p$sd_beta_f),tolerance=.025)
 expect_equal(cor(delta)[1,2],.6,tolerance=.035)
 expect_lt(max(abs(cor(delta)[1:3,4:6])),.06)
})

test_that('fleet identities, weights and SD designs are validated', {
 obs<-multi_data()
 expect_true(check_tam_data(obs))
 bad<-obs;bad$catch$fleet[1]<-NA;expect_error(check_tam_data(bad),'fleet labels')
 bad<-obs;bad$catch<-bad$catch[-1,];expect_error(check_tam_data(bad),'coverage is incomplete')
 bad<-obs;bad$weight$catch_weight<-NA_real_;expect_error(check_tam_data(bad),'catch_weight')
 d<-prepare_tam(obs,F_settings=list(process='iid',mu_form=~fleet+factor(age),sd_form=~fleet))
 expect_equal(nrow(d$process_sd$F$matrix),12L)
 p<-make_par(d);expect_equal(dim(p$log_f),c(7L,12L))
 expect_setequal(d$F_columns$fleet,c('A','B'))
 expect_true(all(.check_tam_structure(d)$status=='pass'))
})

test_that('catch weights alter yield but not historical population likelihoods', {
 obs<-multi_data()
 a<-prepare_tam(obs,ages=1:4)
 obs$weight$catch_weight<-obs$weight$obs*3
 b<-prepare_tam(obs,ages=1:4)
 pa<-make_par(a);pb<-make_par(b)
 oa<-RTMB::MakeADFun(function(p)nll_fun(p,a),pa,silent=TRUE)
 ob<-RTMB::MakeADFun(function(p)nll_fun(p,b),pb,silent=TRUE)
 expect_identical(pa,pb)
 expect_equal(oa$fn(oa$par),ob$fn(ob$par),tolerance=1e-10)
 expect_equal(oa$gr(oa$par),ob$gr(ob$par),tolerance=1e-10)
 ra<-oa$report();rb<-ob$report()
 expect_equal(ra$N,rb$N)
 expect_equal(ra$ssb,rb$ssb)
 expect_equal(rb$total_yield_pred,3*ra$total_yield_pred,ignore_attr=TRUE)
})

test_that('fleet forecasts, modes and information criteria preserve identities', {
 g<-expand.grid(year=2000:2039,age=1:7)
 obs<-multi_data(2000:2039,1:7)
 obs$weight$catch_weight<-1.25*obs$weight$obs
 args<-list(N_settings=list(process='off'),F_settings=list(process='off',mu_form=~fleet+logistic(age,by=fleet)),catch_settings=list(sd_form=~fleet),index_settings=list(q_form=~1))
 d<-do.call(prepare_tam,c(list(data=obs),args));p<-make_par(d)
 p$log_r0<-log(1000);p$log_r[]<-log(1000);p$log_sd_r<-log(.2)
 p$log_mu_f[]<-c(log(.25),log(1.6));p$log_sd_catch[]<-c(log(.15),0);p$log_sd_index[]<-log(.15)
 p[[grep('^F_a50_',names(p),value=TRUE)]][]<-c(3,4)
 p[[grep('^log_F_slope_',names(p),value=TRUE)]][]<-log(c(1,.8))
 set.seed(332);s<-nll_fun(p,d,simulate=TRUE)
 for(type in c('catch','index'))obs[[type]]$obs<-exp(s$log_obs[d$obs_map$type==type])
 f<-do.call(fit_tam,c(list(data=obs,silent=TRUE),args))
 before<-f$parameter_values
 expect_true(f$is_converged)
 expect_setequal(f$pop$F_fleet$fleet,c('A','B'))
 expect_true(all(is.finite(f$pop$F_fleet$se)))
 expect_true(all(is.finite(f$pop$mu_F_fleet$se)))
 expect_equal(attr(logLik(f),'df'),length(f$opt$par))
 changed<-f;changed$dat$obs_map$fleet[1]<-'Other'
 expect_error(AIC(f,changed),'identical|same')
 u<-update(f,start_par=f,proj_settings=list(n_proj=2,catch=c(0,100)),silent=TRUE)
 expect_true(u$is_converged)
 expect_equal(u$rep$total_yield_pred[41:42],c(0,100),tolerance=1e-7,ignore_attr=TRUE)
 expect_true(all(u$pop$F_fleet$upr[u$pop$F_fleet$year==2040]==0))
 for(fleet in u$dat$fleets)expect_equal(u$rep$F_fleet[42,,fleet]/u$rep$F_fleet[40,,fleet],rep(u$rep$F[42,1]/u$rep$F[40,1],7),ignore_attr=TRUE)
 controls<-.hindcast_controls(f,2039,.retro_projection(NULL))
 raw<-subset(obs$catch,year==2039);w<-subset(obs$weight,year==2039)
 expect_equal(controls$target,sum(raw$obs*w$catch_weight[match(raw$age,w$age)]))
 h<-fit_hindcast(f,folds=1,start_from_fit=TRUE,proj_settings=list(F_mult=1),progress=FALSE)
 expect_equal(sum(h$fold_diagnostics$is_converged),2)
 expect_setequal(na.omit(h$scores$fleet),c('A','B'))
 expect_setequal(na.omit(h$mohns_rho$fleet),c('A','B'))
 expect_identical(f$parameter_values,before)
 sim<-sim_tam(f,n=1,par_uncertainty='none',redraw_random=TRUE,progress=FALSE,seed=8)
 expect_setequal(sim$catch$fleet,c('A','B'))
 expect_true(nrow(sim$F_fleet)>0)
})

test_that('fleet warm starts align named states and labels retain leading zeros', {
 obs<-multi_data()
 d<-prepare_tam(obs,F_settings=list(process='iid',mu_form=~fleet+age))
 p<-make_par(d);p$log_f[]<-seq_along(p$log_f)/100
 o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
 f<-structure(list(dat=d,obj=o,parameter_values=o$par),class='tam_fit')
 target<-prepare_tam(obs,years=2000:2005,ages=2:5,F_settings=list(process='iid',mu_form=~fleet+age))
 merged<-.merge_start_par(make_par(target),.start_parameters(f,target))
 expect_equal(merged$log_f,p$log_f[rownames(merged$log_f),colnames(merged$log_f)])
 tab<-.tidy_process_deviations(structure(list(dat=d,rep=o$report(o$par),obj=o,parameter_values=o$par),class='tam_fit'))$F
 expect_equal(nrow(tab),length(p$log_f))
 expect_setequal(tab$fleet,c('A','B'))
 x<-array(1:4,c(2,1,2),dimnames=list(year=2000:2001,age=2,fleet=c('01','02')))
 expect_identical(unique(tidy_array(x)$fleet),c('01','02'))
})

test_that('catch-weight projections average by biological age', {
 obs<-multi_data()
 obs$weight$catch_weight<-obs$weight$obs*(1+(obs$weight$year-2000)/10)
 d<-prepare_tam(obs,ages=1:4,proj_settings=list(n_proj=2,n_mean=2))
 expected<-aggregate(catch_weight~age,subset(obs$weight,year %in% 2005:2006),mean)
 expect_equal(d$W_catch_plus_input[8,],expected$catch_weight[match(4:6,expected$age)],ignore_attr=TRUE)
 expect_equal(d$W_catch[8,1:3],expected$catch_weight[match(1:3,expected$age)],ignore_attr=TRUE)
 expect_equal(d$W_catch[9,],d$W_catch[8,],ignore_attr=TRUE)
})


test_that('warm starts omit coefficients when factor references change', {
 obs<-multi_data()
 d<-prepare_tam(obs,F_settings=list(process='iid',mu_form=~factor(age)))
 p<-make_par(d);p$log_q[]<-2;p$log_mu_f[]<-3
 o<-RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
 f<-structure(list(dat=d,obj=o,parameter_values=o$par),class='tam_fit')
 target<-prepare_tam(obs,ages=2:5,F_settings=list(process='iid',mu_form=~factor(age)))
 starts<-.start_parameters(f,target)
 expect_null(starts$log_q)
 expect_null(starts$log_mu_f)
 expect_true('log_f' %in% names(starts))
})
