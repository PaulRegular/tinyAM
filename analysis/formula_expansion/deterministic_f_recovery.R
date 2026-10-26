# Run from the repository root with the development package loaded.
deterministic_f_recovery <- function(n=20L, seed=1129L,
 output='analysis/formula_expansion/results/deterministic_f.rds') {
 grid<-expand.grid(year=2000:2039,age=1:7)
 data<-list(catch=transform(grid,obs=100,unit='thousand fish'),
  index=transform(grid,obs=100,survey='A',samp_time=.5,unit='thousand fish'),
  weight=transform(grid,obs=age^2/10,unit='kg'),maturity=transform(grid,obs=.5*plogis(age-3)))
 settings<-list(N_settings=list(process='off'),F_settings=list(process='off',mu_form=~logistic(age)),index_settings=list(q_form=~1))
 dat<-do.call(prepare_tam,c(list(data=data),settings));par<-make_par(dat)
 par$log_r0<-log(1000);par$log_r[]<-log(1000);par$log_sd_r<-log(.2)
 par$log_mu_f[]<-log(.4);par$log_sd_catch[]<-par$log_sd_index[]<-log(.15)
 par[grep('^F_a50_',names(par))]<-list(3)
 par[grep('^log_F_slope_',names(par))]<-list(log(.8))
 set.seed(seed);attempts<-list()
 for(i in seq_len(n)) {
  simulated<-nll_fun(par,dat,simulate=TRUE);obs<-dat$obs;warnings<-character()
  for(type in c('catch','index'))obs[[type]]$obs<-exp(simulated$log_obs[dat$obs_map$type==type])
  fit<-tryCatch(withCallingHandlers(do.call(fit_tam,c(list(data=obs,silent=TRUE),settings)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
  attempts[[i]]<-if(inherits(fit,'error'))list(error=conditionMessage(fit),warnings=warnings) else {
   p<-tinyAM:::.start_parameters(fit,fit$dat)
   list(checks=check_tam(fit),warnings=warnings,plateau=exp(p$log_mu_f),midpoint=unname(p[[grep('^F_a50_',names(p),value=TRUE)]]),slope=exp(p[[grep('^log_F_slope_',names(p),value=TRUE)]]))
  }
  dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(seed=seed,truth=par,attempts=attempts),output)
 }
 attr(attempts,"example_fit") <- if (inherits(fit,"tam_fit")) fit else NULL
 invisible(attempts)
}
