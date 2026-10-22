# Local feature trials only: no retained assessment recipe or aggregate export is changed.
fleet_recovery <- function(n=20L,seed=332L,
 output='analysis/formula_expansion/results/fleet_recovery.rds') {
 g<-expand.grid(year=2000:2039,age=1:7)
 data<-list(catch=rbind(transform(g,obs=100,fleet='A',unit='thousand fish'),transform(g,obs=100,fleet='B',unit='thousand fish')),
  index=transform(g,obs=100,survey='S',samp_time=.4,unit='thousand fish'),
  weight=transform(g,obs=age^2/10,catch_weight=age^2/8,unit='kg'),maturity=transform(g,obs=.5*plogis(age-3)))
 settings<-list(N_settings=list(process='off'),F_settings=list(process='off',mu_form=~fleet+logistic(age,by=fleet)),catch_settings=list(sd_form=~fleet),index_settings=list(q_form=~1))
 dat<-do.call(prepare_tam,c(list(data=data),settings));truth<-make_par(dat)
 truth$log_r0<-log(1000);truth$log_r[]<-log(1000);truth$log_sd_r<-log(.2)
 truth$log_mu_f[]<-c(log(.25),log(1.6));truth$log_sd_catch[]<-c(log(.15),0);truth$log_sd_index[]<-log(.15)
 for(name in grep('^F_a50_',names(truth),value=TRUE))truth[[name]][]<-c(3,4)
 for(name in grep('^log_F_slope_',names(truth),value=TRUE))truth[[name]][]<-log(c(1,.8))
 set.seed(seed);attempts<-list()
 for(i in seq_len(n)) {
  sim<-nll_fun(truth,dat,simulate=TRUE);obs<-dat$obs;warnings<-character()
  for(type in c('catch','index'))obs[[type]]$obs<-exp(sim$log_obs[dat$obs_map$type==type])
  fit<-tryCatch(withCallingHandlers(do.call(fit_tam,c(list(data=obs,silent=TRUE),settings)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
  attempts[[i]]<-if(inherits(fit,'error'))list(error=conditionMessage(fit),warnings=warnings) else {
   p<-tinyAM:::.start_parameters(fit,fit$dat)
   list(checks=check_tam(fit),warnings=warnings,
    plateau=exp(c(p$log_mu_f[1],sum(p$log_mu_f))),
    midpoint=unname(p[[grep('^F_a50_',names(p),value=TRUE)]]),
    slope=exp(p[[grep('^log_F_slope_',names(p),value=TRUE)]]))
  }
  dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(seed=seed,truth=truth,attempts=attempts),output)
 }
 attr(attempts,'example_fit')<-if(inherits(fit,'tam_fit'))fit else NULL
 invisible(attempts)
}

fleet_assessment_trial <- function(output='analysis/formula_expansion/results/fleet_assessment') {
 source('analysis/comp_assessments/R/run_assessment.R',local=environment())
 src<-read_assessment('nefsc_summer_flounder_2018',read_database())
 e<-new.env(parent=environment());e$source<-src;e$do_fit<-FALSE;e$silent<-TRUE
 sys.source('analysis/comp_assessments/scripts/translation/stocks/nefsc_summer_flounder_2018.R',e)
 catches<-subset(src$inputs,type=='catch' & measure=='numbers_at_age' & fleet!='Published total' & year %in% e$years)
 catches$age<-pmin(catches$age,max(e$ages))
 catches$obs<-catches$value*vapply(catches$unit,.translation_number_multiplier,numeric(1))/1000
 catch<-aggregate(obs~year+age+fleet,catches,sum,na.action=na.pass)
 catch$unit<-'thousand fish'
 stopifnot(nrow(catch)==length(e$years)*length(e$ages)*4L,!anyNA(catch$obs))
 obs<-e$obs;obs$catch<-catch
 check_tam_data(obs)
 warnings<-character()
 fit<-tryCatch(withCallingHandlers(fit_tam(obs,years=e$years,ages=e$ages,
  N_settings=list(process='iid',init='exp'),
  F_settings=list(process='ar1',mu_form=~0+fleet:factor(age),sd_form=~fleet,mean_ages=4),
  M_settings=list(process='off',mu_supplied=~M_assumption),
  catch_settings=list(sd_form=~1),index_settings=list(q_form=~0+q_key,sd_form=~1),silent=TRUE),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 if(inherits(fit,'error')) {
  result<-list(error=conditionMessage(fit),warnings=warnings)
 } else {
  background<-c('### Summer flounder: local four-fleet feature trial','',print_sources(src$assessment),'',
   'The 2018 ASAP run uses four fleets. This trial uses their reported catch numbers once; the separately published aggregate is excluded. Fleet and aggregate tables differ in some cells, so this is not an identical-observation comparison with the retained single-fishery recipe.',
   '', 'tinyAM uses four independent AR1 F surfaces with shared correlations, fleet-specific process SDs and fleet-by-age mean levels. N retains IID survival deviations and exponential initialization; M retains the published fixed age schedule. The same two Bigelow surveys, q structure and static SSB-weight approximation are retained. Catch weights fall back to those common stock weights; source fleet-specific biology and composition likelihoods are not reproduced.',
   '', 'No rescue fit or recipe change is made. A failed numerical fit remains exploratory.')
  metadata<-as.list(e);metadata$background<-background
  context<-.assessment_context(src,metadata)
  preview<-withCallingHandlers(view_assessment(fit,context,output_file=file.path(output,'dashboard.html'),open_file=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  result<-list(fit=fit,checks=check_tam(fit),warnings=warnings,summary=preview$summary,source=src$assessment)
 }
 saveRDS(result,file.path(output,'attempt.rds'))
 invisible(result)
}
