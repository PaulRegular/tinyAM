# Run from the repository root. Every configured attempt is cached, including
# warnings/errors. No retries, optimizer changes, or selection of successful seeds.
pkgload::load_all(".", quiet = TRUE)

fleet_recovery_configs <- function() {
  full <- data.frame(id="scaa_double_logistic",family="scaa",shape="double_logistic",
    coverage="full",process="off",correlation="shared",seed_group=1L)
  full <- rbind(full,data.frame(id=c("sam_shared","sam_separate"),family="sam",shape="states",
    coverage="full",process="cor_rw",correlation=c("shared","separate"),seed_group=5:6))
  full <- rbind(full,data.frame(id=c("smooth_off","smooth_iid","smooth_ar1","smooth_ar1_weak"),
    family="smooth",shape="smooth",coverage=c("full","full","full","weak_separation"),
    process=c("off","iid","ar1","ar1"),correlation="shared",seed_group=7:10))
  isolated <- data.frame(id=c("logistic","double_logistic","double_logistic_reversed"),
    family="curve",shape=c("logistic","double_logistic","reversed"),coverage="full",seed_group=1:3)
  isolated <- rbind(isolated,data.frame(id=c("smooth_tp","smooth_cr","smooth_by","smooth_tensor"),
    family="smooth",shape=c("tp","cr","by","tensor"),coverage="full",seed_group=7:10))
  list(isolated=isolated,full=full)
}

.fleet_scenario <- function(config, seed) {
  years <- 2000:2011; ages <- 1:10; horizon <- 2L
  grid <- expand.grid(year=2000:2013,age=ages)
  stock <- transform(grid,obs=.2*age^1.5,type="stock",fleet=NA_character_,unit="kg")
  data <- list(catch=rbind(transform(grid,obs=1,fleet="Trawl",unit="thousand fish"),
      transform(grid,obs=1,fleet="Gillnet",unit="thousand fish")),
    index=transform(grid,obs=1,survey="RV",samp_time=.4,unit="fish per tow"),
    weight=rbind(stock,transform(stock,type="catch",fleet="Trawl",obs=.9*obs),
      transform(stock,type="catch",fleet="Gillnet",obs=1.15*obs)),
    maturity=transform(grid,obs=plogis(age-4)))
  form <- ~0+fleet:factor(year)+for_fleet("Trawl",~0+logistic(age))+
    for_fleet("Gillnet",~0+double_logistic(age))
  settings <- list(years=years,ages=ages,N_settings=list(process="off"),
    rec_settings=list(process="iid",mu_form=~1),
    F_settings=list(process=config$process,mu_form=form,mean_ages=3:8),
    M_settings=list(process="off",mu_supplied=~I(.2),mean_ages=3:8),
    catch_settings=list(sd_form=~fleet),index_settings=list(q_form=~logistic(age)))
  if (config$family == "sam") {
    settings$N_settings <- list(process="iid",init="latent")
    settings$rec_settings <- list(process="rw",init="latent")
    settings$F_settings <- list(process="cor_rw",mu_form=NULL,sd_form=~fleet,
      cor_form=if(config$correlation == "shared")~1 else ~fleet,
      age_breaks=list(Trawl=c(1:9,11),Gillnet=c(2:8,11)),mean_ages=3:8)
    data$catch$obs[data$catch$fleet == "Gillnet" & data$catch$age == 1] <- 0
  }
  if (config$family == "smooth") {
    if (config$process != "off") settings$F_settings$sd_form <- ~fleet
    settings$F_settings$mu_form <- if (config$coverage == "weak_separation")
      ~fleet+s(age,k=5,by=fleet)+s(year,k=8,by=fleet) else
      ~fleet+s(age,k=5,by=fleet)+s(year,k=5,by=fleet)
  }
  # Literal study formulas need no caller locals. Keep cached fits/worker payloads
  # independent of this generator's RTMB objects and temporary data.
  for (role in names(settings)) if (is.list(settings[[role]])) {
    for (name in names(settings[[role]])) if (inherits(settings[[role]][[name]],"formula")) {
      environment(settings[[role]][[name]]) <- asNamespace("tinyAM")
    }
  }
  dat <- do.call(tinyAM::prepare_tam,c(list(data=data,proj_settings=list(n_proj=horizon,n_mean=3,F_mult=1)),settings))
  p <- tinyAM::make_par(dat)
  p$log_r0[] <- log(1000)
  if (!is.null(p$log_r)) p$log_r[] <- log(1000)
  if (!is.null(p$rec_beta)) p$rec_beta[] <- log(1000)
  p$log_sd_r[] <- log(.2)
  if (!is.null(p$log_n0)) p$log_n0[] <- log(1000)-.35*(ages[-1]-1)
  if (!is.null(p$log_sd_n)) p$log_sd_n[] <- log(.08)
  p$log_sd_catch[] <- qr.solve(dat$sd_catch_modmat,
    log(ifelse(dat$obs$catch$fleet == "Trawl",.12,.15)))
  p$log_sd_index[] <- log(.12); p$log_q[] <- log(.03)
  for (term in dat$q_terms) {
    p[[term$a50_parameter]][] <- 3
    p[[term$slope_parameter]][] <- log(1.1)
  }
  if (config$family == "sam") {
    states <- dat$F_states
    rates <- ifelse(states$fleet == "Trawl",.25*plogis((states$age-3)/1.5),
      .15*exp(-.5*((states$age-6)/2)^2))
    p$log_f[] <- rep(log(rates),each=nrow(p$log_f))
    p$atanh_rho_f[] <- atanh(if (config$correlation == "shared") .5 else
      ifelse(grepl("Trawl",names(p$atanh_rho_f)),.7,.2))
  } else {
    X <- dat$F_modmat
    mean <- log(ifelse(dat$F_data$fleet == "Trawl",.25,.12))
    if (config$family != "smooth") mean <- mean+.2*sin((dat$F_data$year-2000)/2.5)
    historical <- !dat$is_proj[match(dat$F_data$year,dat$years)]
    p$log_mu_f[] <- qr.solve(X[historical,,drop=FALSE],mean[historical])
    for (term in dat$F_terms) {
      if (term$type == "logistic") {
        p[[term$a50_parameter]][] <- 3; p[[term$slope_parameter]][] <- log(1.2)
      } else if (term$type == "double_logistic") {
        p[[term$a50_up_parameter]][] <- 3; p[[term$scale_up_parameter]][] <- log(1)
        p[[term$a50_down_parameter]][] <- 8; p[[term$scale_down_parameter]][] <- log(1.3)
      } else if (term$type == "smooth") {
        p[[term$sd_parameter]][] <- log(if (term$variable == "age") .2 else .15)
        p$log_mu_f[colnames(term$smooth_fixed)] <- if (term$variable == "age") .2 else .05
      }
    }
  }
  if (config$process != "off") {
    design <- dat$process_sd$F
    desired <- log(ifelse(dat$F_states$fleet == "Trawl",.1,.14))
    p[[design$parameter]][] <- if (design$default)log(.12) else qr.solve(design$matrix,desired)
    if (config$process == "ar1") p$logit_phi_f[] <- stats::qlogis(
      ifelse(grepl("age",names(p$logit_phi_f)),
        if(config$coverage == "weak_separation").7 else .2,
        if(config$coverage == "weak_separation").9 else .4))
  }
  set.seed(seed)
  draw <- tinyAM::nll_fun(p,dat,simulate=TRUE)
  for (name in intersect(names(draw),names(p))) if (name != "log_obs") p[[name]] <- draw[[name]]
  object <- RTMB::MakeADFun(function(p)tinyAM::nll_fun(p,dat),p,silent=TRUE)
  truth <- object$report()
  log_obs <- draw$log_obs
  future <- dat$is_proj[match(dat$obs_map$year,dat$years)]
  finite <- future & is.finite(truth$log_pred)
  set.seed(seed+9000000L)
  log_obs[finite] <- stats::rnorm(sum(finite),truth$log_pred[finite],truth$sd_obs[finite])
  for (type in c("catch","index")) {
    data[[type]] <- dat$obs[[type]]
    data[[type]]$obs <- exp(log_obs[dat$obs_map$type == type])
    data[[type]]$is_proj <- FALSE
  }
  if (config$family == "sam") data$catch$obs[data$catch$fleet == "Gillnet" & data$catch$age == 1] <- 0
  if (config$coverage == "weak") for (type in c("catch","index")) {
    remove <- data[[type]]$year %in% years & data[[type]]$age >= 8 & (data[[type]]$year-2000) %% 4 != 0
    data[[type]]$obs[remove] <- NA_real_
  }
  list(data=data,settings=settings,parameters=p,truth=truth,truth_data=dat,
    seed=seed,forecast_definition="Two years conditional on status-quo terminal fleet F; no new F innovations.")
}

.fleet_isolated <- function(config, seed) {
  set.seed(seed)
  d <- expand.grid(age=1:12,year=2000:2007)
  d$group <- factor(rep(c("A","B"),length.out=nrow(d)))
  d$obs <- 1; d$is_proj <- FALSE
  if (config$family == "curve") {
    type <- if(config$shape == "reversed")"double_logistic" else config$shape
    term <- tinyAM:::.make_q_term(list(type=type,variable="age",by=NULL),d)
    parameters <- tinyAM:::.curve_parameters(term)
    if (type == "logistic") {
      parameters[[term$a50_parameter]][] <- 3; parameters[[term$slope_parameter]][] <- log(1.2)
    } else if (type == "double_logistic") {
      parameters[[term$a50_up_parameter]][] <- if(config$shape == "reversed")8 else 3
      parameters[[term$a50_down_parameter]][] <- if(config$shape == "reversed")3 else 8
      parameters[[term$scale_up_parameter]][] <- log(1)
      parameters[[term$scale_down_parameter]][] <- log(1.3)
    }
    truth <- log(.4)+tinyAM:::.curve_log_selectivity(d$age,parameters,term,1)
    y <- stats::rnorm(nrow(d),truth,.15)
    keep <- rep(TRUE,nrow(d))
    if (config$coverage == "weak") keep <- d$age < 8 | (d$year-2000) %% 4 == 0
    start <- c(list(level=log(.3)),tinyAM:::.curve_parameters(term))
    o <- RTMB::MakeADFun(function(p)-sum(RTMB::dnorm(y[keep],p$level+
      tinyAM:::.curve_log_selectivity(d$age[keep],p,term,1),.15,log=TRUE)),start,silent=TRUE)
    lower <- rep(-Inf,length(o$par)); upper <- rep(Inf,length(o$par))
    opt <- stats::nlminb(o$par,o$fn,o$gr,lower=lower,upper=upper,control=list(eval.max=1000,iter.max=1000))
    o$fn(opt$par)
    gradient <- o$gr(opt$par)
    projected <- tinyAM:::.project_gradient(gradient,opt$par,list(lower=lower,upper=upper))$gradient
    V <- tryCatch(solve(o$he(opt$par)),error=function(e)NULL)
    pd <- !is.null(V) && all(is.finite(V)) && all(eigen(V,symmetric=TRUE,only.values=TRUE)$values > 0)
    estimates <- o$env$parList(x=opt$par)
    prediction <- estimates$level+tinyAM:::.curve_log_selectivity(d$age,estimates,term,1)
    true <- c(list(level=log(.4)),parameters)
    true <- unlist(true,use.names=FALSE)
    se <- if(pd)sqrt(diag(V)) else rep(NA_real_,length(opt$par))
    function_se <- rep(NA_real_,nrow(d))
    if (pd) {
      lengths <- lengths(start); offsets <- cumsum(c(0L,lengths))
      tape <- RTMB::MakeTape(function(v) {
        parameters <- setNames(lapply(seq_along(lengths),function(i)
          v[seq.int(offsets[i]+1L,offsets[i+1L])]),names(start))
        parameters$level+tinyAM:::.curve_log_selectivity(d$age,parameters,term,1)
      },opt$par)
      J <- tape$jacobian(opt$par)
      function_se <- sqrt(pmax(0,rowSums((J %*% V)*J)))
    }
    return(list(opt=opt,max_gradient=max(abs(projected)),pdHess=pd,
      is_converged=opt$convergence == 0 && max(abs(projected)) <= .01 && pd,
      parameters=data.frame(parameter=names(opt$par),truth=true,estimate=opt$par,se=se,
        covered=abs(opt$par-true) <= stats::qnorm(.975)*se),
      function_log_RMSE=sqrt(mean((prediction-truth)^2)),
      function_coverage=mean(abs(prediction-truth) <= stats::qnorm(.975)*function_se),seed=seed))
  }
  formula <- switch(config$shape,tp=~s(age,k=5),cr=~s(age,k=5,bs="cr"),
    by=~group+s(age,k=5,by=group),tensor=~t2(age,year,k=c(4,4)))
  design <- tinyAM:::.parse_q_formula(formula,d)
  p <- c(list(beta=seq_len(ncol(design$q_modmat))/10),tinyAM:::.q_term_parameters(design$q_terms))
  for (term in design$q_terms) {
    p[[term$sd_parameter]][] <- log(.3)
    p[[term$parameter]][] <- stats::rnorm(length(p[[term$parameter]]),0,.3)
  }
  truth <- drop(design$q_modmat %*% p$beta)+tinyAM:::.q_effects(p,list(obs=list(index=d),q_terms=design$q_terms))$contribution
  y <- stats::rnorm(nrow(d),truth,.15)
  start <- c(list(beta=rep(0,ncol(design$q_modmat))),tinyAM:::.q_term_parameters(design$q_terms))
  o <- RTMB::MakeADFun(function(p) {
    effect <- tinyAM:::.q_effects(p,list(obs=list(index=d),q_terms=design$q_terms))
    value <- drop(design$q_modmat %*% p$beta)+effect$contribution
    RTMB::ADREPORT(value)
    effect$nll-sum(RTMB::dnorm(y,value,.15,log=TRUE))
  },start,random=vapply(design$q_terms,`[[`,character(1),"parameter"),silent=TRUE)
  opt <- stats::nlminb(o$par,o$fn,o$gr,control=list(eval.max=1000,iter.max=1000))
  o$fn(opt$par); gradient <- o$gr(opt$par)
  sd <- tryCatch(RTMB::sdreport(o,par.fixed=opt$par,getReportCovariance=FALSE),error=function(e)NULL)
  pd <- !is.null(sd) && isTRUE(sd$pdHess) && all(is.finite(sd$cov.fixed))
  estimate <- o$env$parList(par=o$env$last.par)
  prediction <- drop(design$q_modmat %*% estimate$beta)+
    tinyAM:::.q_effects(estimate,list(obs=list(index=d),q_terms=design$q_terms))$contribution
  truth_par <- c(p$beta,rep(log(.3),length(design$q_terms)))
  se <- if(pd)sqrt(diag(sd$cov.fixed)) else rep(NA_real_,length(opt$par))
  function_se <- if(pd)as.list(sd,"Std. Error",report=TRUE)$value else rep(NA_real_,nrow(d))
  list(opt=opt,max_gradient=max(abs(gradient)),pdHess=pd,
    is_converged=opt$convergence == 0 && max(abs(gradient)) <= .01 && pd,
    parameters=data.frame(parameter=names(opt$par),truth=truth_par,estimate=opt$par,se=se,
      covered=abs(opt$par-truth_par) <= stats::qnorm(.975)*se),
    function_log_RMSE=sqrt(mean((prediction-truth)^2)),
    function_coverage=mean(abs(prediction-truth) <= stats::qnorm(.975)*function_se),seed=seed)
}

.fleet_full <- function(config, seed) {
  example <- .fleet_scenario(config,seed)
  fit <- do.call(tinyAM::fit_tam,c(list(data=example$data,silent=TRUE),example$settings))
  checks <- tinyAM::check_tam(fit)
  out <- list(fit=fit,example=example,diagnostics=checks,
    is_converged=checks$is_converged,max_gradient=checks$max_gradient,
    pdHess=checks$pd_hessian,seed=seed)
  # Keep unsuccessful fits; their numerical estimates do not enter recovery scores.
  if (!checks$is_converged) return(out)
  tryCatch({
  tables <- tinyAM::tidy_pop(fit)
  truth <- tinyAM::tidy_rep(list(dat=example$truth_data,rep=example$truth))
  metrics <- lapply(intersect(c("F","F_fleet","N","ssb","recruitment"),names(tables)),function(name) {
    keys <- intersect(c("year","age","fleet"),names(tables[[name]]))
    d <- merge(tables[[name]],truth[[name]][c(keys,"est")],by=keys,suffixes=c("","_truth"))
    d <- d[!d$is_proj & is.finite(d$est) & d$est_truth > 0,,drop=FALSE]
    interval <- all(c("lwr","upr") %in% names(d))
    data.frame(quantity=name,n=nrow(d),relative_bias=mean((d$est-d$est_truth)/d$est_truth),
      log_RMSE=sqrt(mean((log(d$est)-log(d$est_truth))^2)),
      coverage=if(interval)mean(d$lwr <= d$est_truth & d$upr >= d$est_truth) else NA_real_)
  })
  out$surfaces <- do.call(rbind,metrics)
  est <- tinyAM:::.tam_parameter_summary(fit,"Estimate")
  se <- tinyAM:::.tam_parameter_summary(fit,"Std. Error")
  # Raw fitted scales retain SD/correlation transformations and signed coefficients.
  keep <- names(est)[grepl("sd|phi|rho|a50|scale|slope",names(est))]
  keep <- intersect(keep,names(example$parameters))
  out$parameters <- do.call(rbind,lapply(keep,function(name) {
    target <- example$parameters[[name]]; estimate <- est[[name]]
    if (length(target) != length(estimate)) return(NULL)
    if (!is.null(names(target)) && !identical(names(target),names(estimate))) return(NULL)
    data.frame(parameter=paste(name,if(is.null(names(estimate)))seq_along(estimate) else names(estimate),sep=":"),
      truth=as.vector(target),estimate=as.vector(estimate),se=as.vector(se[[name]]),
      covered=abs(estimate-target) <= stats::qnorm(.975)*se[[name]])
  }))
  projected <- update(fit,start_par=fit,proj_settings=list(n_proj=2,n_mean=3,F_mult=1),silent=TRUE)
  out$forecast_fit <- projected
  out$forecast_diagnostics <- tinyAM::check_tam(projected)
  if (!out$forecast_diagnostics$is_converged) return(out)
  prediction <- tinyAM::tidy_obs_pred(projected)
  validation <- merge(prediction$index[prediction$index$is_proj,c("year","age","survey","pred")],
    example$data$index[c("year","age","survey","obs")],by=c("year","age","survey"))
  history <- example$data$index[example$data$index$year <= max(example$settings$years) &
    is.finite(example$data$index$obs) & example$data$index$obs > 0,,drop=FALSE]
  history <- history[order(history$year,decreasing=TRUE),,drop=FALSE]
  history <- history[!duplicated(history[c("survey","age")]),,drop=FALSE]
  names(history)[names(history) == "obs"] <- "baseline"
  validation <- merge(validation,history[c("survey","age","baseline")],by=c("survey","age"),all.x=TRUE)
  valid <- with(validation,is.finite(obs) & obs > 0 & is.finite(pred) & pred > 0)
  validation$error <- with(validation,log(obs)-log(pred))
  validation$persistence_error <- with(validation,log(obs)-log(baseline))
  validation$horizon <- validation$year-max(example$settings$years)
  out$forecast_cases <- validation
  scores <- lapply(sort(unique(validation$horizon)),function(h) {
    d <- validation[valid & validation$horizon == h,,drop=FALSE]
    # MASE averages age-specific ratios, not a ratio pooled across ages.
    ratios <- vapply(split(d,d$age),function(x) {
      pair <- is.finite(x$error) & is.finite(x$persistence_error)
      den <- sum(abs(x$persistence_error[pair]))
      if (sum(pair) && den > 0)sum(abs(x$error[pair]))/den else NA_real_
    },numeric(1))
    data.frame(horizon=h,n=nrow(d),log_RMSE=sqrt(mean(d$error^2)),
      MASE=if(any(is.finite(ratios)))mean(ratios[is.finite(ratios)]) else NA_real_)
  })
  out$forecast <- do.call(rbind,scores)
  catch <- merge(prediction$catch[prediction$catch$is_proj,c("year","age","fleet","pred")],
    example$data$catch[c("year","age","fleet","obs")],by=c("year","age","fleet"))
  out$composition <- do.call(rbind,lapply(split(catch,interaction(catch$year,catch$fleet)),function(d) {
    if (any(!is.finite(d$obs) | !is.finite(d$pred)) || sum(d$obs) <= 0 || sum(d$pred) <= 0)return(NULL)
    p <- d$pred/sum(d$pred); observed <- d$obs/sum(d$obs)
    data.frame(year=d$year[1],fleet=d$fleet[1],RMSE=sqrt(mean((p-observed)^2)),
      MAE=mean(abs(p-observed)),mean_age_error=sum(d$age*(p-observed)))
  }))
  out
  },error=function(e) {out$evaluation_error <- conditionMessage(e); out})
}

.fleet_attempt <- function(job, directory) {
  folder <- file.path(directory,job$kind,job$config$id)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  path <- file.path(folder,sprintf("%03d.rds",job$replicate))
  if (file.exists(path)) return(path)
  warnings <- messages <- character(); error <- NULL
  elapsed <- system.time(value <- withCallingHandlers(
    tryCatch(if(job$kind == "isolated").fleet_isolated(job$config,job$seed) else
      .fleet_full(job$config,job$seed),error=function(e) {error <<- conditionMessage(e); NULL}),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")},
    message=function(m){messages <<- c(messages,conditionMessage(m)); invokeRestart("muffleMessage")}))[["elapsed"]]
  record <- list(job=job,value=value,error=error,warnings=warnings,messages=messages,
    elapsed=elapsed,session=utils::sessionInfo())
  saveRDS(record,path)
  cat(job$kind,job$config$id,job$replicate,if(!is.null(error))"ERROR" else
    if(isTRUE(value$is_converged))"PASS" else "NUMERICAL FAILURE",round(elapsed,1),"s\n")
  path
}

run_fleet_recovery <- function(workers=3L,directory="analysis/formula_expansion/results/multi_fleet/recovery_supported") {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  script <- normalizePath("analysis/formula_expansion/multi_fleet_recovery.R",winslash="/")
  sources <- c(script,list.files("R",pattern="[.]R$",full.names=TRUE),"NAMESPACE")
  manifest <- list(hashes=tools::md5sum(sources),configurations=fleet_recovery_configs(),
    isolated_replicates=100L,full_replicates=30L,runtime=c(R=as.character(getRversion()),
      RTMB=as.character(packageVersion("RTMB")),mgcv=as.character(packageVersion("mgcv"))))
  file <- file.path(directory,"manifest.rds")
  if (file.exists(file) && !identical(readRDS(file),manifest)) stop("Study sources changed; archive the previous run before starting another.")
  saveRDS(manifest,file)
  jobs <- list()
  for (kind in names(manifest$configurations)) {
    configs <- manifest$configurations[[kind]]
    for (i in seq_len(nrow(configs))) for (r in seq_len(if(kind == "isolated")100L else 30L)) {
      jobs[[length(jobs)+1L]] <- list(kind=kind,config=configs[i,,drop=FALSE],replicate=r,
        seed=(if(kind == "isolated")500000L else 800000L)+1000L*configs$seed_group[i]+r)
    }
  }
  # Interleave configurations so numerical failures cannot select a scenario out.
  jobs <- jobs[order(vapply(jobs,`[[`,integer(1),"replicate"))]
  cluster <- parallel::makePSOCKcluster(workers,outfile=file.path(directory,"workers.log"))
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster,function(root,script,libs) {
    .libPaths(libs); setwd(root); source(script,local=.GlobalEnv); NULL
  },getwd(),script,.libPaths())
  paths <- parallel::parLapplyLB(cluster,jobs,function(job,directory).fleet_attempt(job,directory),directory=directory)
  saveRDS(unlist(paths),file.path(directory,"attempts.rds"))
  summarize_fleet_recovery(directory)
}

summarize_fleet_recovery <- function(directory="analysis/formula_expansion/results/multi_fleet/recovery_supported") {
  paths <- list.files(directory,pattern="^[0-9]{3}[.]rds$",recursive=TRUE,full.names=TRUE)
  summary <- parameters <- surfaces <- forecasts <- compositions <- list()
  for (file in paths) {
    r <- readRDS(file); v <- r$value
    identity <- data.frame(kind=r$job$kind,scenario=r$job$config$id,replicate=r$job$replicate)
    summary[[length(summary)+1L]] <- cbind(identity,success=isTRUE(v$is_converged),
      error=if(is.null(r$error))NA_character_ else r$error,warnings=length(r$warnings),
      elapsed=r$elapsed,max_gradient=if(is.null(v$max_gradient))NA_real_ else v$max_gradient,
      function_log_RMSE=if(is.null(v$function_log_RMSE))NA_real_ else v$function_log_RMSE,
      function_coverage=if(is.null(v$function_coverage))NA_real_ else v$function_coverage)
    if (!isTRUE(v$is_converged)) next
    for (name in c("parameters","surfaces","forecast","composition")) {
      table <- v[[name]]
      if (is.null(table) || !nrow(table)) next
      table <- cbind(identity,table)
      if (name == "parameters") parameters[[length(parameters)+1L]] <- table
      if (name == "surfaces") surfaces[[length(surfaces)+1L]] <- table
      if (name == "forecast") forecasts[[length(forecasts)+1L]] <- table
      if (name == "composition") compositions[[length(compositions)+1L]] <- table
    }
  }
  tables <- list(attempts=do.call(rbind,summary),parameters=do.call(rbind,parameters),
    surfaces=do.call(rbind,surfaces),forecasts=do.call(rbind,forecasts),compositions=do.call(rbind,compositions))
  for (name in names(tables)) if(!is.null(tables[[name]]))utils::write.csv(tables[[name]],file.path(directory,paste0(name,".csv")),row.names=FALSE)
  by <- split(tables$attempts,interaction(tables$attempts$kind,tables$attempts$scenario,drop=TRUE))
  result <- do.call(rbind,lapply(by,function(d) {
    p <- mean(d$success)
    interval <- stats::binom.test(sum(d$success),nrow(d))$conf.int
    data.frame(kind=d$kind[1],scenario=d$scenario[1],attempts=nrow(d),successes=sum(d$success),
      success_rate=p,success_MCSE=sqrt(p*(1-p)/nrow(d)),success_lower=interval[1],
      success_upper=interval[2],errors=sum(!is.na(d$error)),
      warned=sum(d$warnings > 0),function_log_RMSE=mean(d$function_log_RMSE[d$success],na.rm=TRUE),
      function_coverage=mean(d$function_coverage[d$success],na.rm=TRUE))
  }))
  utils::write.csv(result,file.path(directory,"summary.csv"),row.names=FALSE)
  saveRDS(tables,file.path(directory,"tables.rds")); summarize_fleet_metrics(directory); print(result)
  invisible(tables)
}

summarize_fleet_metrics <- function(directory="analysis/formula_expansion/results/multi_fleet/recovery_supported") {
  tables <- readRDS(file.path(directory,"tables.rds"))
  # Each entry being averaged below is one replicate, not an independent
  # age/year cell. MCSE therefore describes between-replicate uncertainty.
  average <- function(x) {
    x <- x[is.finite(x)]
    c(mean=if(length(x))mean(x) else NA_real_,
      MCSE=if(length(x)>1L)stats::sd(x)/sqrt(length(x)) else NA_real_,n=length(x))
  }
  summarize <- function(data,groups,metrics,file) {
    if (is.null(data) || !nrow(data)) return(NULL)
    key <- do.call(interaction,c(data[groups],drop=TRUE,lex.order=TRUE))
    out <- do.call(rbind,lapply(split(data,key),function(d) {
      values <- unlist(lapply(metrics,function(name)stats::setNames(average(d[[name]]),
        paste(name,c("mean","MCSE","n"),sep="_"))))
      cbind(d[1,groups,drop=FALSE],as.data.frame(as.list(values)))
    }))
    utils::write.csv(out,file.path(directory,file),row.names=FALSE); out
  }
  summarize(tables$surfaces,c("scenario","quantity"),c("relative_bias","log_RMSE","coverage"),"surface_summary.csv")
  parameters <- tables$parameters
  parameters$bias <- parameters$estimate-parameters$truth
  summarize(parameters,c("kind","scenario","parameter"),c("truth","estimate","bias","covered"),"parameter_summary.csv")
  summarize(tables$forecasts,c("scenario","horizon"),c("log_RMSE","MASE"),"forecast_summary.csv")
  if (!is.null(tables$compositions)) {
    # Average fleet/horizon cases first, preserving one unit per replicate.
    composition <- stats::aggregate(tables$compositions[c("RMSE","MAE","mean_age_error")],
      tables$compositions[c("scenario","replicate")],mean)
    summarize(composition,"scenario",c("RMSE","MAE","mean_age_error"),"composition_summary.csv")
  }
  invisible(tables)
}

.fleet_N_coverage <- function(fit,example) {
  # Add a report to a private copy of the fitted objective; do not optimize or
  # modify package mathematics, parameters, or the saved fit's objective state.
  fn <- tinyAM::nll_fun
  code <- as.list(body(fn))
  code[[length(code)]] <- quote({study_log_N <- log_N; RTMB::ADREPORT(study_log_N); jnll})
  body(fn) <- as.call(code)
  dat <- fit$dat
  parameters <- fit$obj$env$parList(par=fit$parameter_values)
  o <- RTMB::MakeADFun(function(p)fn(p,dat),parameters,random=fit$obj$env$.random,
    map=fit$parameter_map,silent=TRUE)
  value <- o$fn(fit$opt$par)
  difference <- value-fit$opt$objective
  if (abs(difference) > 1e-7)stop("Added N report changed the fitted objective.")
  sdr <- RTMB::sdreport(o,par.fixed=fit$opt$par,getReportCovariance=FALSE)
  estimate <- as.list(sdr,"Estimate",report=TRUE)$study_log_N
  se <- as.list(sdr,"Std. Error",report=TRUE)$study_log_N
  target <- log(example$truth$N[as.character(dat$years),as.character(dat$ages)])
  list(coverage=mean(abs(estimate-target) <= stats::qnorm(.975)*se),objective_difference=difference)
}

postprocess_fleet_recovery <- function(directory="analysis/formula_expansion/results/multi_fleet/recovery_supported") {
  paths <- list.files(file.path(directory,"full"),pattern="^[0-9]{3}[.]rds$",recursive=TRUE,full.names=TRUE)
  for (file in paths) {
    r <- readRDS(file)
    if (!isTRUE(r$value$is_converged) || !is.null(r$value$N_uncertainty)) next
    warnings <- character()
    result <- withCallingHandlers(tryCatch(.fleet_N_coverage(r$value$fit,r$value$example),
      error=function(e)list(error=conditionMessage(e))),warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
    result$warnings <- warnings; r$value$N_uncertainty <- result
    if (!is.null(result$coverage)) r$value$surfaces$coverage[r$value$surfaces$quantity == "N"] <- result$coverage
    saveRDS(r,file)
  }
  summarize_fleet_recovery(directory)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) && args[1L] == "metrics") summarize_fleet_metrics() else
    if (length(args) && args[1L] == "postprocess") postprocess_fleet_recovery() else
    if (length(args) && args[1L] == "summary") summarize_fleet_recovery() else
    run_fleet_recovery(if(length(args))as.integer(args[1L]) else 3L)
}
