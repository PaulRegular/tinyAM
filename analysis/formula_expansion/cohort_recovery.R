# Run from the repository root after loading tinyAM. Generated attempts stay local.
cohort_recovery <- function(n = 10L, seed = 2718L,
                            output = 'analysis/formula_expansion/results/cohort_recovery.rds') {
  grid <- expand.grid(year = 2000:2029, age = 1:6)
  data <- list(
    catch = transform(grid, obs = 100, unit = 'thousand fish'),
    index = transform(grid, obs = 100, survey = 'A', samp_time = .5, unit = 'thousand fish', temperature = sin((year - 2000) / 3)),
    weight = transform(grid, obs = age^2 / 10, unit = 'kg'),
    maturity = transform(grid, obs = .5 * plogis(age - 3)),
    pop_covariates = transform(grid, temperature = sin((year - 2000) / 3)))
  settings <- list(N_settings = list(process = 'iid', mu_form = ~ 0 + temperature),
    F_settings = list(process = 'iid', mu_form = ~factor(age)),
    M_settings = list(mu_supplied = ~I(.2)),
    index_settings = list(q_form = ~factor(age)))
  dat <- do.call(prepare_tam, c(list(data = data), settings))
  truth <- make_par(dat)
  truth$log_r0 <- log(1000)
  truth$log_r[] <- log(1000)
  truth$log_sd_r <- log(.15)
  truth$mu_n[] <- .12
  truth$log_mu_f[] <- c(log(.25), rep(0, 5))
  truth$log_sd_f <- log(.08)
  truth$log_sd_n <- log(.06)
  truth$log_sd_catch[] <- log(.12)
  truth$log_sd_index[] <- log(.12)
  set.seed(seed)
  attempts <- list()
  for (scenario in c('supported', 'mortality_covariate', 'survey_covariate')) {
    model <- settings
    if (scenario == 'mortality_covariate') model$M_settings$mu_form <- ~0+temperature
    if (scenario == 'survey_covariate') model$index_settings$q_form <- ~factor(age)+temperature
    for (i in seq_len(n)) {
      warnings <- character()
      sim <- nll_fun(truth, dat, simulate = TRUE)
      simulated_par <- modifyList(truth, sim[intersect(names(sim),names(truth))])
      true_log_N <- tinyAM:::.population_states(simulated_par,dat,exp(sim$log_f)+.2)$log_N
      obs <- dat$obs
      for (type in c('catch','index')) obs[[type]]$obs <- exp(sim$log_obs[dat$obs_map$type == type])
      fit <- tryCatch(withCallingHandlers(
        do.call(fit_tam, c(list(data=obs, silent=TRUE), model)),
        warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}), error=identity)
      key <- paste(scenario,i,sep=':')
      if (inherits(fit,'error')) {
        attempts[[key]] <- list(scenario=scenario,replicate=i,error=conditionMessage(fit),warnings=warnings)
      } else {
        p <- tinyAM:::.start_parameters(fit,fit$dat)
        se <- if (is.null(fit$sdrep)) NA_real_ else as.list(fit$sdrep,'Std. Error')$mu_n
        attempts[[key]] <- list(scenario=scenario,replicate=i,warnings=warnings,
          checks=check_tam(fit), estimate=p$mu_n, se=se,
          covered=abs(p$mu_n-.12)<=1.96*se,
          abundance_rmse=sqrt(mean((log(fit$rep$N)-true_log_N)^2)))
      }
      dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
      saveRDS(list(seed=seed,truth=truth,attempts=attempts),output)
    }
  }
  invisible(attempts)
}

# Known-survivor Gaussian check: 100 independent datasets with an exact coefficient.
cohort_gaussian_recovery <- function(n = 100L, seed = 3141L) {
  set.seed(seed)
  x <- sin(seq_len(120) / 7)
  do.call(rbind, lapply(seq_len(n), function(i) {
    response <- .12*x + rnorm(length(x),sd=.06)
    fit <- lm(response ~ 0 + x)
    estimate <- coef(fit)[1]
    se <- sqrt(vcov(fit)[1,1])
    data.frame(replicate=i,estimate=unname(estimate),se=se,covered=abs(estimate-.12)<=1.96*se)
  }))
}
