source("analysis/formula_expansion/recruitment_validation.R")

rec_cross_cases <- c("rw", rec_cases[1:4])
rec_cross_formula <- function(case) if (case == "rw") ~ rw(year) else rec_formula(case)
rec_cross_settings <- function(form, obs) {
  settings <- rec_settings(form)
  if ("f_pressure" %in% names(obs$catch)) settings$F_settings$mu_form <- ~ factor(age) + f_pressure
  settings
}
rec_cross_parameters <- function(generated, fitted) {
  curve <- function(case) strsplit(case, "_")[[1L]][1L]
  parameters <- if (generated$case == fitted) c("log_sd_r", "logit_phi_r") else character()
  if (generated$case != "rw" && curve(generated$case) == curve(fitted)) parameters <- c(parameters, "log_sr_alpha", "log_sr_beta")
  generated$truth[intersect(parameters, names(generated$truth))]
}

rec_cross_simulation <- function(case, stage, seed, sigma = .35, contrast = "natural") {
  set.seed(seed)
  obs <- rec_inputs(if (stage == "isolated") 2000:2079 else 2000:2049)
  if (stage == "full" && contrast == "wider") obs$catch$f_pressure <- 1.5 * sin(2 * pi * (obs$catch$year - 2000) / 25)
  d <- do.call(prepare_tam, c(list(data = obs), rec_cross_settings(rec_cross_formula(case), obs)))
  p <- make_par(d)
  p$log_r0 <- log(1000)
  p$log_sd_r <- log(sigma)
  if (!is.null(p$logit_phi_r)) p$logit_phi_r <- qlogis(.6)
  if (stage == "isolated") {
    S <- exp(seq(log(200), log(5000), length.out = length(d$years)) + rnorm(length(d$years), 0, .1))
    if (contrast == "narrow") S <- exp(log(1000) + rnorm(length(S), 0, .05))
    if (!is.null(d$rec$curve)) {
      p$log_sr_alpha <- log(2)
      p$log_sr_beta <- log(.001)
    }
    mean <- rep(0, length(d$years))
    if (!is.null(d$rec$curve)) mean[d$rec$eligible] <- tinyAM:::.rec_log_curve(
      log(S[d$rec$eligible - 1L]), p, d$rec$curve$type)
    log_R <- numeric(length(d$years))
    log_R[1L] <- p$log_r0
    u <- numeric(length(log_R))
    for (i in d$rec$eligible) {
      expected <- switch(d$rec$type, rw = log_R[i - 1L], iid = mean[i],
        ar1 = mean[i] + if (i == 2L) 0 else .6 * u[i - 1L])
      sd <- if (d$rec$type == "ar1" && i == 2L) sigma / sqrt(1 - .6^2) else sigma
      log_R[i] <- rnorm(1L, expected, sd)
      u[i] <- log_R[i] - mean[i]
    }
    return(list(case = case, obs = d$obs, log_R = log_R, S = S, truth = p))
  }
  p$log_mu_f[] <- c(log(.02), log(c(.05, .12, .2, .25, .3, .32, .35) / .02),
    if (contrast == "wider") 1)
  p$log_f[] <- rep(drop(d$F_modmat %*% p$log_mu_f), length.out = length(p$log_f))
  p$log_sd_f <- log(.15)
  p$log_sd_catch[] <- log(.1)
  p$log_sd_index[] <- c(log(.15), 0)
  p$log_q[] <- c(log(.4), log(.6 / .4))
  if (!is.null(d$rec$curve)) p <- tinyAM:::.initialize_rec_curve(p, d)
  simulated <- nll_fun(p, d, simulate = TRUE)
  shared <- intersect(names(simulated), names(p))
  p[shared] <- simulated[shared]
  report <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$report()
  list(case = case, obs = rec_observations(simulated, d), log_R = log(report$recruitment),
    S = report$ssb, truth = p)
}

rec_cross_isolated <- function(generated, fitted) {
  d <- do.call(prepare_tam, c(list(data = generated$obs), rec_cross_settings(rec_cross_formula(fitted), generated$obs)))
  p <- list(log_sd_r = log(.5))
  if (d$rec$type == "ar1") p$logit_phi_r <- qlogis(.5)
  if (!is.null(d$rec$curve)) {
    S <- generated$S[1L]
    p$log_sr_beta <- -log(S)
    p$log_sr_alpha <- generated$log_R[1L] - log(S) + if (d$rec$curve$type == "bh") log(2) else 1
  }
  obj <- RTMB::MakeADFun(function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mean <- rep(0, length(d$years))
    if (!is.null(d$rec$curve)) mean[d$rec$eligible] <- tinyAM:::.rec_log_curve(
      log(generated$S[d$rec$eligible - 1L]), p, d$rec$curve$type)
    RTMB::ADREPORT(mean)
    tinyAM:::.rec_nll(generated$log_R, mean, p, d)
  }, p, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  obj$fn(opt$par)
  sdr <- tryCatch(RTMB::sdreport(obj), error = function(e) NULL)
  obj$fn(opt$par)
  estimates <- obj$env$parList()
  mock_fit <- list(dat = d, sdrep = sdr,
    rep = list(rec_log_parent = log(generated$S[d$rec$eligible - 1L])))
  list(diagnostics = rec_diagnostics(obj, opt, sdr),
    advisories = tinyAM:::.rec_fit_advisories(mock_fit),
    parameters = rec_parameter_table(estimates,
      if (is.null(sdr)) list() else as.list(sdr, "Std. Error"),
      rec_cross_parameters(generated, fitted)),
    estimates = estimates)
}

rec_cross_full <- function(generated, fitted, baseline) {
  fit <- if (fitted == "rw") baseline else do.call(fit_tam,
    c(list(data = generated$obs, start_par = tinyAM:::.tam_parameter_summary(baseline, "Estimate"),
      silent = TRUE, add_osa_res = FALSE), rec_cross_settings(rec_cross_formula(fitted), generated$obs)))
  checks <- check_tam(fit, detailed = TRUE)
  estimates <- tinyAM:::.tam_parameter_summary(fit, "Estimate")
  list(diagnostics = rec_diagnostics(fit$obj, fit$opt, fit$sdrep),
    advisories = checks$advisories, numerical = checks$numerical,
    curvature = checks$curvature,
    parameters = rec_parameter_table(estimates, as.list(fit$sdrep, "Std. Error"),
      rec_cross_parameters(generated, fitted)),
    estimates = estimates, curve = fit$rec$curve, pairs = fit$rec$pairs, fit = fit)
}

rec_cross_dataset <- function(stage, truth, i, save_fit = FALSE, sigma = .35,
                              contrast = "natural", fitted_cases = rec_cross_cases, seed = NULL) {
  if (is.null(seed)) seed <- 840000L + match(truth, rec_cross_cases) * 1000L + i + if (stage == "full") 100000L else 0L
  generation_warnings <- character()
  generated <- tryCatch(withCallingHandlers(rec_cross_simulation(truth, stage, seed, sigma, contrast),
    warning = function(w) { generation_warnings <<- c(generation_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) list(error = conditionMessage(e)))
  baseline_warnings <- character()
  baseline <- if (stage == "full" && is.null(generated$error)) tryCatch(withCallingHandlers(
    do.call(fit_tam, c(list(data = generated$obs, silent = TRUE, add_osa_res = FALSE), rec_cross_settings(~ rw(year), generated$obs))),
    warning = function(w) { baseline_warnings <<- c(baseline_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) list(error = conditionMessage(e))) else NULL
  rows <- list()
  for (fitted in fitted_cases) {
    warnings <- c(generation_warnings, baseline_warnings)
    began <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers({
      if (!is.null(generated$error)) stop(generated$error)
      if (!is.null(baseline$error)) stop(baseline$error)
      if (stage == "isolated") rec_cross_isolated(generated, fitted) else rec_cross_full(generated, fitted, baseline)
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(error = conditionMessage(e)))
    diagnostics <- if (is.null(result$error)) result$diagnostics else data.frame(optimizer = NA_integer_,
      message = result$error, objective = NA_real_, gradient = NA_real_, pdHess = FALSE, success = FALSE)
    issues <- result$advisories$issue
    result$attempt <- cbind(data.frame(stage, truth, fitted, replicate = i, seed), diagnostics,
      curve_warning = any(issues %in% c("stock_recruit_support", "stock_recruit_uncertainty")),
      persistence_warning = any(issues %in% c("recruitment_AR1_uncertainty", "high_process_correlation")),
      parameter_correlation = "parameter_correlation" %in% issues,
      warm_start_converged = if (is.null(baseline$is_converged)) NA else baseline$is_converged,
      elapsed = proc.time()[["elapsed"]] - began, warnings = paste(unique(warnings), collapse = " | "))
    if (save_fit && !is.null(result$fit)) saveRDS(result$fit, file.path("analysis/formula_expansion/results",
      paste0("recruitment_cross_", truth, "_", fitted, ".rds")))
    result$fit <- NULL
    rows[[fitted]] <- result
  }
  list(results = rows, generation = generated[c("log_R", "S", "truth", "error")],
    baseline = if (!is.null(baseline$obj)) rec_diagnostics(baseline$obj, baseline$opt, baseline$sdrep) else baseline)
}
