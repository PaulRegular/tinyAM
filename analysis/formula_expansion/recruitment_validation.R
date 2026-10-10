rec_cases <- c("bh_iid", "bh_ar1", "ricker_iid", "ricker_ar1", "cov_iid", "cov_rw")

rec_formula <- function(case) {
  parts <- strsplit(case, "_")[[1L]]
  curve <- if (parts[1L] == "cov") "temperature" else paste0(parts[1L], "(ssb)")
  stats::as.formula(paste("~", curve, "+", paste0(parts[2L], "(year)")))
}

rec_inputs <- function(years = 2000:2049) {
  d <- expand.grid(year = years, age = 1:8)
  catch <- transform(d, obs = 100)
  index <- rbind(transform(d, obs = 100, survey = "A", samp_time = .3),
                 transform(d, obs = 100, survey = "B", samp_time = .8))
  weight <- transform(d, obs = age^2 / 10)
  maturity <- transform(d, obs = .5 * plogis(2 * (age - 3)),
                         temperature = sin((year - min(year)) / 4))
  maturity$temperature[maturity$age != 1L] <- NA
  list(catch = catch, index = index, weight = weight, maturity = maturity)
}

rec_settings <- function(form) {
  list(N_settings = list(process = "off", init = "exp", rec_form = form),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.2)),
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ survey, sd_form = ~ survey, fill_missing = FALSE))
}

rec_diagnostics <- function(obj, opt, sdr) {
  gradient <- max(abs(obj$gr(opt$par)))
  data.frame(optimizer = opt$convergence, message = if (is.null(opt$message)) "" else opt$message,
    objective = opt$objective, gradient = gradient, pdHess = isTRUE(sdr$pdHess),
    success = opt$convergence == 0 && isTRUE(sdr$pdHess) && is.finite(gradient) && gradient < .01)
}

rec_parameter_table <- function(est, se, truth) {
  do.call(rbind, lapply(names(truth), function(name) {
    e <- est[[name]]
    s <- se[[name]]
    if (is.null(e)) return(NULL)
    if (is.null(s)) s <- rep(NA_real_, length(e))
    data.frame(parameter = paste(name, seq_along(e), sep = ":"), truth = unname(truth[[name]]),
      estimate = unname(e), lower = unname(e - qnorm(.975) * s), upper = unname(e + qnorm(.975) * s))
  }))
}

rec_observations <- function(simulated, dat) {
  obs <- dat$obs
  for (type in c("catch", "index")) obs[[type]]$obs <- exp(simulated$log_obs[dat$obs_map$type == type])
  obs
}

rec_isolated_recovery <- function(case, design, seed) {
  set.seed(seed)
  obs <- rec_inputs(2000:2079)
  d <- do.call(prepare_tam, c(list(data = obs), rec_settings(rec_formula(case))))
  p <- make_par(d)
  truth <- list(log_sd_r = log(.35))
  if (d$rec$type == "ar1") truth$logit_phi_r <- qlogis(.6)
  S <- exp(seq(log(200), log(5000), length.out = length(d$years)) + rnorm(length(d$years), 0, .1))
  if (design == "narrow SSB") S <- exp(log(1000) + rnorm(length(S), 0, .05))
  if (design == "correlated covariate") {
    d$rec$data$temperature <- as.numeric(scale(log(S)))
    d$rec$matrix <- cbind(temperature = d$rec$data$temperature)
  }
  if (is.null(d$rec$curve)) truth$rec_beta <- if (d$rec$type == "rw") .3 else c(log(1000), .3) else {
    truth$log_sr_alpha <- log(2)
    truth$log_sr_beta <- log(.001)
    if (design == "correlated covariate") truth$rec_beta <- .3
  }
  if (is.null(d$rec$curve)) {
    p$rec_beta <- setNames(truth$rec_beta, colnames(d$rec$matrix))
    mean <- tinyAM:::.rec_mean(p, d)
  } else {
    p$log_sr_alpha <- truth$log_sr_alpha
    p$log_sr_beta <- truth$log_sr_beta
    if (design == "correlated covariate") p$rec_beta <- c(temperature = .3)
    mean <- tinyAM:::.rec_mean(p, d)
    mean[d$rec$eligible] <- mean[d$rec$eligible] + tinyAM:::.rec_log_curve(
      log(S[d$rec$eligible - d$rec$curve$lag]), p, d$rec$curve$type)
  }
  u <- numeric(length(d$years))
  if (d$rec$type == "rw") u[1L] <- log(1000) - mean[1L]
  for (i in d$rec$eligible) {
    expected <- switch(d$rec$type, iid = 0, rw = u[i - 1L], ar1 = .6 * u[i - 1L])
    sd <- if (d$rec$type == "ar1" && i == d$rec$eligible[1L]) .35 / sqrt(1 - .6^2) else .35
    u[i] <- rnorm(1, expected, sd)
  }
  log_R <- mean + u
  estimate <- lapply(truth, function(x) x + .1)
  if (!is.null(estimate$rec_beta)) names(estimate$rec_beta) <- colnames(d$rec$matrix)
  obj <- RTMB::MakeADFun(function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mean <- tinyAM:::.rec_mean(p, d)
    if (!is.null(d$rec$curve)) mean[d$rec$eligible] <- mean[d$rec$eligible] +
      tinyAM:::.rec_log_curve(log(S[d$rec$eligible - d$rec$curve$lag]), p, d$rec$curve$type)
    RTMB::ADREPORT(mean)
    tinyAM:::.rec_nll(log_R, mean, p, d)
  }, estimate, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  obj$fn(opt$par)
  sdr <- tryCatch(RTMB::sdreport(obj), error = function(e) NULL)
  obj$fn(opt$par)
  fitted <- obj$env$parList()
  se <- if (is.null(sdr)) list() else as.list(sdr, "Std. Error")
  predicted <- tinyAM:::.rec_mean(fitted, d)
  if (!is.null(d$rec$curve)) predicted[d$rec$eligible] <- predicted[d$rec$eligible] +
    tinyAM:::.rec_log_curve(log(S[d$rec$eligible - d$rec$curve$lag]), fitted, d$rec$curve$type)
  prediction_se <- if (is.null(sdr)) rep(NA_real_, length(mean)) else
    as.list(sdr, "Std. Error", report = TRUE)$mean
  list(diagnostics = rec_diagnostics(obj, opt, sdr),
    parameters = rec_parameter_table(fitted, se, truth),
    curve = data.frame(year = d$years[d$rec$eligible], ssb = S[d$rec$eligible - if (is.null(d$rec$curve)) 0 else d$rec$curve$lag],
      truth = exp(mean[d$rec$eligible]), estimate = exp(predicted[d$rec$eligible]),
      lower = exp(predicted[d$rec$eligible] - qnorm(.975) * prediction_se[d$rec$eligible]),
      upper = exp(predicted[d$rec$eligible] + qnorm(.975) * prediction_se[d$rec$eligible])),
    curve_rmse = sqrt(mean((predicted[d$rec$eligible] - mean[d$rec$eligible])^2)))
}

rec_full_recovery <- function(case, seed) {
  set.seed(seed)
  obs <- rec_inputs()
  settings <- rec_settings(rec_formula(case))
  d <- do.call(prepare_tam, c(list(data = obs), settings))
  p <- make_par(d)
  p$log_r0 <- log(1000)
  p$log_mu_f[] <- c(log(.02), log(c(.05, .12, .2, .25, .3, .32, .35) / .02))
  p$log_f[] <- rep(drop(d$F_modmat %*% p$log_mu_f), length.out = length(p$log_f))
  p$log_sd_f <- log(.15)
  p$log_sd_r <- log(.35)
  p$log_sd_catch[] <- log(.1)
  p$log_sd_index[] <- c(log(.15), 0)
  p$log_q[] <- c(log(.4), log(.6 / .4))
  if (!is.null(p$logit_phi_r)) p$logit_phi_r[] <- qlogis(.6)
  if (is.null(d$rec$curve)) p$rec_beta[] <- if (d$rec$type == "rw") .3 else c(log(1000), .3) else
    p <- tinyAM:::.initialize_rec_curve(p, d)
  truth <- p[intersect(c("log_sd_r", "logit_phi_r", "log_sr_alpha", "log_sr_beta", "rec_beta"), names(p))]
  generated <- nll_fun(p, d, simulate = TRUE)
  p[names(generated)[names(generated) %in% names(p)]] <- generated[names(generated) %in% names(p)]
  report <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$report()
  obs <- rec_observations(generated, d)
  baseline_settings <- rec_settings(~ rw(year))
  baseline <- do.call(fit_tam, c(list(data = obs, silent = TRUE, add_osa_res = FALSE), baseline_settings))
  warm_start <- rec_diagnostics(baseline$obj, baseline$opt, baseline$sdrep)
  if (!baseline$is_converged) {
    warm_start$success <- FALSE
    warm_start$message <- paste("Prespecified RW warm start failed:", warm_start$message)
    return(list(diagnostics = warm_start, warm_start = warm_start))
  }
  start <- as.list(baseline$sdrep, "Estimate")
  if (is.null(d$rec$curve)) start$rec_beta <- setNames(if (d$rec$type == "rw") 0 else c(mean(log(baseline$rep$recruitment)), 0), colnames(d$rec$matrix))
  fit <- do.call(fit_tam, c(list(data = obs, start_par = start, silent = TRUE, add_osa_res = FALSE), settings))
  estimates <- tinyAM:::.tam_parameter_summary(fit, "Estimate")
  se <- as.list(fit$sdrep, "Std. Error")
  recruitment <- data.frame(year = d$years, ssb = report$ssb, truth = report$recruitment,
                            estimate = fit$rep$recruitment)
  curve <- if (is.null(d$rec$curve)) {
    data.frame(year = d$years, truth = exp(tinyAM:::.rec_mean(p, d)),
               estimate = exp(tinyAM:::.rec_mean(estimates, d)))
  } else {
    tab <- tidy_recruitment(fit)$curve
    data.frame(ssb = tab$ssb, truth = exp(tinyAM:::.rec_log_curve(log(tab$ssb), p, d$rec$curve$type)),
               estimate = tab$est, lower = tab$lwr, upper = tab$upr)
  }
  list(diagnostics = rec_diagnostics(fit$obj, fit$opt, fit$sdrep), warm_start = warm_start,
    parameters = rec_parameter_table(estimates, se, truth),
    curve = curve, recruitment = recruitment,
    curve_rmse = sqrt(mean((log(curve$estimate) - log(curve$truth))^2)),
    fit = fit)
}
