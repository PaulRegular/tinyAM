q_specs <- c("iid", "rw", "ar1", "intercept", "numeric_by", "logistic")

q_formula <- function(spec) {
  switch(spec,
    iid = ~ survey + iid(year, by = survey),
    rw = ~ survey + rw(year, by = survey),
    ar1 = ~ survey + ar1(year, by = survey),
    intercept = ~ survey + (1 | q_group),
    numeric_by = ~ survey + rw(year, by = scaled_age),
    logistic = ~ survey + logistic(age, by = survey))
}

q_grid <- function(sparse = FALSE, full = FALSE) {
  d <- expand.grid(year = if (full) 2000:2029 else 2000:2039,
    age = if (sparse) 1L else 1:8, survey = c("early", "late"))
  d$survey <- factor(d$survey)
  d$q_group <- interaction(d$year, d$survey, drop = TRUE)
  d$scaled_age <- d$age / 8
  d$obs <- 1
  d$is_proj <- FALSE
  d
}

q_truth_parameters <- function(p, dat, link) {
  baseline <- if (link == "log") log(c(.55, .8)) else qlogis(c(.55, .8))
  nm <- if (link == "log") "log_q" else "logit_q"
  p[[nm]][] <- c(baseline[1], diff(baseline))
  for (term in dat$q_terms) {
    if (term$type == "logistic") {
      p[[paste0("q_a50_", term$id)]][] <- c(3.5, 4.5)
      p[[paste0("log_q_slope_", term$id)]][] <- log(c(1.2, 1))
    } else {
      p[[term$sd_parameter]][] <- log(.18)
      if (term$type == "ar1") p[[term$phi_parameter]][] <- qlogis(.65)
    }
  }
  p
}

q_log_curve <- function(p, dat, effects = tinyAM:::.q_effects(p, dat)) {
  link <- dat$index_settings$q_link
  predictor <- drop(dat$q_modmat %*% p[[if (link == "log") "log_q" else "logit_q"]]) +
    effects$contribution
  if (link == "logit") predictor <- -RTMB::logspace_add(0, -predictor)
  predictor + effects$log_selectivity
}

q_parameter_intervals <- function(estimates, errors, truth) {
  names <- grep("^(log_sd_q_|logit_phi_q_|q_a50_|log_q_slope_)", names(truth), value = TRUE)
  do.call(rbind, lapply(names, function(nm) {
    transform <- if (startsWith(nm, "logit_")) plogis else
      if (startsWith(nm, "log_")) exp else identity
    data.frame(parameter = sub("^log(it)?_", "", nm), group = seq_along(truth[[nm]]),
      truth = transform(truth[[nm]]), estimate = transform(estimates[[nm]]),
      lower = transform(estimates[[nm]] - qnorm(.975) * errors[[nm]]),
      upper = transform(estimates[[nm]] + qnorm(.975) * errors[[nm]]))
  }))
}

q_capture <- function(fun, ...) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers(fun(...), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  result$warnings <- unique(warnings)
  result
}

q_isolated_recovery <- function(spec, link, sparse, seed) {
  set.seed(seed)
  d <- q_grid(sparse)
  if (sparse && spec == "logistic") {
    d <- q_grid(FALSE)
    d <- d[d$age <= 3, ]
  }
  dat <- c(list(obs = list(index = d), index_settings = list(q_link = link),
    sd_index_modmat = model.matrix(~ 1, d)), tinyAM:::.parse_q_formula(q_formula(spec), d))
  tinyAM:::.check_q_terms(dat)
  p <- tinyAM:::.q_term_parameters(dat$q_terms)
  nm <- if (link == "log") "log_q" else "logit_q"
  p[[nm]] <- setNames(numeric(ncol(dat$q_modmat)), colnames(dat$q_modmat))
  p$log_observation_sd <- log(.1)
  p <- q_truth_parameters(p, dat, link)
  simulated <- tinyAM:::.q_effects(p, dat, simulate = TRUE)
  truth <- q_log_curve(p, dat, simulated)
  y <- as.numeric(truth) + rnorm(nrow(d), 0, .1)
  p[names(simulated$parameters)] <- simulated$parameters
  original <- p
  p[names(tinyAM:::.q_term_parameters(dat$q_terms))] <- tinyAM:::.q_term_parameters(dat$q_terms)
  obj <- RTMB::MakeADFun(function(p) {
    effects <- tinyAM:::.q_effects(p, dat)
    log_q <- q_log_curve(p, dat, effects)
    RTMB::ADREPORT(log_q)
    effects$nll - sum(RTMB::dnorm(y, log_q, exp(p$log_observation_sd), log = TRUE))
  }, p, random = tinyAM:::.q_random_parameters(dat), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  obj$fn(opt$par)
  gradient <- max(abs(obj$gr(opt$par)))
  sdr <- RTMB::sdreport(obj, par.fixed = opt$par, getReportCovariance = FALSE)
  estimates <- as.list(sdr, "Estimate")
  errors <- as.list(sdr, "Std. Error")
  curves <- as.numeric(as.list(sdr, "Estimate", report = TRUE)$log_q)
  curve_se <- as.numeric(as.list(sdr, "Std. Error", report = TRUE)$log_q)
  list(parameters = q_parameter_intervals(estimates, errors, original),
    diagnostics = data.frame(optimizer = opt$convergence, message = opt$message,
      gradient = gradient, pdHess = isTRUE(sdr$pdHess),
      converged = opt$convergence == 0 && isTRUE(sdr$pdHess) && is.finite(gradient) && gradient <= .01,
      q_rmse = sqrt(mean((exp(curves) - exp(truth))^2)),
      q_coverage = mean(abs(curves - truth) <= qnorm(.975) * curve_se),
      observation_sd = exp(estimates$log_observation_sd)),
    curve = data.frame(d, truth = exp(truth), estimate = exp(curves),
      lower = exp(curves - qnorm(.975) * curve_se), upper = exp(curves + qnorm(.975) * curve_se)))
}

q_full_observations <- function(sparse = FALSE) {
  grid <- expand.grid(year = 2000:2029, age = 1:8)
  index <- q_grid(sparse, full = TRUE)
  index$samp_time <- ifelse(index$survey == "early", .25, .75)
  index$relative_sd <- .05
  list(catch = transform(grid, obs = 1, relative_sd = .08), index = index,
    weight = transform(grid, obs = .2 * age^1.5),
    maturity = transform(grid, obs = plogis(age - 4)))
}

q_full_settings <- function(spec, link) {
  list(N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.25)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd, fill_missing = FALSE),
    index_settings = list(q_form = q_formula(spec), q_link = link, sd_form = ~ 1,
      sd_supplied = NULL, fill_missing = FALSE))
}

q_full_recovery <- function(spec, link, seed, sparse = FALSE) {
  set.seed(seed)
  settings <- q_full_settings(spec, link)
  dat <- do.call(prepare_tam, c(list(data = q_full_observations(sparse)), settings))
  p <- make_par(dat)
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.12)
  p$log_sd_f <- log(.08)
  p$log_mu_f[] <- c(log(.15), .09 * seq_len(7))
  p$log_f[] <- matrix(drop(dat$F_modmat %*% p$log_mu_f), 30, 8)
  p$log_sd_index[] <- log(.1)
  p <- q_truth_parameters(p, dat, link)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  keep <- intersect(names(p), names(simulated))
  p[keep] <- simulated[keep]
  actual_q <- as.numeric(exp(q_log_curve(p, dat)))
  obs <- dat$obs
  values <- split(exp(simulated$log_obs), dat$obs_map$type)
  obs$catch$obs <- values$catch
  obs$index$obs <- values$index
  start <- p
  start[names(tinyAM:::.q_term_parameters(dat$q_terms))] <- tinyAM:::.q_term_parameters(dat$q_terms)
  fit <- do.call(fit_tam, c(list(data = obs, start_par = start, silent = TRUE,
    add_osa_res = FALSE, grad_tol = .01), settings))
  estimates <- tinyAM:::.tam_parameter_summary(fit, "Estimate")
  errors <- tinyAM:::.tam_parameter_summary(fit, "Std. Error")
  q <- fit$obs_pred$index
  list(parameters = q_parameter_intervals(estimates, errors, p),
    diagnostics = data.frame(optimizer = fit$opt$convergence, message = fit$opt$message,
      gradient = fit$diagnostics$max_gradient, pdHess = isTRUE(fit$sdrep$pdHess),
      converged = fit$is_converged, q_rmse = sqrt(mean((q$q - actual_q)^2)),
      q_coverage = mean(q$q_lwr <= actual_q & q$q_upr >= actual_q),
      observation_sd = exp(estimates$log_sd_index)),
    curve = data.frame(q[, c("year", "age", "survey")], truth = actual_q,
      estimate = q$q, lower = q$q_lwr, upper = q$q_upr), fit = fit)
}
