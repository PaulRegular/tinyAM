sd_formula <- function(type, grouped = FALSE) {
  baseline <- if (grouped) "survey" else "1"
  term <- if (grouped) sprintf("%s(age, by = survey)", type) else sprintf("%s(age)", type)
  stats::as.formula(paste("~", baseline, "+", term))
}

sd_interval <- function(estimate, se, truth, parameter, link = "log") {
  transform <- if (link == "log") exp else plogis
  data.frame(parameter = parameter, truth = truth, estimate = transform(estimate),
    lower = transform(estimate - qnorm(.975) * se),
    upper = transform(estimate + qnorm(.975) * se))
}

sd_diagnostics <- function(obj, opt, sdr) {
  gradient <- max(abs(obj$gr(opt$par)))
  pdHess <- isTRUE(sdr$pdHess)
  data.frame(optimizer = opt$convergence, message = opt$message,
    objective = opt$objective, gradient = gradient, pdHess = pdHess,
    success = opt$convergence == 0 && pdHess && is.finite(gradient) && gradient < .01)
}

sd_simple_recovery <- function(type, repetitions, seed) {
  set.seed(seed)
  d <- expand.grid(age = 1:15, year = seq_len(repetitions))
  d$obs <- 1
  d$is_proj <- FALSE
  compiled <- tinyAM:::.parse_sd_formula(sd_formula(type), d, "catch")
  dat <- list(obs = list(catch = d), sd_catch_terms = compiled$terms)
  p <- tinyAM:::.q_term_parameters(compiled$terms)
  term <- compiled$terms[[1]]
  p[[term$sd_parameter]][] <- log(.35)
  if (type == "ar1") p[[term$phi_parameter]][] <- qlogis(.65)
  truth <- tinyAM:::.sd_effects(p, dat, "catch", TRUE)$contribution
  y <- rnorm(nrow(d), 0, exp(log(.25) + truth))
  p[[term$parameter]][] <- 0
  p[[term$sd_parameter]][] <- log(.2)
  p$baseline <- log(.25)
  obj <- RTMB::MakeADFun(function(p) {
    e <- tinyAM:::.sd_effects(p, dat, "catch")
    e$nll - sum(RTMB::dnorm(y, 0, exp(p$baseline + e$contribution), log = TRUE))
  }, p, random = term$parameter, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1500))
  obj$fn(opt$par)
  sdr <- tryCatch(RTMB::sdreport(obj, getReportCovariance = FALSE), error = function(e) NULL)
  se <- if (is.null(sdr)) rep(NA_real_, length(opt$par)) else sqrt(diag(sdr$cov.fixed))
  names(se) <- names(opt$par)
  table <- rbind(sd_interval(opt$par[term$sd_parameter], se[term$sd_parameter], .35, "process SD"),
    sd_interval(opt$par["baseline"], se["baseline"], .25, "baseline SD"))
  if (type == "ar1") table <- rbind(table, sd_interval(opt$par[term$phi_parameter],
    se[term$phi_parameter], .65, "correlation", "logit"))
  obj$fn(opt$par)
  fitted <- obj$env$parList()
  predicted <- exp(fitted$baseline + tinyAM:::.sd_effects(fitted, dat, "catch")$contribution)
  list(diagnostics = sd_diagnostics(obj, opt, sdr), parameters = table,
    curve = data.frame(age = d$age[1:15], truth = exp(log(.25) + truth[1:15]),
      estimate = predicted[1:15]), log_sd_rmse = sqrt(mean((log(predicted) - log(.25) - truth)^2)))
}

sd_tail_recovery <- function(type, seed) {
  set.seed(seed)
  d <- expand.grid(age = 1:15, year = 1:30)
  d$obs <- 1
  d$is_proj <- FALSE
  truth <- exp(log(.18) + .025 * (d$age - 8)^2)
  y <- rnorm(nrow(d), 0, truth)
  form <- switch(type, common = ~ 1, quadratic = ~ age + I(age^2), sd_formula(type))
  compiled <- tinyAM:::.parse_sd_formula(form, d, "catch")
  dat <- list(obs = list(catch = d), sd_catch_terms = compiled$terms)
  p <- c(list(beta = c(log(.3), rep(0, ncol(compiled$matrix) - 1L))),
    tinyAM:::.q_term_parameters(compiled$terms))
  obj <- RTMB::MakeADFun(function(p) {
    e <- tinyAM:::.sd_effects(p, dat, "catch")
    e$nll - sum(RTMB::dnorm(y, 0, exp(drop(compiled$matrix %*% p$beta) + e$contribution), log = TRUE))
  }, p, random = unlist(lapply(compiled$terms, `[[`, "parameter")), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1500))
  obj$fn(opt$par)
  sdr <- tryCatch(RTMB::sdreport(obj, getReportCovariance = FALSE), error = function(e) NULL)
  obj$fn(opt$par)
  fitted <- obj$env$parList()
  predicted <- exp(drop(compiled$matrix %*% fitted$beta) + tinyAM:::.sd_effects(fitted, dat, "catch")$contribution)
  list(diagnostics = sd_diagnostics(obj, opt, sdr), parameters = data.frame(),
    curve = data.frame(age = 1:15, truth = truth[1:15], estimate = predicted[1:15]),
    log_sd_rmse = sqrt(mean((log(predicted) - log(truth))^2)))
}

sd_observations <- function() {
  grid <- expand.grid(year = 2000:2029, age = 1:15)
  index <- rbind(transform(grid, survey = "early", samp_time = .25, obs = 1),
    transform(grid, survey = "late", samp_time = .75, obs = 1))
  index$survey <- factor(index$survey)
  list(catch = transform(grid, obs = 1), index = index,
    weight = transform(grid, obs = .2 * age^1.5), maturity = transform(grid, obs = plogis(age - 4)))
}

sd_full_recovery <- function(type, component, seed) {
  set.seed(seed)
  settings <- list(N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.25)),
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ survey, sd_form = ~ survey, fill_missing = FALSE))
  settings[[paste0(component, "_settings")]]$sd_form <- sd_formula(type, component == "index")
  dat <- do.call(prepare_tam, c(list(data = sd_observations()), settings))
  p <- make_par(dat)
  term <- dat[[paste0("sd_", component, "_terms")]][[1]]
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.12)
  p$log_sd_f <- log(.1)
  p$log_mu_f[] <- c(log(.15), .04 * 1:14)
  p$log_f[] <- matrix(drop(dat$F_modmat %*% p$log_mu_f), 30, 15)
  p$log_q[] <- c(log(.7), log(.9 / .7))
  p$log_sd_catch[] <- log(.2)
  p$log_sd_index[] <- c(log(.25), 0)
  p[[term$sd_parameter]][] <- log(.35)
  if (type == "ar1") p[[term$phi_parameter]][] <- qlogis(.65)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  p[intersect(names(p), names(simulated))] <- simulated[intersect(names(p), names(simulated))]
  truth <- exp(drop(dat[[paste0("sd_", component, "_modmat")]] %*%
    p[[paste0("log_sd_", component)]]) + tinyAM:::.sd_effects(p, dat, component)$contribution)
  obs <- dat$obs
  values <- split(exp(simulated$log_obs), dat$obs_map$type)
  obs$catch$obs <- values$catch
  obs$index$obs <- values$index
  start <- make_par(dat)
  start$log_r0 <- log(1e6)
  start$log_mu_f <- p$log_mu_f
  start$log_f[] <- matrix(drop(dat$F_modmat %*% start$log_mu_f), 30, 15)
  start$log_q <- p$log_q
  start$log_sd_r <- start$log_sd_f <- log(.1)
  start$log_sd_catch[] <- log(.2)
  start$log_sd_index[] <- c(log(.25), 0)
  fit <- do.call(fit_tam, c(list(data = obs, silent = TRUE, start_par = start, grad_tol = .01), settings))
  table <- fit$fixed_par
  select <- function(parameter, truth, label) {
    row <- table[table$par == sub("^log(it)?_", "", parameter), ]
    data.frame(parameter = label, truth = truth, estimate = row$est, lower = row$lwr, upper = row$upr)
  }
  parameters <- select(term$sd_parameter, .35, "process SD")
  if (type == "ar1") parameters <- rbind(parameters, select(term$phi_parameter, .65, "correlation"))
  predicted <- fit$rep$sd_obs[dat$obs_map$type == component]
  curves <- unique(data.frame(age = dat$obs[[component]]$age,
    group = if (component == "index") dat$obs$index$survey else "catch", truth, estimate = predicted))
  list(diagnostics = data.frame(optimizer = fit$opt$convergence, message = fit$opt$message,
      objective = fit$opt$objective, gradient = fit$diagnostics$max_gradient,
      pdHess = isTRUE(fit[["sdrep"]]$pdHess), success = fit$is_converged),
    parameters = parameters, curve = curves, fit = fit,
    log_sd_rmse = sqrt(mean((log(predicted) - log(truth))^2)))
}
