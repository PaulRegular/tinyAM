# Run from the repository root: Rscript analysis/formula_expansion/m_random_slopes.R [100] [30] [3]
# Dashboard rendering requires Pandoc, as with vis_tam().
args <- commandArgs(trailingOnly = TRUE)
n_isolated <- if (length(args)) as.integer(args[1L]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2L]) else 30L
workers <- if (length(args) > 2L) as.integer(args[3L]) else 3L
stopifnot(n_isolated > 0, n_full > 0, workers > 0)
pkgload::load_all(".", quiet = TRUE)
repo <- normalizePath(".", winslash = "/")
output <- file.path(repo, "analysis/formula_expansion/results/m_random_slopes")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
scenarios <- list(
  informative = list(groups = 8L, persistence = .6, persistent_x = FALSE),
  four_groups = list(groups = 4L, persistence = .6, persistent_x = FALSE),
  persistent = list(groups = 4L, persistence = .9, persistent_x = TRUE))

make_case <- function(config, seed) {
  set.seed(seed)
  years <- 2000:2039
  ages <- seq.int(2L, 2L * config$groups + 2L)
  grid <- expand.grid(year = years, age = ages)
  if (config$persistent_x) {
    x <- as.numeric(stats::arima.sim(list(ar = .95), n = length(years)))
  } else x <- stats::rnorm(length(years))
  x <- as.numeric(scale(x)) # Explicitly standardize study covariates, never user data.
  covariates <- grid
  covariates$capelin <- x[match(grid$year, years)]
  covariates$M_age_group <- cut_ages(grid$age, seq.int(min(ages), max(ages), by = 2L))
  data <- list(
    catch = transform(grid, obs = 1, relative_sd = .08, unit = "thousand fish"),
    index = rbind(transform(grid, obs = 1, survey = "early", samp_time = .25,
      relative_sd = .08, unit = "thousand fish"),
      transform(grid, obs = 1, survey = "late", samp_time = .75,
        relative_sd = .08, unit = "thousand fish")),
    weight = transform(grid, obs = .08 * age^1.5, unit = "kg"),
    maturity = transform(grid, obs = plogis(age - 5)), pop_covariates = covariates)
  settings <- list(years = years, ages = ages,
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "ar1", mu_supplied = ~ I(.3), first_dev_year = min(years),
      age_breaks = seq.int(min(ages), max(ages), by = 2L),
      mu_form = ~ 0 + capelin + iid(M_age_group, by = capelin)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd),
    index_settings = list(q_form = ~ survey, sd_form = ~ 0, sd_supplied = ~ relative_sd))
  dat <- do.call(prepare_tam, c(list(data = data), settings))
  term <- dat$M_terms[[1L]]
  p <- make_par(dat)
  p$mu_m[] <- -.25
  p[[term$sd_parameter]][] <- log(.15)
  p$log_sd_m[] <- log(.06)
  p$logit_phi_m[] <- qlogis(c(.25, config$persistence))
  list(dat = dat, par = p, term = term, settings = settings, config = config)
}

intervals <- function(fixed, covariance, case) {
  specifications <- list(
    beta = list(name = "mu_m", index = 1L, truth = -.25, link = identity),
    slope_sd = list(name = case$term$sd_parameter, index = 1L, truth = .15, link = exp),
    residual_sd = list(name = "log_sd_m", index = 1L, truth = .06, link = exp),
    phi_age = list(name = "logit_phi_m", index = 1L, truth = .25, link = plogis),
    phi_year = list(name = "logit_phi_m", index = 2L, truth = case$config$persistence, link = plogis))
  do.call(rbind, lapply(names(specifications), function(parameter) {
    s <- specifications[[parameter]]
    i <- which(names(fixed) == s$name)[s$index]
    variance <- if (length(i) && !is.na(i) && is.matrix(covariance)) covariance[i, i] else NA_real_
    se <- if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
    value <- if (length(fixed) && length(i) && !is.na(i)) fixed[i] else NA_real_
    data.frame(parameter = parameter, truth = s$truth, estimate = s$link(value),
      lower = s$link(value - qnorm(.975) * se), upper = s$link(value + qnorm(.975) * se))
  }))
}

starting_values <- function(p, term) {
  p$mu_m[] <- 0
  p[[term$parameter]][] <- 0
  p[[term$sd_parameter]][] <- log(.1)
  p$log_sd_m[] <- log(.1)
  p$logit_phi_m[] <- qlogis(c(.3, .5))
  p
}

isolated_trial <- function(case) {
  dat <- case$dat
  p <- case$par
  term <- case$term
  ny <- length(dat$years)
  rows <- dat$pop_data$age %in% dat$M_settings$age_block_start
  z <- tinyAM:::.q_term_design(term, nrow(dat$pop_data))[rows, , drop = FALSE]
  x <- dat$pop_data$capelin[rows]
  nb <- ncol(z)
  phi <- plogis(p$logit_phi_m)
  covariance <- .06^2 / prod(1 - phi^2) * kronecker(
    phi[1]^abs(outer(seq_len(nb), seq_len(nb), `-`)),
    phi[2]^abs(outer(seq_len(ny), seq_len(ny), `-`)))
  b <- rnorm(nb, 0, .15)
  residual <- drop(t(chol(covariance)) %*% rnorm(nrow(z)))
  truth <- -.25 * x + drop(z %*% b) + residual
  y <- log(.3) + truth
  p <- starting_values(p[c("mu_m", term$parameter, term$sd_parameter, "log_sd_m", "logit_phi_m")], term)
  obj <- RTMB::MakeADFun(function(p) {
    effect <- tinyAM:::.mean_effects(p, dat, "M")
    eta <- matrix(y - log(.3) - p$mu_m * x - effect$contribution[rows], ny, nb)
    effect$nll - dprocess_ar1(eta, phi = plogis(p$logit_phi_m), sd = exp(p$log_sd_m))
  }, p, random = term$parameter, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  obj$fn(opt$par)
  sdrep <- RTMB::sdreport(obj, getReportCovariance = FALSE)
  gradient <- max(abs(obj$gr(opt$par)))
  numerical <- opt$convergence == 0L && isTRUE(sdrep$pdHess) && is.finite(gradient) && gradient <= .01
  fitted <- obj$env$parList()
  estimated <- fitted[[term$parameter]]
  list(parameters = intervals(sdrep$par.fixed, sdrep$cov.fixed, case),
    diagnostics = data.frame(is_converged = numerical, optimizer = opt$convergence,
      message = opt$message, objective = opt$objective, gradient = gradient, pdHess = sdrep$pdHess),
    slopes = data.frame(group = seq_along(b), truth = b, estimate = estimated),
    mean = data.frame(truth = -.25 * x + drop(z %*% b),
      estimate = as.numeric(fitted$mu_m) * x + drop(z %*% estimated)))
}

full_trial <- function(case) {
  dat <- case$dat
  p <- case$par
  term <- case$term
  p$log_r0 <- log(1000)
  p$log_sd_r[] <- log(.12)
  p$log_sd_f[] <- log(.08)
  p$log_mu_f[] <- c(log(.15), .04 * seq_len(length(dat$ages) - 1L))
  p$log_f[] <- matrix(drop(dat$F_modmat %*% p$log_mu_f), length(dat$years), length(dat$ages))
  p$log_q[] <- c(log(.7), log(.9 / .7))
  simulated <- nll_fun(p, dat, simulate = TRUE)
  keep <- intersect(names(p), names(simulated))
  p[keep] <- simulated[keep]
  data <- dat$obs
  values <- split(exp(simulated$log_obs), dat$obs_map$type)
  data$catch$obs <- values$catch
  data$index$obs <- values$index
  start <- starting_values(p, term)
  start$log_sd_r[] <- start$log_sd_f[] <- log(.1)
  fit <- do.call(fit_tam, c(list(data = data, silent = TRUE, start_par = start), case$settings))
  sdrep <- if (is.list(fit$sdrep)) fit$sdrep else list()
  levels <- fit$formula_effects$levels[[term$id]]
  truth <- as.numeric(p[[term$parameter]])
  list(parameters = intervals(sdrep$par.fixed, sdrep$cov.fixed, case),
    diagnostics = data.frame(is_converged = fit$is_converged, optimizer = fit$opt$convergence,
      message = fit$opt$message, objective = fit$opt$objective,
      gradient = fit$diagnostics$max_gradient, pdHess = isTRUE(sdrep$pdHess)),
    slopes = data.frame(group = seq_along(truth), truth = truth, estimate = levels$est,
      lower = levels$lwr, upper = levels$upr),
    mortality = data.frame(truth = as.numeric(p$log_m),
      estimate = as.numeric(log(fit$rep$M[, as.character(dat$M_settings$age_block_start)]))),
    example = list(fit = fit, truth = p))
}

attempt <- function(scenario, config, kind, replicate) {
  file <- file.path(output, paste(scenario, kind, replicate, "rds", sep = "."))
  # Checkpoints include failed fits; never silently retry an unsuccessful result.
  if (file.exists(file)) return(readRDS(file))
  warnings <- messages <- character()
  seed <- 20261019L + match(scenario, names(scenarios)) * 10000L +
    (kind == "full") * 1000L + replicate
  began <- proc.time()[3L]
  result <- withCallingHandlers(tryCatch({
    case <- make_case(config, seed)
    if (kind == "full") full_trial(case) else isolated_trial(case)
  }, error = function(e) list(error = conditionMessage(e))),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
  message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") })
  result$run <- data.frame(scenario = scenario, kind = kind, replicate = replicate, seed = seed,
    elapsed_seconds = unname(proc.time()[3L] - began), error = if (is.null(result$error)) "" else result$error,
    warnings = paste(warnings, collapse = " | "), messages = paste(messages, collapse = " | "))
  if (!is.null(result$example)) {
    if (replicate == 1L) saveRDS(result$example, file.path(output, paste0(scenario, "_example.rds")))
    result$example <- NULL
  }
  saveRDS(result, file)
  result
}

old_plan <- future::plan()
future::plan(future::multisession, workers = workers)
jobs <- rbind(expand.grid(scenario = names(scenarios), kind = "isolated", replicate = seq_len(n_isolated)),
  expand.grid(scenario = names(scenarios), kind = "full", replicate = seq_len(n_full)))
results <- furrr::future_map(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  scenario <- as.character(job$scenario)
  checkpoint <- file.path(output, paste(scenario, job$kind, job$replicate, "rds", sep = "."))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  pkgload::load_all(repo, quiet = TRUE)
  attempt(scenario, scenarios[[scenario]], as.character(job$kind), job$replicate)
}, .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
future::plan(old_plan)
runs <- do.call(rbind, lapply(results, function(x) cbind(x$run,
  if (is.null(x$diagnostics)) data.frame(is_converged = FALSE, optimizer = NA, message = "",
    objective = NA, gradient = NA, pdHess = FALSE) else x$diagnostics)))
parameters <- do.call(rbind, lapply(results, function(x) {
  if (is.null(x$parameters)) return(NULL)
  cbind(x$run[c("scenario", "kind", "replicate")], x$parameters,
    is_converged = x$diagnostics$is_converged, pdHess = x$diagnostics$pdHess)
}))
saveRDS(list(runs = runs, parameters = parameters, results = results,
  scenarios = scenarios, n_isolated = n_isolated, n_full = n_full,
  revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  script_md5 = unname(tools::md5sum(file.path(repo, "analysis/formula_expansion/m_random_slopes.R"))),
  session = sessionInfo()),
  file.path(output, "recovery.rds"))
write.csv(runs, file.path(output, "attempts.csv"), row.names = FALSE)
write.csv(parameters, file.path(output, "parameters.csv"), row.names = FALSE)
print(aggregate(is_converged ~ scenario + kind, runs, mean))

validate_workflows <- function(example_file, output) {
  example <- readRDS(example_file)
  fit <- example$fit
  original <- serialize(fit, NULL)
  checks <- list()
  step <- function(name, expression) {
    warnings <- character()
    error <- ""
    value <- withCallingHandlers(tryCatch(force(expression), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = !nzchar(error),
      error = error, warnings = paste(warnings, collapse = " | "))
    value
  }
  step("tidying and diagnostics", {
    check <- check_tam(fit)
    stopifnot(all(c("M_random_slopes_AR1", "M_random_slope_groups") %in% check$advisories$issue))
    tabs <- tidy_tam(model_list = list(Slopes = fit))
    effects <- tabs$formula_effects$levels[[names(fit$dat$M_terms)[1L]]]
    stopifnot(nrow(effects) == 4L, all(is.finite(effects$est)),
      all(is.finite(effects$lwr)), all(effects$se_scale == "reported"))
    TRUE
  })
  step("unevaluated partial update", {
    call <- update(fit, M_settings = list(mean_ages = 5:10), evaluate = FALSE)
    stopifnot(identical(call$M_settings$mu_form, fit$dat$M_settings$mu_form))
    TRUE
  })
  projected <- step("warm start and projections", {
    p <- update(fit, proj_settings = list(n_proj = 2, n_mean = 2, F_mult = 1),
      start_par = as.list(fit$sdrep, "Estimate"), silent = TRUE)
    state <- names(p$dat$M_terms)[1L]
    parameters <- as.list(p$sdrep, "Estimate")
    stopifnot(length(parameters[[p$dat$M_terms[[state]]$parameter]]) == 4L,
      nrow(parameters$log_m) == 42L,
      all(p$rep$M[, 1L] == p$rep$M[, 2L]))
    d <- p$dat$pop_data
    stopifnot(identical(d$capelin[d$year == 2041], d$capelin[d$year == 2039]))
    p
  })
  retro <- step("retrospective fold", {
    r <- fit_retro(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
    parameters <- as.list(r$fits[[1]][["2038"]]$sdrep, "Estimate")
    stopifnot(length(r$fits[[1]]) == 2L,
      nrow(parameters$log_m) == 39L,
      length(parameters[[fit$dat$M_terms[[1L]]$parameter]]) == 4L)
    r
  })
  hindcast <- step("hindcast fold", {
    h <- fit_hindcast(fit, folds = 1, start_from_fit = TRUE, proj_settings = list(F_mult = 1), progress = FALSE)
    stopifnot(nrow(h$fold_diagnostics) == 2L)
    scores <- subset(h$scores, type == "index" & level == "overall" & is.na(horizon))
    if (isTRUE(h$fits[[1]][["2038"]]$is_converged)) {
      stopifnot(sum(h$fits[[1]][["2038"]]$dat$is_proj) == 1L, is.finite(scores$RMSE))
    }
    h
  })
  if (!is.null(hindcast)) step("hindcast numerical convergence", {
    stopifnot(all(hindcast$fold_diagnostics$is_converged))
    TRUE
  })
  step("dashboard", {
    models <- list(Slopes = fit)
    if (inherits(projected, "tam_fit")) models$Projection <- projected
    file <- file.path(output, "dashboard.html")
    vis_tam(model_list = models, output_file = file, open_file = FALSE)
    stopifnot(file.exists(file), file.info(file)$size > 0)
    TRUE
  })
  step("original fit unchanged", {
    stopifnot(identical(original, serialize(fit, NULL)))
    TRUE
  })
  checks <- do.call(rbind, checks)
  write.csv(checks, file.path(output, "workflow_checks.csv"), row.names = FALSE)
  validation <- list(checks = checks,
    retrospective = if (!is.null(retro)) retro$fold_diagnostics else NULL,
    hindcast = if (!is.null(hindcast)) hindcast$fold_diagnostics else NULL)
  saveRDS(validation, file.path(output, "workflow_checks.rds"))
  saveRDS(validation, file.path(output, paste0("workflow_checks_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")))
  checks
}
print(validate_workflows(file.path(output, "four_groups_example.rds"), output))
