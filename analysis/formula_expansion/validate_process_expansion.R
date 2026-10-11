# Gaussian recovery and full population recovery ----
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
n_isolated <- if (length(args)) as.integer(args[1]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2]) else 30L
root <- "analysis/formula_expansion/results/process_expansion"
dir.create(root, recursive = TRUE, showWarnings = FALSE)

capture_attempt <- function(fun) {
  warnings <- character()
  started <- proc.time()[3]
  result <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  result$warnings <- warnings
  result$elapsed <- unname(proc.time()[3] - started)
  result
}

isolated <- function(component, grouped, seed) {
  set.seed(seed)
  design <- if (grouped) model.matrix(~ factor(rep(c("young", "old"), each = 4))) else matrix(1, 8, 1)
  beta <- if (grouped) c(log(.20), log(.08 / .20)) else log(.15)
  sd <- drop(exp(design %*% beta))
  x <- matrix(0, 81, 8)
  rho <- .6
  process <- if (component == "F") "cor_rw" else "iid"
  x <- if (process == "cor_rw") tinyAM:::.rprocess_cor_rw(x, sd, rho) else
    tinyAM:::.rprocess_scaled(x, sd, "iid")
  p <- list(beta = rep(log(.12), ncol(design)))
  if (component == "F") p$atanh_rho <- 0
  obj <- RTMB::MakeADFun(function(p) {
    scale <- exp(drop(design %*% p$beta))
    if (process == "cor_rw") -tinyAM:::.dprocess_cor_rw(x, scale, tanh(p$atanh_rho)) else
      -tinyAM:::.dprocess_scaled(x, scale, "iid")
  }, p, silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
  obj$fn(opt$par)
  sdr <- RTMB::sdreport(obj)
  estimates <- sdr$par.fixed
  se <- sqrt(diag(sdr$cov.fixed))
  truth <- c(beta, if (component == "F") atanh(rho))
  list(numerical = opt$convergence == 0 && sdr$pdHess && max(abs(obj$gr(opt$par))) < .01,
    optimizer = opt$convergence, pdHess = sdr$pdHess, gradient = max(abs(obj$gr(opt$par))),
    recovery = data.frame(parameter = c(paste0("beta", seq_along(beta)), if (component == "F") "atanh_rho"),
      truth = truth, estimate = unname(estimates), se = unname(se),
      covered = abs(estimates - truth) <= qnorm(.975) * se))
}

observations <- function() {
  grid <- expand.grid(year = 2000:2029, age = 1:8)
  index <- rbind(transform(grid, survey = "early", samp_time = .25),
    transform(grid, survey = "late", samp_time = .75))
  index$survey <- factor(index$survey)
  index$obs <- 1
  index$relative_sd <- .08
  list(catch = transform(grid, obs = 1, relative_sd = .08), index = index,
    weight = transform(grid, obs = .2 * age^1.5, age_group = factor(ifelse(age <= 4, "young", "old"))),
    maturity = transform(grid, obs = plogis(age - 4)))
}

full <- function(component, seed, keep_fit = FALSE) {
  set.seed(seed)
  settings <- list(N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "rw"),
    M_settings = list(process = "off", mu_supplied = ~ I(.25)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd, fill_missing = FALSE),
    index_settings = list(q_form = ~ survey, sd_form = ~ 0, sd_supplied = ~ relative_sd, fill_missing = FALSE),
    ssb_settings = list(spawn_time = .25))
  settings[[paste0(component, "_settings")]]$process <- if (component == "F") "cor_rw" else "iid"
  settings[[paste0(component, "_settings")]]$sd_form <- ~ age_group
  if (component == "M") {
    settings$M_settings$age_breaks <- 1:8
    settings$M_settings$first_dev_year <- 2000
  }
  dat <- do.call(prepare_tam, c(list(data = observations()), settings))
  p <- make_par(dat)
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.15)
  if (!is.null(p$log_sd_f)) p$log_sd_f <- log(.1)
  p$log_f[] <- rep(log(.05 + .025 * (1:8)), each = 30)
  p$log_q[] <- c(log(.7), log(.9 / .7))
  name <- dat$process_sd[[component]]$parameter
  truth <- c(log(.2), log(.1 / .2))
  p[[name]][] <- truth
  if (component == "F") p$atanh_rho_f <- atanh(.6)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  p[intersect(names(p), names(simulated))] <- simulated[intersect(names(p), names(simulated))]
  obs <- dat$obs
  values <- split(exp(simulated$log_obs), dat$obs_map$type)
  obs$catch$obs <- values$catch
  obs$index$obs <- values$index
  start <- make_par(dat)
  start$log_r0 <- log(8e5)
  start$log_f[] <- log(.2)
  start[[name]][] <- c(log(.15), 0)
  start$log_q[] <- c(log(.7), 0)
  fit <- do.call(fit_tam, c(list(data = obs, silent = TRUE, start_par = start), settings))
  tab <- fit$fixed_par[fit$fixed_par$par == name, ]
  result <- list(numerical = fit$is_converged, optimizer = fit$opt$convergence,
    pdHess = if (is.null(fit$sdrep)) FALSE else fit$sdrep$pdHess,
    gradient = max(abs(fit$gradient)), advisories = fit$diagnostics$advisories,
    recovery = data.frame(parameter = tab$coef, truth = truth, estimate = tab$est, se = tab$se,
      covered = tab$lwr <= truth & tab$upr >= truth))
  if (component == "F") {
    tab <- fit$fixed_par[fit$fixed_par$par == "rho_f", ]
    result$recovery <- rbind(result$recovery, data.frame(parameter = "rho_f", truth = .6,
      estimate = tab$est, se = tab$se, covered = tab$lwr <= .6 & tab$upr >= .6))
  }
  if (keep_fit) saveRDS(fit, file.path(root, paste0("example_", component, ".rds")))
  result
}

attempts <- list()
for (component in c("F", "N", "M")) for (grouped in c(FALSE, TRUE)) for (i in seq_len(n_isolated)) {
  key <- paste("isolated", component, if (grouped) "grouped" else "common", i, sep = "_")
  attempts[[key]] <- capture_attempt(function() isolated(component, grouped, 18000 + length(attempts)))
}
saveRDS(attempts, file.path(root, "recovery.rds"))
for (component in c("F", "N", "M")) for (i in seq_len(n_full)) {
  key <- paste("full", component, "grouped", i, sep = "_")
  cli::cli_inform("Recovery: {key}")
  attempts[[key]] <- capture_attempt(function() full(component, 28000 + length(attempts), i == 1L))
  saveRDS(attempts, file.path(root, "recovery.rds"))
}
rows <- lapply(names(attempts), function(key) {
  a <- attempts[[key]]
  data.frame(case = key, numerical = isTRUE(a$numerical),
    optimizer = if (is.null(a$optimizer)) NA else a$optimizer,
    pdHess = isTRUE(a$pdHess), gradient = if (is.null(a$gradient)) NA else a$gradient,
    error = if (is.null(a$error)) "" else a$error, warnings = paste(a$warnings, collapse = " | "), elapsed = a$elapsed)
})
utils::write.csv(do.call(rbind, rows), file.path(root, "attempts.csv"), row.names = FALSE)
cli::cli_inform("All recovery attempts retained in {root}.")
