# Generating models and recovery ----
f_curve_inputs <- function(years = 2000:2029, ages = 1:8) {
  grid <- expand.grid(year = years, age = ages)
  index <- rbind(
    transform(grid, survey = "early", samp_time = .25),
    transform(grid, survey = "late", samp_time = .75)
  )
  index$survey <- factor(index$survey)
  list(
    catch = transform(grid, obs = 1, relative_sd = .08),
    index = transform(index, obs = 1, relative_sd = .08),
    weight = transform(grid, obs = .2 * age^1.5),
    maturity = transform(grid, obs = plogis(age - 4), spawn_time = .25)
  )
}

f_curve_settings <- function(curve, process) {
  list(
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(
      process = process,
      mu_form = if (curve == "mono") ~ mono(age) else ~ logistic(age)
    ),
    M_settings = list(process = "off", mu_supplied = ~ I(.25)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd),
    index_settings = list(q_form = ~ survey, sd_form = ~ 0,
                          sd_supplied = ~ relative_sd))
}

f_curve_truth <- function(dat, curve, stress = "none") {
  p <- make_par(dat)
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.25)
  p$log_sd_f <- log(.15)
  p$log_mu_f[] <- log(if (curve == "mono") .08 else .5)
  if (curve == "mono") {
    p$dmu_f[] <- c(.7, .5, .3, if (stress == "plateau") rep(0, 4) else c(.15, .1, .05, .05))
  } else {
    term <- dat$F_terms[[1]]
    p[[term$a50_parameter]][] <- if (stress == "limited_age") 1 else 4.5
    p[[term$slope_parameter]][] <- log(1.1)
  }
  if (dat$F_settings$process == "ar1") p$logit_phi_f[] <- qlogis(c(.3, .5))
  p$log_q[] <- c(log(.7), log(.8 / .7))
  p$log_f[] <- matrix(
    drop(dat$F_modmat %*% p$log_mu_f) +
      tinyAM:::.mean_effects(p, dat, "F")$contribution,
    length(dat$years), length(dat$ages)
  )
  p
}

f_curve_recovery <- function(estimate, error, truth) {
  names <- intersect(c("log_mu_f", "dmu_f", grep("^(F_a50_|log_F_slope_)", names(truth), value = TRUE)), names(estimate))
  do.call(rbind, lapply(names, function(name) {
    data.frame(
      parameter = paste(name, seq_along(truth[[name]]), sep = ":"),
      truth = as.numeric(truth[[name]]), estimate = as.numeric(estimate[[name]]),
      se = as.numeric(error[[name]]),
      covered = abs(estimate[[name]] - truth[[name]]) <= qnorm(.975) * error[[name]]
    )
  }))
}

f_curve_isolated <- function(curve, process, stress, seed) {
  set.seed(seed)
  dat <- do.call(prepare_tam, c(list(data = f_curve_inputs(2000:2049)), f_curve_settings(curve, process)))
  truth <- f_curve_truth(dat, curve, stress)
  residual <- if (process == "iid") {
    matrix(rnorm(50 * 8, 0, .15), 50, 8)
  } else {
    rprocess_ar1(50, 8, c(.3, .5), .15)
  }
  observed <- truth$log_f + residual
  full <- make_par(dat)
  keys <- c(
    "log_mu_f", "dmu_f", "log_sd_f", "logit_phi_f",
    grep("^(F_a50_|log_F_slope_)", names(full), value = TRUE)
  )
  p <- full[intersect(keys, names(full))]
  p$log_mu_f[] <- log(.3)
  p$log_sd_f <- log(.2)
  obj <- RTMB::MakeADFun(function(p) {
    mean <- matrix(drop(dat$F_modmat %*% p$log_mu_f) +
      tinyAM:::.mean_effects(p, dat, "F")$contribution, 50, 8)
    if (process == "iid") {
      -sum(RTMB::dnorm(observed - mean, 0, exp(p$log_sd_f), log = TRUE))
    } else {
      -dprocess_ar1(observed - mean, plogis(p$logit_phi_f), exp(p$log_sd_f))
    }
  }, p, silent = TRUE)
  lower <- rep(-Inf, length(obj$par))
  lower[names(obj$par) == "dmu_f"] <- 0
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower)
  obj$fn(opt$par)
  sdr <- RTMB::sdreport(obj)
  gradient <- obj$gr(opt$par)
  gradient[opt$par == lower & gradient > 0] <- 0
  estimates <- as.list(sdr, "Estimate")
  errors <- as.list(sdr, "Std. Error")
  fitted <- matrix(drop(dat$F_modmat %*% estimates$log_mu_f) +
    tinyAM:::.mean_effects(estimates, dat, "F")$contribution, 50, 8)
  list(
    numerical = opt$convergence == 0 && sdr$pdHess && max(abs(gradient)) < .01,
    optimizer = opt$convergence, message = opt$message, pdHess = sdr$pdHess,
    gradient = max(abs(gradient)), curve_rmse = sqrt(mean((exp(fitted) - exp(truth$log_f))^2)),
    recovery = f_curve_recovery(estimates, errors, truth)
  )
}

f_curve_full <- function(curve, process, stress, seed, keep = NULL) {
  set.seed(seed)
  settings <- f_curve_settings(curve, process)
  dat <- do.call(prepare_tam, c(list(data = f_curve_inputs()), settings))
  truth <- f_curve_truth(dat, curve, stress)
  sim <- nll_fun(truth, dat, simulate = TRUE)
  data <- dat$obs
  values <- split(exp(sim$log_obs), dat$obs_map$type)
  data$catch$obs <- values$catch
  data$index$obs <- values$index
  start <- make_par(dat)
  start$log_r0 <- log(8e5)
  start$log_mu_f[] <- log(.3)
  start$log_f[] <- log(.2)
  start$log_q[] <- c(log(.7), 0)
  fit <- do.call(fit_tam, c(list(data = data, start_par = start, silent = TRUE), settings))
  estimates <- tinyAM:::.tam_parameter_summary(fit, "Estimate")
  errors <- tinyAM:::.tam_parameter_summary(fit, "Std. Error")
  profile <- tidy_pop(fit)$mu_F
  true_mean <- as.numeric(exp(truth$log_f))
  result <- list(
    numerical = fit$is_converged, optimizer = fit$opt$convergence,
    message = fit$opt$message, pdHess = !is.null(fit$sdrep) && fit$sdrep$pdHess,
    gradient = fit$diagnostics$max_gradient,
    curve_rmse = sqrt(mean((profile$est - true_mean)^2)),
    curve_coverage = mean(profile$lwr <= true_mean & profile$upr >= true_mean),
    recovery = f_curve_recovery(estimates, errors, truth), advisories = fit$diagnostics$advisories
  )
  if (!is.null(keep)) saveRDS(fit, keep)
  result
}

# Every attempt is retained; there are no retries or model changes.
run_f_curve_validation <- function(n_isolated = 100L, n_full = 30L, n_stress = 20L,
                                   root = "analysis/formula_expansion/results/f_mean_curves") {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  attempts <- list()
  for (kind in c("isolated", "full")) {
    for (curve in c("mono", "logistic")) {
      for (process in c("iid", "ar1")) {
        for (stress in c("none", if (curve == "mono") "plateau" else "limited_age")) {
          n <- if (stress != "none") n_stress else if (kind == "isolated") n_isolated else n_full
          for (i in seq_len(n)) {
            key <- paste(kind, curve, process, stress, i, sep = "_")
            warnings <- character()
            started <- proc.time()[3]
            fun <- if (kind == "isolated") f_curve_isolated else f_curve_full
            cli::cli_inform("F curve recovery: {key}")
            call <- list(curve, process, stress, 35000 + length(attempts))
            if (kind == "full" && i == 1L) call$keep <- file.path(root, paste0(key, ".rds"))
            result <- tryCatch(withCallingHandlers(do.call(fun, call), warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }), error = function(e) list(error = conditionMessage(e)))
            result$warnings <- warnings
            result$elapsed <- unname(proc.time()[3] - started)
            result$seed <- call[[4L]]
            result$kind <- kind
            result$curve <- curve
            result$process <- process
            result$stress <- stress
            attempts[[key]] <- result
            saveRDS(attempts, file.path(root, "recovery.rds"))
          }
        }
      }
    }
  }
  rows <- lapply(names(attempts), function(key) {
    a <- attempts[[key]]
    value <- function(name, missing = NA) if (is.null(a[[name]])) missing else a[[name]]
    data.frame(
      case = key, kind = a$kind, curve = a$curve, process = a$process, stress = a$stress,
      numerical = isTRUE(a$numerical), optimizer = value("optimizer"), message = value("message", ""),
      pdHess = isTRUE(a$pdHess), gradient = value("gradient"), curve_rmse = value("curve_rmse"),
      curve_coverage = value("curve_coverage"), error = value("error", ""),
      warnings = paste(a$warnings, collapse = " | "), elapsed = a$elapsed
    )
  })
  write.csv(do.call(rbind, rows), file.path(root, "attempts.csv"), row.names = FALSE)
  invisible(attempts)
}

view_f_curve_examples <- function(root = "analysis/formula_expansion/results/f_mean_curves",
                                  open_file = TRUE) {
  paths <- list.files(root, "^full_.*none_1[.]rds$", full.names = TRUE)
  if (!length(paths)) cli::cli_abort("Run the F mean-curve validation before viewing its examples.")
  fits <- setNames(lapply(paths, readRDS), sub("^full_|_none_1[.]rds$", "", basename(paths)))
  vis_tam(
    model_list = fits, output_file = file.path(root, "dashboard.html"),
    open_file = open_file, background = c(
      "### F mean-curve recovery examples", "",
      "Four first-replicate examples use independently generated observations.",
      "They show recovery for each generating model; the report summarizes all replicates."
    ),
    render_args = list(quiet = TRUE)
  )
}

if (sys.nframe() == 0L) {
  pkgload::load_all(quiet = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  run_f_curve_validation(
    if (length(args)) as.integer(args[1]) else 100L,
    if (length(args) > 1L) as.integer(args[2]) else 30L,
    if (length(args) > 2L) as.integer(args[3]) else 20L
  )
}
