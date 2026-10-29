# Run from the repository root with R 4.4 and the pinned SAM installation:
# Rscript analysis/sam_equivalence/compare_haddock.R model.RData SAM-library
# Optional third argument: output directory (default below; ignored by Git).
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Supply the official haddock model.RData path and a SAM library directory.")
.libPaths(c(normalizePath(args[2L]), .libPaths()))
if (!requireNamespace("stockassessment", quietly = TRUE) || !requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install stockassessment at the pinned revision and pkgload before running.")
}
revision <- utils::packageDescription("stockassessment")$RemoteSha
stopifnot(identical(revision, "1cc464b80f6f2a4a38b35cce276f13e521761da3"))
pkgload::load_all(".", quiet = TRUE)
output <- if (length(args) >= 3L) args[3L] else "analysis/sam_equivalence/results/haddock"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stopifnot(unname(tools::md5sum(args[1L])) == "3ae172debf40dbc59d202335d9d08371")
warnings <- character()
tryCatch(withCallingHandlers({
  env <- new.env()
  load(args[1L], env)
  native <- Filter(function(x) is.list(x) && all(c("conf", "data", "pl") %in% names(x)), as.list(env))[[1L]]
  cf <- native$conf; src <- native$data
  stopifnot(cf$minAge == 0, cf$maxAge == 8, cf$initState == 0,
    all(cf$logNMeanAssumption == 0), all(cf$fracMixN == 0), cf$fracMixF == 0,
    all(cf$fracMixObs == 0), all(as.character(cf$obsCorStruct) == "ID"),
    all(src$propF == 0), all(src$propM == 0), cf$corFlag == 2,
    identical(as.numeric(cf$keyLogFsta[1L, ]), 0:8 + 0),
    cf$stockWeightModel == 0, cf$catchWeightModel == 0)
  years <- as.integer(src$years); ages <- 0:8
  stopifnot(identical(years, 1972:2026))
  grid <- expand.grid(year = years, age = ages)
  aux <- as.data.frame(src$aux)
  names(aux)[1:3] <- c("year", "fleet", "age")
  cells <- cbind(aux$fleet, aux$age + 1L)
  key <- function(x) factor(ifelse(x < 0, NA, paste0("K", x)))
  catch <- grid
  rows <- match(paste(grid$year, grid$age), paste(aux$year[aux$fleet == 1], aux$age[aux$fleet == 1]))
  catch$obs <- exp(src$logobs[aux$fleet == 1])[rows]
  catch$sd_key <- key(cf$keyVarObs[1L, catch$age + 1L])
  catch$unit <- "thousand fish"
  i <- which(aux$fleet > 1)
  index <- data.frame(year = aux$year[i], age = aux$age[i], obs = exp(src$logobs[i]),
    survey = ifelse(aux$fleet[i] == 2L, "IBTS Q1", "IBTS Q3+Q4"),
    samp_time = src$sampleTimes[aux$fleet[i]])
  index$q_key <- key(cf$keyLogFpar[cells[i, , drop = FALSE]])
  index$sd_key <- key(cf$keyVarObs[cells[i, , drop = FALSE]])
  index$q_power_group <- key(cf$keyQpow[cells[i, , drop = FALSE]])
  index$relative_sd <- ifelse(is.finite(src$weight[i]), 1 / sqrt(src$weight[i]), 1)
  catch_weights <- matrix(src$catchMeanWeight, nrow = dim(src$catchMeanWeight)[1L])
  stopifnot(ncol(catch_weights) == length(ages))
  if (nrow(catch_weights) == length(years) - 1L) {
    # The survey-only final year has no catch-weight input in SAM. This padding
    # satisfies tinyAM's complete grid but enters neither likelihood nor SSB.
    catch_weights <- rbind(catch_weights, catch_weights[nrow(catch_weights), ])
    cat("2026 catch weights retain 2025 solely for unused catch reporting; 2026 yield is not compared.\n")
  }
  weight <- rbind(transform(grid, obs = as.vector(src$stockMeanWeight), type = "stock", unit = "kg"),
    transform(grid, obs = as.vector(catch_weights), type = "catch", unit = "kg"))
  maturity <- transform(grid, obs = as.vector(src$propMat), spawn_time = 0)
  # The youngest age has no older-cohort survival residual. Avoid a design
  # column for its separate recruitment SD, which rec_settings already estimates.
  survival_keys <- cf$keyVarLogN; survival_keys[1L] <- survival_keys[2L]
  pop <- transform(grid, M = as.vector(src$natMor), N_sd_key = key(survival_keys[age + 1L]))
  obs <- list(catch = catch, index = index, weight = weight, maturity = maturity, pop_covariates = pop)
  check_tam_data(obs)
  settings <- list(N_settings = list(process = "iid", init = "latent", sd_form = ~0 + N_sd_key),
    rec_settings = list(process = "rw", init = "latent"),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(process = "off", mu_supplied = ~M),
    catch_settings = list(sd_form = ~0 + sd_key),
    index_settings = list(q_form = ~0 + q_key + power(N, by = q_power_group),
      sd_form = ~0 + sd_key, sd_supplied = ~relative_sd))
  dat <- do.call(prepare_tam, c(list(data = obs, years = years, ages = ages), settings))
  par <- make_par(dat)
  n_sd <- dat$process_sd$N$parameter
  key_index <- function(x) as.integer(sub(".*K", "", x)) + 1L
  par$log_r0 <- native$pl$logN[1L, 1L]
  par$log_n0[] <- native$pl$logN[-1L, 1L]
  par$log_r[] <- native$pl$logN[1L, -1L]
  par$log_n[] <- t(native$pl$logN[-1L, -1L])
  par$log_f[] <- t(native$pl$logF)
  par$log_sd_r <- native$pl$logSdLogN[1L]
  par[[n_sd]][] <- native$pl$logSdLogN[key_index(names(par[[n_sd]]))]
  par$log_sd_f <- native$pl$logSdLogFsta[1L]
  # SAM toInterval(theta, -1, 1, 2) equals tanh(theta).
  par$atanh_rho_f <- native$pl$itrans_rho[1L]
  par$log_q[] <- native$pl$logFpar[key_index(names(par$log_q))]
  par$log_q_power[] <- native$pl$logQpow[key_index(names(par$log_q_power))]
  for (nm in c("log_sd_catch", "log_sd_index")) par[[nm]][] <- native$pl$logSdLogObs[key_index(names(par[[nm]]))]
  tiny_joint <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  sam_parameters <- stockassessment::defpar(src, cf)
  for (nm in names(sam_parameters)) {
    stopifnot(length(sam_parameters[[nm]]) == length(native$pl[[nm]]))
    sam_parameters[[nm]][] <- native$pl[[nm]]
  }
  sam_parameters$missing <- numeric(sum(is.na(src$logobs)))
  sam_data <- c(src, list(forecast = list(), referencepoints = list()), cf,
    list(simFlag = c(1L, 1L), reportingLevel = 0L, resFlag = 0L))
  sam_joint <- TMB::MakeADFun(sam_data, sam_parameters, DLL = "stockassessment", silent = TRUE,
    map = list(logFScaleMSY = factor(NA), implicitFunctionDelta = factor(NA), splinePenalty = factor(NA)))
  sam <- TMB::MakeADFun(sam_data, sam_parameters, DLL = "stockassessment", silent = TRUE,
    random = c("logN", "logF"),
    map = list(logFScaleMSY = factor(NA), implicitFunctionDelta = factor(NA), splinePenalty = factor(NA)))
  pack <- function(pl, template) {
    for (nm in unique(names(template))) template[names(template) == nm] <- as.vector(pl[[nm]])
    template
  }
  to_sam <- function(p) {
    out <- sam_parameters
    out$logN[1L, 1L] <- p$log_r0; out$logN[-1L, 1L] <- p$log_n0
    out$logN[1L, -1L] <- p$log_r; out$logN[-1L, -1L] <- t(p$log_n)
    out$logF[] <- t(p$log_f)
    out$logSdLogN[1L] <- p$log_sd_r
    out$logSdLogN[key_index(names(p[[n_sd]]))] <- p[[n_sd]]
    out$logSdLogFsta[] <- p$log_sd_f
    out$itrans_rho[] <- p$atanh_rho_f
    out$logFpar[key_index(names(p$log_q))] <- p$log_q
    out$logQpow[key_index(names(p$log_q_power))] <- p$log_q_power
    for (nm in c("log_sd_catch", "log_sd_index")) out$logSdLogObs[key_index(names(p[[nm]]))] <- p[[nm]]
    out
  }
  markers <- par
  k <- 0L
  for (nm in names(markers)) {
    markers[[nm]][] <- k + seq_along(markers[[nm]])
    k <- k + length(markers[[nm]])
  }
  permutation <- as.integer(pack(to_sam(markers), sam_joint$par))
  stopifnot(identical(sort(permutation), seq_along(tiny_joint$par)))
  ordered <- order(permutation)
  checks <- list()
  for (attempt in 1:2) {
    x <- tiny_joint$par
    if (attempt == 2L) { set.seed(761); x <- x + rnorm(length(x), sd = .001) }
    p <- tiny_joint$env$parList(x)
    y <- pack(to_sam(p), sam_joint$par)
    gradient <- sam_joint$gr(y)[ordered]
    hessian <- sam_joint$he(y)[ordered, ordered]
    checks[[attempt]] <- data.frame(point = c("native", "perturbed")[attempt],
      joint = tiny_joint$fn(x) - sam_joint$fn(y),
      gradient = max(abs(as.numeric(tiny_joint$gr(x)) - as.numeric(gradient))),
      hessian = max(abs(tiny_joint$he(x) - hessian)))
  }
  checks <- do.call(rbind, checks)
  sam_report <- sam_joint$report(pack(to_sam(par), sam_joint$par))
  tiny_report <- tiny_joint$report(tiny_joint$par)
  fleet <- ifelse(dat$obs_map$type == "catch", 1L,
    match(dat$obs_map$survey, c("IBTS Q1", "IBTS Q3+Q4")) + 1L)
  observation_order <- match(paste(aux$year, aux$fleet, aux$age),
    paste(dat$obs_map$year, fleet, dat$obs_map$age))
  stopifnot(!anyNA(observation_order))
  prediction_difference <- max(abs(tiny_report$log_pred[observation_order] - sam_report$predObs))
  cat("Maximum log-prediction difference:", prediction_difference, "\n")
  print(checks)
  write.csv(checks, file.path(output, "joint_checks.csv"), row.names = FALSE)
  ran <- c("log_f", "log_r", "log_r0", "log_n0", "log_n")
  tiny <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, random = ran, silent = TRUE)
  marginal <- data.frame(tinyAM = tiny$fn(tiny$par), SAM = sam$fn(sam$par),
    n_fixed = length(tiny$par), n_random = sum(tiny$env$lrandom()), n_observed = length(dat$observed),
    max_log_prediction_difference = prediction_difference)
  print(marginal)
  write.csv(marginal, file.path(output, "marginal_checks.csv"), row.names = FALSE)
  stopifnot(marginal$n_fixed == 37L, marginal$n_random == 990L,
    marginal$n_observed == 1153L, max(abs(checks$joint)) < 1e-6,
    max(checks$gradient) < 1e-6, max(checks$hessian) < 1e-5,
    prediction_difference < 1e-8, abs(marginal$tinyAM - marginal$SAM) < 1e-5)
  fitted <- do.call(fit_tam, c(list(data = obs, years = years, ages = ages,
    start_par = par, silent = TRUE), settings))
  fitted_par <- fitted$obj$env$parList(par = fitted$parameter_values)
  fitted_outer <- pack(to_sam(fitted_par), sam$par)
  fit_checks <- data.frame(
    objective_difference = fitted$opt$objective - native$opt$objective,
    max_fixed_difference = max(abs(fitted_outer - sam$par)),
    max_relative_N = max(abs(fitted$rep$N / exp(t(native$pl$logN)) - 1)),
    max_relative_F = max(abs(fitted$rep$F / exp(t(native$pl$logF)) - 1)),
    max_gradient = max(abs(fitted$gradient)), converged = fitted$is_converged)
  write.csv(fit_checks, file.path(output, "fit_checks.csv"), row.names = FALSE)
  print(fit_checks)
  saveRDS(list(fit = fitted, native = native, joint_checks = checks,
    marginal_checks = marginal, input_md5 = tools::md5sum(args[1L]),
    fit_checks = fit_checks, SAM_revision = revision), file.path(output, "comparison.rds"))
  print(check_tam(fitted))
  cat("Fitted objective difference:", fitted$opt$objective - native$opt$objective, "\n")
  vis_tam(fitted, output_file = normalizePath(file.path(output, "haddock.html"), mustWork = FALSE),
    open_file = FALSE, render_args = list(quiet = TRUE),
    background = "### Exact-likelihood SAM comparison\nSource-aligned settings; see the comparison report for numerical checks and limitations.")
}, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
}), error = function(e) {
  writeLines(unique(warnings), file.path(output, "warnings.txt"))
  writeLines(conditionMessage(e), file.path(output, "failure.txt"))
  stop(e)
})
writeLines(unique(warnings), file.path(output, "warnings.txt"))
