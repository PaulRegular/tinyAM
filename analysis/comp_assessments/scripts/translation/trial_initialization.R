# Source-verified SAM boundaries; generated fits and dashboards stay local.
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
root <- "analysis/comp_assessments"
cases <- c(
  ices_bluewhiting_northeast_atlantic_2026 = "ices_bluewhiting_northeast_atlantic_2026/model.RData",
  ices_cod_northeast_arctic_2026 = "ices_cod_northeast_arctic_2026/baserun_model.RData",
  ices_haddock_north_sea_2026 = "ices_haddock_north_sea_2026/run_model.RData",
  ices_herring_western_baltic_2026 = "ices_herring_western_baltic_2026/model.RData",
  ices_saithe_north_sea_2026 = "ices_north_sea_saithe_2026/model.RData",
  ices_whiting_north_sea_2026 = "ices_whiting_north_sea_2026/NSwhiting_2026n_updateIBTSQ1_model.RData"
)
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested)) cases <- cases[requested]
stopifnot(!anyNA(cases))
out <- file.path(root, "results", "cache", "initialization_trials")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
database <- read_database()
revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
baseline_init <- setNames(rep("free", length(cases)), names(cases))
baseline_init["ices_whiting_north_sea_2026"] <- "exp"
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(tryCatch(expr, error = identity), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
records <- list()
for (id in names(cases)) {
  native <- new.env()
  native_file <- file.path(root, "source_cache", cases[[id]])
  load(native_file, native)
  stopifnot(native$fit$conf$initState == 0, length(native$fit$pl$initN) == 0)
  baseline <- capture(run_assessment(id, database, silent = TRUE))
  result <- baseline$value
  if (inherits(result, "error") || !inherits(result$fit, "tam_fit")) {
    records[[id]] <- data.frame(assessment_id = id, trial = "baseline", initialization = "baseline",
      converged = FALSE, pdHess = FALSE, gradient = NA_real_, raw_gradient = NA_real_, optimizer = NA_integer_,
      message = if (inherits(result, "error")) conditionMessage(result) else result$diagnostics$reason,
      warnings = paste(baseline$warnings, collapse = " | "), revision = revision,
      native_file = cases[[id]])
    next
  }
  if (result$fit$dat$N_settings$init != baseline_init[[id]]) {
    original <- capture(update(result$fit, N_settings = list(init = baseline_init[[id]]), silent = TRUE))
    if (inherits(original$value, "error")) stop(original$value)
    result$fit <- original$value
    baseline$warnings <- c(baseline$warnings, original$warnings)
  }
  start <- tinyAM:::.tam_parameter_summary(result$fit, "Estimate")
  start$log_n0 <- setNames(log(result$fit$rep$N[1, -1]), as.character(result$fit$dat$ages[-1]))
  candidate <- capture(update(result$fit, N_settings = list(init = "latent"),
    start_par = start, silent = TRUE))
  fits <- list(baseline = result$fit, latent = candidate$value)
  if (inherits(candidate$value, "tam_fit") && isTRUE(candidate$value$is_converged)) {
    literal <- capture(update(result$fit, N_settings = list(init = "latent"),
      start_par = result$fit$refit_args$start_par, silent = TRUE))
    fits$script_start <- literal$value
  }
  for (name in names(fits)) {
    fit <- fits[[name]]
    good <- inherits(fit, "tam_fit")
    records[[paste(id, name)]] <- data.frame(
      assessment_id = id, trial = name, initialization = if (good) fit$dat$N_settings$init else "latent",
      converged = good && isTRUE(fit$is_converged),
      pdHess = good && isTRUE(fit$sdrep$pdHess),
      gradient = if (good) fit$diagnostics$max_gradient else NA_real_,
      raw_gradient = if (good) fit$diagnostics$raw_max_gradient else NA_real_,
      optimizer = if (good) fit$opt$convergence else NA_integer_,
      message = if (good) fit$opt$message else conditionMessage(fit),
      warnings = paste(switch(name, baseline = baseline$warnings,
        latent = candidate$warnings, script_start = literal$warnings), collapse = " | "),
      revision = revision, native_file = cases[[id]])
    if (!good) next
    saveRDS(fit, file.path(out, paste0(id, "_", name, ".rds")))
    common <- capture(.assessment_comparison(fit, result))
    if (!inherits(common$value, "error")) {
      write.csv(common$value$summary, file.path(out, paste0(id, "_", name, "_comparison.csv")), row.names = FALSE)
    }
  }
  if (inherits(candidate$value, "tam_fit")) {
    dashboard <- capture(view_assessment(candidate$value, result,
      output_file = file.path(out, paste0(id, "_latent.html")), open_file = FALSE))
    if (inherits(dashboard$value, "error")) {
      writeLines(conditionMessage(dashboard$value), file.path(out, paste0(id, "_dashboard_error.txt")))
    }
  }
  write.csv(do.call(rbind, records), file.path(out, "diagnostics.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, records), file.path(out, "diagnostics.csv"), row.names = FALSE)
