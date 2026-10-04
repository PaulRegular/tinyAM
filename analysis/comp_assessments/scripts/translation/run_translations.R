root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "database_to_tam_ref.R"))
source(file.path(root, "R", "audit_assumptions.R"))

database <- read_committed_database()
assessments <- database$assessments[
  database$assessments$is_current & database$assessments$is_applied, , drop = FALSE
]
requested_assessments <- commandArgs(trailingOnly = TRUE)
if (length(requested_assessments)) {
  unknown <- setdiff(requested_assessments, assessments$assessment_id)
  if (length(unknown)) stop("Unknown assessment_id: ", paste(unknown, collapse = ", "))
  assessments <- assessments[match(requested_assessments, assessments$assessment_id), , drop = FALSE]
}
readiness <- read.csv(file.path(root, "results", "fit_readiness.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
stock_scripts <- file.path(root, "scripts", "translation", "stocks")
output_root <- file.path(root, "results", "translations")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

diagnostic_row <- function(assessment_id, status, reason = "", fit = NULL,
                           elapsed = NA_real_) {
  gradient <- if (is.null(fit)) numeric() else
    tryCatch(fit$obj$gr(fit$opt$par), error = function(e) numeric())
  data.frame(
    assessment_id = assessment_id,
    status = status,
    is_converged = if (is.null(fit)) NA else isTRUE(fit$is_converged),
    optimizer_code = if (is.null(fit)) NA_integer_ else fit$opt$convergence,
    optimizer_message = if (is.null(fit) || is.null(fit$opt$message)) "" else fit$opt$message,
    objective = if (is.null(fit)) NA_real_ else fit$opt$objective,
    max_abs_gradient = if (length(gradient)) max(abs(gradient), na.rm = TRUE) else NA_real_,
    sdreport_success = if (is.null(fit)) NA else !is.null(fit$sdrep),
    positive_definite_hessian = if (is.null(fit)) NA else isTRUE(fit$sdrep$pdHess),
    n_fixed_effects = if (is.null(fit)) NA_integer_ else length(fit$opt$par),
    n_random_effects = if (is.null(fit)) NA_integer_ else length(fit$obj$env$random),
    elapsed_seconds = elapsed,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

fit_diagnostics <- function(fit, assessment_id, elapsed) {
  diagnostic_row(assessment_id,
                 if (isTRUE(fit$is_converged)) "converged" else "not_converged",
                 fit = fit, elapsed = elapsed)
}

percent_differences <- function(fit, reference, scales = c(ssb = 1e-3,
    recruitment = 1e-3, N = 1e-3, F = 1, M = 1)) {
  if (!is.null(reference$comparison_scales)) {
    scales[names(reference$comparison_scales)] <- reference$comparison_scales
  }
  rows <- list()
  for (metric in names(scales)) {
    source <- reference$pop[[metric]]
    tiny <- fit$pop[[metric]]
    if (is.null(source) || is.null(tiny)) next
    keys <- if (metric %in% c("N", "F", "M")) c("year", "age") else "year"
    source <- source[c(keys, "est")]
    tiny <- tiny[c(keys, "est")]
    names(source)[ncol(source)] <- "source"
    names(tiny)[ncol(tiny)] <- "tinyAM"
    common <- merge(source, tiny, by = keys)
    if (!nrow(common)) next
    if (!"age" %in% names(common)) common$age <- NA_integer_
    common$tinyAM <- common$tinyAM * scales[[metric]]
    common$metric <- metric
    common$percent_difference <- ifelse(
      is.finite(common$source) & common$source != 0,
      100 * (common$tinyAM - common$source) / common$source,
      NA_real_
    )
    rows[[metric]] <- common[c("metric", "year", "age", "source", "tinyAM",
                               "percent_difference")]
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

diagnostics <- lapply(assessments$assessment_id, function(assessment_id) {
  script <- file.path(stock_scripts, paste0(assessment_id, ".R"))
  out_dir <- file.path(output_root, assessment_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(script)) {
    blocker <- readiness$missing_items[match(assessment_id, readiness$assessment_id)]
    if (is.na(blocker) || !nzchar(blocker) || blocker == "none") {
      blocker <- "The source-to-tinyAM model mapping has not been specified."
    }
    return(diagnostic_row(assessment_id, "blocked_no_stock_spec",
                          blocker))
  }

  source_data <- read_committed_assessment(assessment_id, database)
  stock_env <- new.env(parent = environment())
  sys.source(script, envir = stock_env)
  translated <- tryCatch(stock_env$translate_stock(source_data), error = identity)
  if (inherits(translated, "error")) {
    return(diagnostic_row(assessment_id, "translation_failed",
                          conditionMessage(translated)))
  }

  writeLines(database$commit, file.path(out_dir, "database_revision.txt"))
  writeLines(translated$background, file.path(out_dir, "background.md"))
  saveRDS(translated$obs, file.path(out_dir, "obs.rds"))
  write.csv(attr(translated$obs, "translation")$source_provenance,
            file.path(out_dir, "source_provenance.csv"), row.names = FALSE, na = "")
  write.csv(audit_assumptions(assessment_id, source_data$assumptions),
            file.path(out_dir, "assumption_audit.csv"), row.names = FALSE, na = "")
  settings_output <- capture.output(dput(translated$settings))
  writeLines(sub("[ \\t]+$", "", settings_output),
             file.path(out_dir, "fit_settings.R"))

  started <- Sys.time()
  fit_args <- c(list(obs = translated$obs, years = translated$years,
                     ages = translated$ages, silent = TRUE), translated$settings)
  if (!is.null(translated$start_par)) fit_args$start_par <- translated$start_par
  fit <- tryCatch(do.call(tinyAM::fit_tam, fit_args),
                  error = identity)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fit, "error")) {
    return(diagnostic_row(assessment_id, "fit_failed", conditionMessage(fit),
                          elapsed = elapsed))
  }
  diagnostics <- fit_diagnostics(fit, assessment_id, elapsed)
  if (!isTRUE(fit$is_converged)) {
    write.csv(diagnostics, file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE, na = "")
    return(diagnostics)
  }

  comparison_outputs <- if (is.null(translated$comparison_outputs)) {
    database$outputs
  } else {
    translated$comparison_outputs
  }
  reference <- database_to_tam_ref(
    assessment_id, comparison_outputs, obs = translated$obs,
    years = translated$years, ages = translated$ages,
    terminal_year = source_data$assessment$terminal_year[[1]],
    age_plus_group = translated$age_plus_group,
    comparison_scales = translated$comparison_scales,
    template = fit
  )
  saveRDS(fit, file.path(out_dir, "tinyAM_fit.rds"))
  saveRDS(reference, file.path(out_dir, "assessment_outputs.rds"))
  write.csv(diagnostics, file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE, na = "")

  differences <- percent_differences(fit, reference)
  write.csv(differences, file.path(out_dir, "percent_differences.csv"),
            row.names = FALSE, na = "")
  summary <- do.call(rbind, lapply(split(differences, differences$metric), function(x) {
    valid <- x[is.finite(x$percent_difference), , drop = FALSE]
    terminal_year <- if (nrow(valid)) max(valid$year) else NA_integer_
    terminal <- if (nrow(valid)) {
      valid[valid$year == terminal_year, , drop = FALSE]
    } else valid
    data.frame(
      metric = x$metric[[1]],
      n = nrow(valid),
      mean_absolute_percent_difference = if (nrow(valid)) mean(abs(valid$percent_difference)) else NA_real_,
      median_absolute_percent_difference = if (nrow(valid)) stats::median(abs(valid$percent_difference)) else NA_real_,
      terminal_year = terminal_year,
      terminal_mean_percent_difference = if (nrow(terminal)) mean(terminal$percent_difference) else NA_real_,
      trend_correlation = if (nrow(valid) > 1L) {
        tryCatch(stats::cor(valid$source, valid$tinyAM), error = function(e) NA_real_)
      } else NA_real_
    )
  }))
  write.csv(summary, file.path(out_dir, "comparison_summary.csv"), row.names = FALSE, na = "")

  dashboard <- tryCatch({
    tinyAM::vis_tam(
      model_list = list(Accepted = reference, tinyAM = fit),
      background = translated$background,
      output_file = file.path(out_dir, "dashboard.html"),
      open_file = FALSE
    )
    NULL
  }, error = identity)
  if (inherits(dashboard, "error")) {
    diagnostics$status <- "converged_dashboard_failed"
    diagnostics$reason <- conditionMessage(dashboard)
  }
  write.csv(diagnostics, file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE, na = "")
  diagnostics
})

diagnostics <- do.call(rbind, diagnostics)
diagnostics_path <- file.path(output_root, "fit_diagnostics.csv")
if (length(requested_assessments) && file.exists(diagnostics_path)) {
  previous <- read.csv(diagnostics_path, stringsAsFactors = FALSE)
  if (identical(names(previous), names(diagnostics))) {
    previous <- previous[!previous$assessment_id %in% diagnostics$assessment_id, , drop = FALSE]
    diagnostics <- rbind(previous, diagnostics)
  }
}
write.csv(diagnostics, diagnostics_path, row.names = FALSE, na = "")
print(diagnostics, row.names = FALSE)
