.assessment_root <- file.path("analysis", "comp_assessments")
for (.assessment_helper in c(
  "read_database.R", "database_to_tam_obs.R", "database_to_tam_ref.R",
  "audit_assumptions.R"
)) {
  sys.source(file.path(.assessment_root, "R", .assessment_helper),
             envir = environment())
}
rm(.assessment_helper)

.assessment_diagnostics <- function(assessment_id, database, status,
                                    fit = NULL, elapsed = NA_real_, reason = "") {
  gradient <- if (is.null(fit)) numeric() else
    tryCatch(fit$obj$gr(fit$opt$par), error = function(e) numeric())
  data.frame(
    assessment_id = assessment_id,
    database_revision = if (length(database$commit)) database$commit[[1]] else NA_character_,
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

.assessment_percent_differences <- function(fit, reference,
                                             scales = c(ssb = 1e-3,
                                               recruitment = 1e-3, N = 1e-3,
                                               F = 1, M = 1)) {
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

.assessment_comparison_summary <- function(differences, assessment_id) {
  if (!nrow(differences)) return(data.frame())
  summary <- do.call(rbind, lapply(split(differences, differences$metric), function(x) {
    valid <- x[is.finite(x$percent_difference), , drop = FALSE]
    terminal_year <- if (nrow(valid)) max(valid$year) else NA_integer_
    terminal <- if (nrow(valid)) valid[valid$year == terminal_year, , drop = FALSE] else valid
    data.frame(
      assessment_id = assessment_id,
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
  rownames(summary) <- NULL
  summary
}

run_assessment <- function(assessment_id, database = NULL, fit = TRUE,
                            dashboard = FALSE, cache = FALSE) {
  if (is.null(database)) database <- read_database()
  if (length(fit) != 1L || is.na(fit) || !is.logical(fit) ||
      length(dashboard) != 1L || is.na(dashboard) || !is.logical(dashboard) ||
      length(cache) != 1L || is.na(cache) || !is.logical(cache)) {
    stop("fit, dashboard, and cache must each be one TRUE/FALSE value.", call. = FALSE)
  }

  source <- read_assessment(assessment_id, database)
  audit <- audit_assumptions(assessment_id, source$assumptions)
  result <- list(
    assessment_id = assessment_id,
    source = source,
    obs = NULL,
    settings = NULL,
    fit = NULL,
    ref = NULL,
    audit = audit,
    diagnostics = NULL,
    differences = NULL,
    summary = NULL,
    background = NULL,
    dashboard_file = NULL
  )

  script <- file.path(.assessment_root, "scripts", "translation", "stocks",
                      paste0(assessment_id, ".R"))
  if (!file.exists(script)) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "blocked_no_stock_spec",
      reason = "No stock translation recipe is available."
    )
    return(result)
  }

  stock_env <- new.env(parent = environment())
  sys.source(script, envir = stock_env)
  translated <- tryCatch(stock_env$translate_stock(source), error = identity)
  if (inherits(translated, "error")) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "translation_failed",
      reason = conditionMessage(translated)
    )
    return(result)
  }
  result$obs <- translated$obs
  result$settings <- translated$settings
  result$background <- translated$background

  if (!fit) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "not_fitted", reason = "fit = FALSE."
    )
    return(result)
  }

  started <- Sys.time()
  fit_args <- c(list(obs = translated$obs, years = translated$years,
                     ages = translated$ages, silent = TRUE), translated$settings)
  if (!is.null(translated$start_par)) fit_args$start_par <- translated$start_par
  fitted <- tryCatch(do.call(tinyAM::fit_tam, fit_args), error = identity)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fitted, "error")) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "fit_failed", elapsed = elapsed,
      reason = conditionMessage(fitted)
    )
    return(result)
  }
  result$fit <- fitted
  result$diagnostics <- .assessment_diagnostics(
    assessment_id, database,
    if (isTRUE(fitted$is_converged)) "converged" else "not_converged",
    fit = fitted, elapsed = elapsed
  )
  if (!isTRUE(fitted$is_converged)) return(result)

  outputs <- if (is.null(translated$comparison_outputs)) {
    database$outputs
  } else {
    translated$comparison_outputs
  }
  result$ref <- tryCatch(database_to_tam_ref(
    assessment_id, outputs, obs = translated$obs,
    years = translated$years, ages = translated$ages,
    terminal_year = source$assessment$terminal_year[[1]],
    age_plus_group = translated$age_plus_group,
    comparison_scales = translated$comparison_scales,
    template = fitted
  ), error = identity)
  if (inherits(result$ref, "error")) {
    result$diagnostics$status <- "reference_failed"
    result$diagnostics$reason <- conditionMessage(result$ref)
    result$ref <- NULL
    return(result)
  }

  result$differences <- .assessment_percent_differences(fitted, result$ref)
  result$summary <- .assessment_comparison_summary(
    result$differences, assessment_id
  )

  cache_dir <- file.path(.assessment_root, "results", "cache", assessment_id)
  if (cache) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(fitted, file.path(cache_dir, "fit.rds"))
  }
  if (dashboard) {
    result$dashboard_file <- if (cache) {
      file.path(cache_dir, "dashboard.html")
    } else {
      tempfile(pattern = paste0(assessment_id, "_"), fileext = ".html")
    }
    tinyAM::vis_tam(
      model_list = list(Accepted = result$ref, tinyAM = fitted),
      background = result$background,
      output_file = result$dashboard_file,
      open_file = interactive()
    )
  }
  result
}
