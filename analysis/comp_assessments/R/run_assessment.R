.assessment_root <- file.path("analysis", "comp_assessments")
for (.assessment_helper in c(
  "translation_registry.R", "read_database.R", "database_to_tam_obs.R", "database_to_tam_ref.R",
  "audit_assumptions.R", "compare_assessment.R"
)) {
  sys.source(file.path(.assessment_root, "R", .assessment_helper),
    envir = environment()
  )
}
rm(.assessment_helper)


.assessment_catch_reporting <- function(fit, reporting) {
  if (is.null(reporting)) {
    return(fit)
  }
  catch <- fit$obs_pred$catch
  weights <- reporting$weights
  key <- function(x) paste(x$year, x$age)
  index <- match(key(catch), key(weights))
  if (anyDuplicated(key(weights)) || anyNA(index) ||
    any(!is.finite(weights$weight) | weights$weight <= 0)) {
    cli::cli_abort("Catch reporting requires one positive catch weight for every predicted year and age.")
  }
  count_scale <- tinyAM:::.unit_scale(unique(catch$unit), "numbers")
  if (length(count_scale) != 1L || !is.finite(count_scale)) {
    cli::cli_abort("Catch reporting requires known, consistent abundance units in the fitted catch table.")
  }
  predicted <- tapply(catch$pred * count_scale * weights$weight[index] / 1000, catch$year, sum)
  totals <- reporting$totals
  if (anyDuplicated(totals$year) || any(!is.finite(totals$yield) | totals$yield < 0)) {
    cli::cli_abort("Catch reporting requires unique years and non-negative total yields in kg.")
  }
  # Reporting only: source total biomass and all-age predictions use catch weights.
  fit$pop$total_yield$est <- totals$yield[match(fit$pop$total_yield$year, totals$year)] / 1000
  fit$pop$total_yield_pred$est <- as.numeric(predicted[as.character(fit$pop$total_yield_pred$year)])
  fit$pop$total_yield$unit <- fit$pop$total_yield_pred$unit <- "tonnes"
  fit$rep$total_yield <- fit$pop$total_yield$est
  fit$rep$total_yield_pred <- fit$pop$total_yield_pred$est
  fit
}

.assessment_diagnostics <- function(assessment_id, database, status,
                                    fit = NULL, elapsed = NA_real_, reason = "") {
  opt <- if (is.list(fit$opt)) fit$opt else NULL
  sdrep_ok <- is.list(fit$sdrep) && !inherits(fit$sdrep, "try-error")
  if (!is.null(fit) && is.null(opt)) {
    status <- "fit_failed"
    reason <- paste("Optimizer did not return a result:", paste(as.character(fit$opt), collapse = " "))
  }
  checks <- if (inherits(fit, "tam_fit")) tinyAM::check_tam(fit) else NULL
  attempts <- attr(fit, "initialization_attempts")
  if (!is.null(attempts)) {
    reason <- paste(c(reason, paste0(attempts$initialization, ": ",
      ifelse(attempts$is_converged, "converged", "failed numerical checks"),
      "; gradient=", signif(attempts$max_abs_gradient, 3),
      "; ", attempts$notes)), collapse = " | ")
  }
  data.frame(
    assessment_id = assessment_id,
    database_revision = if (length(database$commit)) database$commit[[1]] else NA_character_,
    database_source = if (length(database$database_source)) database$database_source[[1]] else NA_character_,
    database_ref = if (length(database$database_ref)) database$database_ref[[1]] else NA_character_,
    database_selection = if (length(database$database_selection)) database$database_selection[[1]] else NA_character_,
    database_key = if (length(database$database_key)) database$database_key[[1]] else
      if (length(database$commit)) database$commit[[1]] else NA_character_,
    tinyAM_revision = .database_head("."),
    translation_tier = if (length(database$translation_tier)) database$translation_tier[[1]] else NA_character_,
    status = status,
    is_converged = if (is.null(checks)) NA else checks$is_converged,
    optimizer_code = if (is.null(opt)) NA_integer_ else opt$convergence,
    optimizer_message = if (is.null(opt) || is.null(opt$message)) "" else opt$message,
    objective = if (is.null(opt)) NA_real_ else opt$objective,
    max_abs_gradient = if (is.null(checks)) NA_real_ else checks$max_gradient,
    raw_max_abs_gradient = if (is.null(checks)) NA_real_ else checks$raw_max_gradient,
    sdreport_success = if (is.null(fit)) NA else sdrep_ok,
    positive_definite_hessian = if (is.null(fit)) NA else sdrep_ok && isTRUE(fit$sdrep$pdHess),
    n_fixed_effects = if (is.null(opt)) NA_integer_ else length(opt$par),
    n_random_effects = if (is.null(fit)) NA_integer_ else length(fit$obj$env$random),
    elapsed_seconds = elapsed,
    reason = reason,
    stringsAsFactors = FALSE
  )
}


.assessment_comparison_summary <- function(differences, assessment_id) {
  if (!nrow(differences)) {
    return(data.frame())
  }
  summary <- do.call(rbind, lapply(split(differences, differences$metric), function(x) {
    valid <- x[is.finite(x$source) & is.finite(x$tinyAM), , drop = FALSE]
    percent <- valid$percent_difference[is.finite(valid$percent_difference)]
    terminal_year <- if (nrow(valid)) max(valid$year) else NA_integer_
    terminal <- if (nrow(valid)) valid[valid$year == terminal_year, , drop = FALSE] else valid
    trends <- if (nrow(valid)) split(valid, ifelse(is.na(valid$age), "all", valid$age)) else list()
    has_variation <- function(value) {
      diff(range(value)) > sqrt(.Machine$double.eps) * max(1, max(abs(value)))
    }
    correlations <- vapply(trends, function(z) {
      if (nrow(z) > 1L && has_variation(z$source) && has_variation(z$tinyAM)) {
        stats::cor(z$source, z$tinyAM)
      } else {
        NA_real_
      }
    }, numeric(1))
    trend_correlation <- if (any(is.finite(correlations))) {
      mean(correlations[is.finite(correlations)])
    } else {
      NA_real_
    }
    data.frame(
      assessment_id = assessment_id,
      metric = x$metric[[1]],
      unit = if ("unit" %in% names(x)) x$unit[[1]] else "",
      comparison_status = if ("comparison_status" %in% names(x)) x$comparison_status[[1]] else "matched",
      definition = if ("definition" %in% names(x)) x$definition[[1]] else "",
      reason = if ("reason" %in% names(x)) x$reason[[1]] else "",
      n = nrow(valid),
      mean_absolute_difference = if (nrow(valid)) mean(abs(valid$tinyAM - valid$source)) else NA_real_,
      mean_absolute_percent_difference = if (length(percent)) mean(abs(percent)) else NA_real_,
      median_absolute_percent_difference = if (length(percent)) stats::median(abs(percent)) else NA_real_,
      terminal_year = terminal_year,
      terminal_mean_percent_difference = if (any(is.finite(terminal$percent_difference))) {
        mean(terminal$percent_difference[is.finite(terminal$percent_difference)])
      } else {
        NA_real_
      },
      trend_correlation = trend_correlation
    )
  }))
  rownames(summary) <- NULL
  summary
}


print_sources <- function(assessment) {
  fields <- c(
    assessment_url = "Assessment report",
    framework_url = "Framework / methodology",
    data_url = "Assessment data",
    model_url = "Native model / fitted run",
    repository_url = "Assessment repository"
  )
  if (is.null(assessment) || is.null(names(assessment))) {
    return(character())
  }

  urls <- vapply(names(fields), function(field) {
    if (!field %in% names(assessment) ||
      length(assessment[[field]]) != 1L) {
      return(NA_character_)
    }
    value <- assessment[[field]][[1L]]
    if (length(value) != 1L || is.na(value)) {
      return(NA_character_)
    }
    as.character(value)
  }, character(1))

  keep <- !is.na(urls) & nzchar(trimws(urls))
  if (!any(keep)) {
    return(character())
  }

  available_urls <- urls[keep]
  available_labels <- unname(fields[keep])
  unique_urls <- unique(available_urls)
  links <- vapply(unique_urls, function(url) {
    labels <- available_labels[available_urls == url]
    paste0("- [", paste(labels, collapse = " / "), "](", url, ")")
  }, character(1), USE.NAMES = FALSE)

  c("#### Assessment documentation", links)
}

.assessment_stock <- function(source, do_fit, silent, tier = NULL) {
  entry <- .assessment_translation(source$assessment$assessment_id[[1]], tier)
  stock_env <- new.env(parent = environment(.assessment_stock))
  stock_env$source <- source
  stock_env$translation_tier <- entry$tier[[1]]
  stock_env$do_fit <- do_fit
  stock_env$silent <- silent
  stock_env$fit_stage <- "translation"
  stock_env$fit_started <- Sys.time()
  script <- entry$script[[1]]
  error <- tryCatch(
    {
      sys.source(script, envir = stock_env)
      NULL
    },
    error = identity
  )
  list(env = stock_env, error = error)
}

.assessment_context <- function(source, translated) {
  source$translation_tier <- translated$translation_tier
  keys <- c(
    "age_plus_group", "comparison_scales", "comparison_outputs",
    "comparison_aggregates", "comparison_age_groups", "comparison_definitions",
    "catch_reporting"
  )
  list(
    assessment_id = source$assessment$assessment_id[[1]], source = source,
    translation_tier = translated$translation_tier,
    background = translated$background,
    comparison = stats::setNames(lapply(keys, function(key) translated[[key]]), keys)
  )
}

.assessment_comparison <- function(fit, assessment) {
  if (!inherits(fit, "tam_fit")) cli::cli_abort("Supply a fitted tinyAM model as {.arg fit}.")
  if (!is.list(assessment) || is.null(assessment$source$assessment) ||
    !is.list(assessment$comparison)) {
    cli::cli_abort("Load the assessment context with {.file run_stock.R} or {.fn run_assessment}.")
  }
  if (!all(c("dat", "pop", "rep", "obs_pred", "fixed_par", "random_par") %in% names(fit)) ||
    !any(vapply(fit$obs_pred, function(x) {
      is.data.frame(x) &&
        "pred" %in% names(x) && any(is.finite(x$pred))
    }, logical(1)))) {
    cli::cli_abort("This fit has no usable reporting predictions. Inspect {.code check_tam(fit)} and refit before previewing.")
  }
  comparison <- assessment$comparison
  plus <- comparison$age_plus_group
  source_ages <- suppressWarnings(as.numeric(c(
    assessment$source$inputs$age,
    assessment$source$outputs$age
  )))
  if (!is.null(plus)) {
    if (!is.numeric(plus) || length(plus) != 1L || !is.finite(plus) ||
      plus != trunc(plus) || !plus %in% c(source_ages, fit$dat$ages)) {
      cli::cli_abort("The comparison plus age must be a recorded source age or the fitted plus age.")
    }
    # The supplied fit determines the terminal group, including for age subsets.
    plus <- max(fit$dat$ages)
  }
  mappings <- comparison$comparison_age_groups
  if (!is.null(mappings) && (!is.list(mappings) || (length(mappings) &&
    (is.null(names(mappings)) || anyNA(names(mappings)) || any(!nzchar(names(mappings))) ||
      anyDuplicated(names(mappings)))))) {
    cli::cli_abort("{.code assessment$comparison$comparison_age_groups} must be a uniquely named list of metrics.")
  }
  for (metric in names(mappings)) {
    groups <- comparison$comparison_age_groups[[metric]]
    if (!is.list(groups) || (length(groups) &&
      (is.null(names(groups)) || anyNA(names(groups)) || any(!nzchar(names(groups))) ||
        anyDuplicated(names(groups))))) {
      cli::cli_abort("Comparison age groups for {metric} must be a uniquely named list.")
    }
    valid <- vapply(groups, function(ages) {
      is.numeric(ages) && length(ages) > 0L &&
        all(is.finite(ages)) && all(ages == trunc(ages)) &&
        all(ages %in% c(source_ages, fit$dat$ages)) &&
        !anyDuplicated(ages)
    }, logical(1))
    if (!all(valid) || anyDuplicated(unlist(groups, use.names = FALSE))) {
      cli::cli_abort("Comparison age groups for {metric} must map to distinct fitted ages. Update {.code assessment$comparison$comparison_age_groups} for this age range.")
    }
  }
  fit <- .assessment_catch_reporting(fit, comparison$catch_reporting)
  outputs <- if (is.null(comparison$comparison_outputs)) {
    assessment$source$outputs
  } else {
    comparison$comparison_outputs
  }
  obs <- if (is.null(fit$refit_args$data)) fit$dat$obs else fit$refit_args$data
  ref <- database_to_tam_ref(
    assessment$assessment_id, outputs,
    obs = obs,
    years = fit$dat$years, ages = fit$dat$ages,
    terminal_year = assessment$source$assessment$terminal_year[[1]],
    age_plus_group = plus, comparison_scales = comparison$comparison_scales,
    template = fit, assumptions = assessment$source$assumptions,
    comparison_aggregates = comparison$comparison_aggregates,
    comparison_age_groups = comparison$comparison_age_groups,
    comparison_definitions = comparison$comparison_definitions
  )
  differences <- .assessment_percent_differences(ref)
  historical <- is.na(differences$year) | differences$year %in% fit$dat$years[!fit$dat$is_proj]
  list(
    fit = fit, ref = ref, differences = differences,
    summary = .assessment_comparison_summary(differences[historical, ], assessment$assessment_id)
  )
}

# Preview a supplied fit without re-running its stock script or exporting results.
view_assessment <- function(fit, assessment, output_file = NULL, open_file = TRUE) {
  if (!is.logical(open_file) || length(open_file) != 1L || is.na(open_file)) {
    cli::cli_abort("{.arg open_file} must be TRUE or FALSE.")
  }
  if (!is.null(output_file) && (!is.character(output_file) || length(output_file) != 1L ||
    is.na(output_file) || !nzchar(output_file))) {
    cli::cli_abort("{.arg output_file} must be NULL or one HTML file path.")
  }
  source_key <- if (!is.null(assessment$source$database_key)) {
    .assessment_cache_key(assessment$source)
  } else NA_character_
  fit_key <- fit$assessment_provenance$database_key
  if (length(fit_key) && !is.na(source_key) && !identical(fit_key, source_key)) {
    cli::cli_abort("The supplied fit and assessment use different SAIO source revisions. Supply the original matching assessment context or refit from the selected data; a preview must not relabel historical estimates.")
  }
  result <- .assessment_comparison(fit, assessment)
  result$diagnostics <- .assessment_diagnostics(assessment$assessment_id,
    assessment$source, "preview",
    fit = result$fit
  )
  result$diagnostics$status <- if (isTRUE(result$diagnostics$is_converged)) {
    "converged"
  } else {
    "not_converged"
  }
  if (!isTRUE(result$diagnostics$is_converged)) {
    cli::cli_warn("This fit has not passed numerical convergence checks. The dashboard is exploratory; inspect {.code check_tam(fit)} before interpreting estimates or uncertainty.")
  }
  result$dashboard_file <- if (is.null(output_file)) {
    tempfile("assessment_preview_", fileext = ".html")
  } else {
    output_file
  }
  provenance <- paste0("### Source provenance\n\nAssessment source: ",
    if (is.na(source_key)) "unrecorded" else source_key,
    ". Fit source: ", if (length(fit_key)) fit_key else
      "unrecorded; supplied fit provenance has not been verified", ".")
  tinyAM::vis_tam(
    model_list = list(Official = result$ref, tinyAM = result$fit),
    background = c(if (!isTRUE(result$diagnostics$is_converged)) c(
      "### Troubleshooting fit", "Numerical checks did not pass. Inspect Diagnostics before interpreting estimates or uncertainty.", ""),
      assessment$background, "", provenance), output_file = result$dashboard_file,
    open_file = open_file
  )
  invisible(result)
}

# Rendering is deliberately separate from fitting and numerical eligibility.
render_assessment_dashboards <- function(runs, open_file = FALSE) {
  if (!is.list(runs) || any(!vapply(runs, is.list, logical(1)))) {
    cli::cli_abort("Supply the results returned by run_assessments().")
  }
  if (!is.logical(open_file) || length(open_file) != 1L || is.na(open_file)) {
    cli::cli_abort("open_file must be TRUE or FALSE.")
  }
  if (!length(runs)) return(data.frame(assessment_id = character(), status = character(),
    file = character(), error = character()))
  tiers <- unique(vapply(runs, function(x) x$translation_tier, character(1)))
  keys <- unique(vapply(runs, .assessment_result_key, character(1)))
  if (length(tiers) != 1L || length(keys) != 1L) {
    cli::cli_abort("Render one translation tier and SAIO snapshot at a time.")
  }
  tier <- tiers[[1]]
  database <- list(database_key = keys[[1]])
  manifest <- lapply(names(runs), function(id) {
    result <- runs[[id]]
    row <- data.frame(assessment_id = id, status = "skipped_no_fit", file = NA_character_,
      error = "No usable fitted model was returned.", stringsAsFactors = FALSE)
    if (!inherits(result$fit, "tam_fit")) {
      reason <- result$diagnostics$reason
      if (length(reason) && !is.na(reason[[1]]) && nzchar(reason[[1]])) row$error <- reason[[1]]
      return(row)
    }
    row$file <- file.path(.assessment_cache(id, tier, database), "dashboard.html")
    dir.create(dirname(row$file), recursive = TRUE, showWarnings = FALSE)
    troubleshooting <- !isTRUE(result$fit$is_converged)
    cli::cli_inform("Rendering {id}{if (troubleshooting) ' (troubleshooting: numerical checks did not pass)' else ''}.")
    error <- tryCatch({
      view_assessment(result$fit, result, output_file = row$file, open_file = FALSE)
      NULL
    }, error = conditionMessage)
    row$status <- if (!is.null(error)) "render_failed" else if (troubleshooting) "troubleshooting" else "rendered"
    row$error <- if (is.null(error)) NA_character_ else error
    if (identical(row$status, "troubleshooting")) cli::cli_inform("Troubleshooting dashboard: {.file {row$file}}")
    row
  })
  manifest <- do.call(rbind, manifest)
  folder <- .assessment_cache_root(tier, database)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(manifest, file.path(folder, "dashboard_manifest.csv"), row.names = FALSE, na = "")
  index <- .index_assessment_dashboards(tier, database)
  cli::cli_inform(c("Dashboards: {sum(manifest$status == 'rendered')} rendered; {sum(manifest$status == 'troubleshooting')} for troubleshooting; {sum(manifest$status == 'render_failed')} failed; {sum(manifest$status == 'skipped_no_fit')} skipped.",
    "i" = "Fits and numerical comparison eligibility are unchanged. Inspect dashboard_manifest.csv and fit_diagnostics.csv before interpreting results."))
  if (open_file) utils::browseURL(index)
  manifest
}

.assessment_log <- function(path, text) {
  if (!is.null(path)) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text, "\n",
    file = path, append = TRUE)
}

.assessment_save_rds <- function(object, file) {
  # Keep the previous complete checkpoint if serialization is interrupted.
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(file), "-"), tmpdir = dirname(file))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, file)) cli::cli_abort("Could not replace checkpoint {.file {file}}.")
  invisible(file)
}

run_assessment <- function(assessment_id, database = NULL, fit = TRUE,
                           cache = FALSE, silent = TRUE,
                           log_file = NULL, tier = NULL) {
  entry <- .assessment_translation(assessment_id, tier, required = FALSE)
  tier <- if (is.null(entry)) "curated" else entry$tier[[1]]
  if (is.null(database)) database <- read_database()
  database$translation_tier <- tier
  if (length(fit) != 1L || is.na(fit) || !is.logical(fit) ||
    length(cache) != 1L || is.na(cache) || !is.logical(cache)) {
    stop("fit and cache must each be one TRUE/FALSE value.", call. = FALSE)
  }

  source <- read_assessment(assessment_id, database)
  audit <- audit_assumptions(assessment_id, source$assumptions)
  result <- list(
    assessment_id = assessment_id,
    translation_tier = tier,
    source = source,
    obs = NULL,
    fit = NULL,
    ref = NULL,
    audit = audit,
    diagnostics = NULL,
    differences = NULL,
    summary = NULL,
    background = NULL,
    dashboard_file = NULL
  )

  script <- if (is.null(entry)) "" else entry$script[[1]]
  if (!file.exists(script)) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "blocked_no_stock_spec",
      reason = "No stock translation recipe is available."
    )
    return(result)
  }

  .assessment_log(log_file, if (fit) "Translating inputs and fitting model." else "Translating inputs.")
  stock <- .assessment_stock(source, do_fit = fit, silent = silent, tier = tier)
  stock_env <- stock$env
  sourced <- stock$error
  translated <- as.list(stock_env)
  result$obs <- translated$obs
  context <- .assessment_context(source, translated)
  result$background <- context$background
  result$comparison <- context$comparison
  elapsed <- as.numeric(difftime(Sys.time(), stock_env$fit_started, units = "secs"))
  if (inherits(sourced, "error")) {
    attempts <- sourced$initialization_attempts
    result$initialization_attempts <- attempts
    failure_reason <- conditionMessage(sourced)
    if (!is.null(attempts)) {
      failure_reason <- paste(failure_reason, paste(attempts$initialization,
        attempts$optimizer_message, attempts$notes, sep = ": ", collapse = " | "))
    }
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, paste0(stock_env$fit_stage, "_failed"),
      elapsed = elapsed, reason = failure_reason
    )
    return(result)
  }
  if (!fit) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "not_fitted",
      reason = "fit = FALSE."
    )
    return(result)
  }
  fitted <- stock_env$fit
  .assessment_log(log_file, "Model returned; checking diagnostics.")
  fitted <- .assessment_catch_reporting(fitted, translated$catch_reporting)
  fitted$assessment_provenance <- list(
    database_revision = database$commit, database_source = database$database_source,
    database_ref = database$database_ref, database_selection = database$database_selection,
    database_key = .assessment_cache_key(database), translation_tier = tier,
    tinyAM_revision = .database_head("."))
  result$fit <- fitted
  result$diagnostics <- .assessment_diagnostics(
    assessment_id, database,
    if (isTRUE(fitted$is_converged)) "converged" else "not_converged",
    fit = fitted, elapsed = elapsed
  )
  cache_dir <- .assessment_cache(assessment_id, tier, database)
  if (cache) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    .assessment_log(log_file, "Saving fitted object.")
    .assessment_save_rds(fitted, file.path(cache_dir, "fit.rds"))
  }
  if (!isTRUE(fitted$is_converged)) {
    return(result)
  }

  .assessment_log(log_file, "Building reference comparisons.")
  comparison <- tryCatch(.assessment_comparison(fitted, result), error = identity)
  if (inherits(comparison, "error")) {
    result$diagnostics$status <- "reference_failed"
    result$diagnostics$reason <- conditionMessage(comparison)
    return(result)
  }
  result[c("ref", "differences", "summary")] <- comparison[c("ref", "differences", "summary")]

  result
}

.failed_assessment_run <- function(assessment_id, database, error) {
  result <- list(
    assessment_id = assessment_id,
    translation_tier = database$translation_tier,
    source = NULL,
    obs = NULL,
    fit = NULL,
    ref = NULL,
    audit = NULL,
    diagnostics = .assessment_diagnostics(
      assessment_id, database, "runner_failed",
      reason = conditionMessage(error)
    ),
    differences = NULL,
    summary = NULL,
    background = NULL,
    dashboard_file = NULL
  )
  result
}

.assessment_worker <- function(database, repo_root, fit, cache = FALSE,
                               timeout = Inf, tier = NULL, isolated = TRUE, tick = NULL) {
  id <- database$assessment$assessment_id[[1]]
  if (!is.null(tick)) tick(amount = 0, message = paste(id, "started"))
  log_file <- NULL
  if (cache) {
    cache_dir <- .assessment_cache(id, tier, database)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(cache_dir, "run.log")
    writeLines(character(), log_file)
  }
  .assessment_log(log_file, paste(id, "started; timeout", timeout, "seconds."))
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  begun <- Sys.time()
  result <- tryCatch({
    # A subprocess deadline can also interrupt native optimizer code.
    if (is.finite(timeout)) {
      callr::r(function(id, database, repo_root, fit, cache, log_file, tier) {
        pkgload::load_all(repo_root, quiet = TRUE)
        runner <- new.env(parent = globalenv())
        sys.source(file.path(repo_root, "analysis", "comp_assessments", "R",
          "run_assessment.R"), envir = runner)
        runner$run_assessment(id, database = database, fit = fit,
          cache = cache, log_file = log_file, tier = tier)
      }, args = list(id, database, repo_root, fit, cache, log_file, tier),
      wd = repo_root, libpath = .libPaths(), timeout = timeout,
      stdout = if (is.null(log_file)) NULL else paste0(log_file, ".output"),
      stderr = "2>&1")
    } else if (isolated) {
      setwd(repo_root)
      pkgload::load_all(repo_root, quiet = TRUE)
      runner <- new.env(parent = globalenv())
      sys.source(file.path(repo_root, "analysis", "comp_assessments", "R",
        "run_assessment.R"), envir = runner)
      runner$run_assessment(id, database = database, fit = fit,
        cache = cache, log_file = log_file, tier = tier)
    } else {
      run_assessment(id, database = database, fit = fit,
        cache = cache, log_file = log_file, tier = tier)
    }
  }, error = function(e) {
    result <- .failed_assessment_run(id, database, e)
    result$diagnostics$elapsed_seconds <- as.numeric(difftime(Sys.time(), begun, units = "secs"))
    if (inherits(e, "callr_timeout_error")) {
      result$diagnostics$status <- "timed_out"
      result$diagnostics$reason <- paste("Time limit of", timeout,
        "seconds exceeded; model settings were not changed.")
    }
    result
  })
  if (cache) {
    .assessment_log(log_file, "Saving completed assessment checkpoint.")
    .assessment_save_rds(result, file.path(cache_dir, "result.rds"))
  }
  .assessment_log(log_file, paste(id, "finished:", result$diagnostics$status))
  if (!is.null(tick)) tick(message = paste(id, result$diagnostics$status))
  result
}

run_assessments <- function(assessment_ids = NULL, database = NULL,
                            fit = TRUE, cache = FALSE,
                            save_results = FALSE, timeout = Inf, progress = TRUE,
                            tier = "curated") {
  tier <- match.arg(tier, c("curated", "preliminary"))
  if (is.null(database)) database <- read_database()
  if (!is.logical(fit) || length(fit) != 1L || is.na(fit) ||
    !is.logical(cache) || length(cache) != 1L || is.na(cache) ||
    !is.logical(save_results) || length(save_results) != 1L || is.na(save_results)) {
    stop("fit, cache, and save_results must each be one TRUE/FALSE value.",
      call. = FALSE
    )
  }
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    cli::cli_abort("{.arg timeout} must be positive seconds or Inf.")
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    cli::cli_abort("{.arg progress} must be TRUE or FALSE.")
  }
  if (is.finite(timeout) && !requireNamespace("callr", quietly = TRUE)) {
    cli::cli_abort("Install callr to use an assessment time limit.")
  }

  database$translation_tier <- tier
  eligible <- assessment_translations(tier)$assessment_id
  eligible <- eligible[eligible %in% database$assessments$assessment_id]
  if (is.null(assessment_ids)) {
    assessment_ids <- eligible
  } else {
    if (!is.character(assessment_ids) || !length(assessment_ids) ||
      anyNA(assessment_ids) || any(!nzchar(assessment_ids)) ||
      anyDuplicated(assessment_ids)) {
      stop("assessment_ids must be a non-empty character vector without missing or duplicate values.",
        call. = FALSE
      )
    }
    unknown <- setdiff(assessment_ids, intersect(database$assessments$assessment_id,
      assessment_translations(tier)$assessment_id))
    if (length(unknown)) {
      stop("Unknown or unregistered assessment_id in the selected tier: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  if (!length(assessment_ids)) return(.assessment_result_tables(list(), save_results, tier))
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!requireNamespace("pkgload", quietly = TRUE)) cli::cli_abort("Install pkgload to load the current tinyAM checkout in assessment workers.")
  use_worker <- is.finite(timeout) || !inherits(future::plan(), "sequential")
  jobs <- lapply(assessment_ids, read_assessment, database = database)
  runs <- progressr::with_progress({
    tick <- progressr::progressor(steps = length(assessment_ids))
    furrr::future_map(jobs, .assessment_worker, repo_root = repo_root, fit = fit,
      cache = cache, timeout = timeout, tier = tier, isolated = use_worker, tick = tick,
      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))
  }, enable = progress, handlers = progressr::make_progression_handler(
    "assessment", interval = 0, enable = progress,
    reporter = list(update = function(config, state, progression, ...) {
      if (length(progression$message)) cli::cli_inform("{progression$message}")
    })
  ))
  names(runs) <- assessment_ids

  .assessment_result_tables(runs, save_results, tier)
}

# Rebuild complete-run tables without fitting or changing comparison eligibility.
.assessment_result_tables <- function(runs, save_results = FALSE, tier = "curated") {
  tiers <- unique(vapply(runs, function(x) {
    if (is.null(x$translation_tier)) tier else x$translation_tier
  }, character(1)))
  if (save_results && any(tiers != tier)) cli::cli_abort("Cannot save mixed translation tiers into one result set.")
  if (save_results && length(unique(vapply(runs, .assessment_result_key, character(1)))) > 1L) {
    cli::cli_abort("Cannot save mixed SAIO revisions into one current result set; keep historical batches separate.")
  }
  diagnostics <- if (length(runs)) {
    do.call(rbind, lapply(runs, `[[`, "diagnostics"))
  } else {
    data.frame()
  }
  summaries <- Filter(
    function(x) is.data.frame(x) && nrow(x),
    lapply(runs, `[[`, "summary")
  )
  comparison_summary <- if (length(summaries)) {
    do.call(rbind, summaries)
  } else {
    data.frame(
      assessment_id = character(), metric = character(), unit = character(),
      comparison_status = character(), definition = character(), reason = character(),
      n = integer(), mean_absolute_difference = numeric(),
      mean_absolute_percent_difference = numeric(),
      median_absolute_percent_difference = numeric(),
      terminal_year = integer(), terminal_mean_percent_difference = numeric(),
      trend_correlation = numeric()
    )
  }
  rownames(diagnostics) <- NULL
  rownames(comparison_summary) <- NULL

  # Assemble annual population comparisons across stocks.
  population_metrics <- c("ssb", "recruitment", "abundance", "F_bar")
  series <- Filter(Negate(is.null), lapply(runs, function(x) {
    if (!is.data.frame(x$differences)) return(NULL)
    historical_years <- x$fit$dat$years[!x$fit$dat$is_proj]
    d <- subset(
      x$differences,
      metric %in% population_metrics & year %in% historical_years
    )
    if (!nrow(d)) return(NULL)
    data.frame(assessment_id = x$assessment_id, d)
  }))

  population_timeseries <- if (length(series)) {
    do.call(rbind, series)
  } else {
    data.frame(
      assessment_id = character(),
      metric = character(),
      year = integer(),
      age = character(),
      source = numeric(),
      tinyAM = numeric(),
      absolute_difference = numeric(),
      percent_difference = numeric(),
      unit = character(),
      comparison_status = character(),
      definition = character(),
      reason = character()
    )
  }
  rownames(population_timeseries) <- NULL

  if (save_results) {
    results_dir <- .assessment_results(tier)
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(diagnostics, file.path(results_dir, "fit_diagnostics.csv"),
      row.names = FALSE, na = ""
    )
    write.csv(comparison_summary,
      file.path(results_dir, "comparison_summary.csv"),
      row.names = FALSE, na = ""
    )
    write.csv(
      population_timeseries,
      file.path(results_dir, "population_timeseries.csv"),
      row.names = FALSE,
      na = ""
    )
  }
  attr(runs, "diagnostics") <- diagnostics
  attr(runs, "comparison_summary") <- comparison_summary
  attr(runs, "population_timeseries") <- population_timeseries
  runs
}
