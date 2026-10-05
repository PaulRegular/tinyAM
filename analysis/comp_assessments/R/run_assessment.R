.assessment_root <- file.path("analysis", "comp_assessments")
for (.assessment_helper in c(
  "read_database.R", "database_to_tam_obs.R", "database_to_tam_ref.R",
  "audit_assumptions.R"
)) {
  sys.source(file.path(.assessment_root, "R", .assessment_helper),
             envir = environment())
}
rm(.assessment_helper)

.assessment_fit_call <- function(args) {
  args$obs <- quote(obs)
  if (!is.null(args$start_par)) args$start_par <- quote(start_par)
  as.call(c(list(quote(fit_tam)), args))
}

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
                            dashboard = FALSE, cache = FALSE, silent = TRUE) {
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
  sourced <- tryCatch({
    sys.source(script, envir = stock_env)
    NULL
  }, error = identity)
  if (inherits(sourced, "error")) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "translation_script_failed",
      reason = conditionMessage(sourced)
    )
    return(result)
  }
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
                     ages = translated$ages, silent = silent), translated$settings)
  if (!is.null(translated$start_par)) {
    fit_args$start_par <- translated$start_par
  } else if (!is.null(translated$warm_start_settings)) {
    warm_settings <- utils::modifyList(
      translated$settings, translated$warm_start_settings
    )
    warm_args <- c(list(obs = translated$obs, years = translated$years,
                        ages = translated$ages, silent = silent), warm_settings)
    warm_fit <- tryCatch(do.call(tinyAM::fit_tam, warm_args), error = identity)
    if (inherits(warm_fit, "error") || !isTRUE(warm_fit$is_converged)) {
      result$diagnostics <- .assessment_diagnostics(
        assessment_id, database, "warm_start_failed",
        elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")),
        reason = if (inherits(warm_fit, "error")) {
          conditionMessage(warm_fit)
        } else {
          "The preliminary fit did not converge."
        }
      )
      return(result)
    }
    fit_args$start_par <- as.list(warm_fit$sdrep, "Estimate")
  }
  fitted <- tryCatch(do.call(tinyAM::fit_tam, fit_args), error = identity)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fitted, "error")) {
    result$diagnostics <- .assessment_diagnostics(
      assessment_id, database, "fit_failed", elapsed = elapsed,
      reason = conditionMessage(fitted)
    )
    return(result)
  }
  fitted$call <- .assessment_fit_call(fit_args)
  result$fit <- fitted
  result$diagnostics <- .assessment_diagnostics(
    assessment_id, database,
    if (isTRUE(fitted$is_converged)) "converged" else "not_converged",
    fit = fitted, elapsed = elapsed
  )
  cache_dir <- file.path(.assessment_root, "results", "cache", assessment_id)
  if (cache) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(fitted, file.path(cache_dir, "fit.rds"))
  }
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

.failed_assessment_run <- function(assessment_id, database, error) {
  result <- list(
    assessment_id = assessment_id,
    source = NULL,
    obs = NULL,
    settings = NULL,
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

.parallel_assessment_worker <- function(assessment_id, database, repo_root, fit) {
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  tryCatch(
    {
      setwd(repo_root)
      pkgload::load_all(repo_root, quiet = TRUE)
      runner_env <- new.env(parent = globalenv())
      sys.source(file.path(repo_root, "analysis", "comp_assessments", "R",
                           "run_assessment.R"), envir = runner_env)
      runner_env$run_assessment(
        assessment_id, database = database, fit = fit,
        dashboard = FALSE, cache = FALSE
      )
    },
    error = function(e) .failed_assessment_run(assessment_id, database, e)
  )
}

run_assessments <- function(assessment_ids = NULL, database = NULL,
                            parallel = FALSE, workers = 4L, fit = TRUE,
                            dashboard = FALSE, cache = FALSE,
                            save_results = FALSE) {
  if (is.null(database)) database <- read_database()
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel) ||
      !is.logical(fit) || length(fit) != 1L || is.na(fit) ||
      !is.logical(dashboard) || length(dashboard) != 1L || is.na(dashboard) ||
      !is.logical(cache) || length(cache) != 1L || is.na(cache) ||
      !is.logical(save_results) || length(save_results) != 1L || is.na(save_results)) {
    stop("parallel, fit, dashboard, cache, and save_results must each be one TRUE/FALSE value.",
         call. = FALSE)
  }
  if (length(workers) != 1L || is.na(workers) || !is.numeric(workers) ||
      !is.finite(workers) || workers < 1 || workers != as.integer(workers)) {
    stop("workers must be one positive integer.", call. = FALSE)
  }

  eligible <- database$assessments$assessment_id[
    !is.na(database$assessments$is_current) & database$assessments$is_current &
      !is.na(database$assessments$is_applied) & database$assessments$is_applied
  ]
  if (is.null(assessment_ids)) {
    assessment_ids <- eligible
  } else {
    if (!is.character(assessment_ids) || !length(assessment_ids) ||
        anyNA(assessment_ids) || any(!nzchar(assessment_ids)) ||
        anyDuplicated(assessment_ids)) {
      stop("assessment_ids must be a non-empty character vector without missing or duplicate values.",
           call. = FALSE)
    }
    unknown <- setdiff(assessment_ids, database$assessments$assessment_id)
    if (length(unknown)) {
      stop("Unknown assessment_id: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  safe_run <- function(id) {
    tryCatch(
      run_assessment(id, database = database, fit = fit,
                     dashboard = FALSE, cache = FALSE),
      error = function(e) .failed_assessment_run(id, database, e)
    )
  }
  if (parallel) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install pkgload to run parallel workers from the current tinyAM checkout.",
           call. = FALSE)
    }
    repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = as.integer(workers))
    runs <- furrr::future_map(
      assessment_ids, .parallel_assessment_worker,
      database = database, repo_root = repo_root, fit = fit,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    runs <- lapply(assessment_ids, safe_run)
  }
  names(runs) <- assessment_ids

  for (id in assessment_ids) {
    result <- runs[[id]]
    if (!is.null(result$fit) && cache) {
      cache_dir <- file.path(.assessment_root, "results", "cache", id)
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(result$fit, file.path(cache_dir, "fit.rds"))
    }
    if (dashboard && !is.null(result$fit) && isTRUE(result$fit$is_converged) &&
        !is.null(result$ref)) {
      result$dashboard_file <- if (cache) {
        file.path(.assessment_root, "results", "cache", id, "dashboard.html")
      } else {
        tempfile(pattern = paste0(id, "_"), fileext = ".html")
      }
      tinyAM::vis_tam(
        model_list = list(Accepted = result$ref, tinyAM = result$fit),
        background = result$background,
        output_file = result$dashboard_file,
        open_file = FALSE
      )
      runs[[id]] <- result
    }
  }

  diagnostics <- if (length(runs)) {
    do.call(rbind, lapply(runs, `[[`, "diagnostics"))
  } else {
    data.frame()
  }
  summaries <- Filter(function(x) is.data.frame(x) && nrow(x),
                      lapply(runs, `[[`, "summary"))
  comparison_summary <- if (length(summaries)) {
    do.call(rbind, summaries)
  } else {
    data.frame(assessment_id = character(), metric = character(), n = integer(),
               mean_absolute_percent_difference = numeric(),
               median_absolute_percent_difference = numeric(),
               terminal_year = integer(), terminal_mean_percent_difference = numeric(),
               trend_correlation = numeric())
  }
  rownames(diagnostics) <- NULL
  rownames(comparison_summary) <- NULL
  if (save_results) {
    results_dir <- file.path(.assessment_root, "results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(diagnostics, file.path(results_dir, "fit_diagnostics.csv"),
              row.names = FALSE, na = "")
    write.csv(comparison_summary,
              file.path(results_dir, "comparison_summary.csv"),
              row.names = FALSE, na = "")
  }
  attr(runs, "diagnostics") <- diagnostics
  attr(runs, "comparison_summary") <- comparison_summary
  runs
}
