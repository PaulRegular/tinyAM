# Run from the tinyAM root. Set TINYAM_PLD_SIMSURVEY to the pinned source
# checkout when generating the whole-area study.
Sys.setenv(OMP_NUM_THREADS = 1L, OPENBLAS_NUM_THREADS = 1L, MKL_NUM_THREADS = 1L)
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1L] else "report"
study <- "analysis/formula_expansion"
results <- file.path(study, "results")
out <- if (length(args) >= 2L) args[2L] else file.path(results,
  if (phase == "controlled") "pld-controlled" else "pld-whole-area")
if (!phase %in% c("controlled", "audit", "spatial", "report")) {
  stop("Use simulate_pld.R {controlled|audit|spatial|report} [results directory].")
}
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out, winslash = "/")
source(file.path(study, "pld_simple.R"))
if (phase == "report") {
  controlled <- if (length(args) >= 3L) args[3L] else file.path(results, "pld-simple-refined")
  path <- function(file) normalizePath(file, winslash = "/", mustWork = TRUE)
  pilot <- function(name) {
    file <- file.path(results, name, "summary.rds")
    if (file.exists(file)) path(file) else ""
  }
  spatial <- file.path(out, "spatial-summary.rds")
  if (!file.exists(spatial)) spatial <- file.path(out, "spatial-audit.rds")
  report_inputs <- list(result = path(file.path(controlled, "summary.rds")),
    spatial_result = if (file.exists(spatial)) path(spatial) else "",
    default_result = pilot("pld-simple"), tight_result = pilot("pld-simple-tight"))
  html <- rmarkdown::render(file.path(study, "pld_validation_report.Rmd"),
    params = report_inputs, output_dir = out, quiet = TRUE)
  cat("Report:", html, "\n")
  quit(save = "no")
}

pkg <- Sys.getenv("TINYAM_PLD_PACKAGE_SOURCE", ".")
pkgload::load_all(pkg, quiet = TRUE, helpers = FALSE)
source(file.path(study, "pld_validation_fit.R"))
data.table::setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(warn = 1)

checkpoint <- function(file, evaluate) {
  file <- file.path(out, file)
  if (file.exists(file)) return(readRDS(file))
  warnings <- character()
  began <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(evaluate(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  ok <- !inherits(value, "error")
  result <- list(ok = ok, value = if (ok) value else NULL,
    error = if (ok) "" else conditionMessage(value), warnings = unique(warnings),
    elapsed = proc.time()[["elapsed"]] - began)
  temporary <- paste0(file, ".partial")
  saveRDS(result, temporary)
  if (!file.rename(temporary, file)) stop("Failed to checkpoint ", file)
  cat(sprintf("%s %.2fs %s\n", basename(file), result$elapsed,
    if (ok) "saved" else paste("error", result$error)))
  flush.console()
  invisible(result)
}


if (phase == "controlled") {
  optimizer <- if (length(args) >= 3L) match.arg(args[3L], c("default", "tight", "newton")) else "newton"
  files <- sort(c("DESCRIPTION", "NAMESPACE",
    list.files(file.path(pkg, "R"), "[.]R$", full.names = TRUE)))
  files[files %in% c("DESCRIPTION", "NAMESPACE")] <- file.path(pkg, files[files %in% c("DESCRIPTION", "NAMESPACE")])
  sources <- c(tools::md5sum(files), tools::md5sum(file.path(study,
    c("pld_simple.R", "pld_validation_fit.R", "simulate_pld.R"))))
  manifest <- list(sources = sources, optimizer = optimizer, R = as.character(getRversion()),
    RTMB = as.character(utils::packageVersion("RTMB")))
  stopifnot(identical(checkpoint("manifest.rds", function() manifest)$value, manifest))
  calibration <- checkpoint("calibration.rds", pld_simple_calibration)$value
  stopifnot(nrow(calibration$values) == 2L)
  print(calibration$values)

  stages <- c("exact", "finite", "pulse", "constant_m")
  planned <- c(exact = 32L, finite = 32L, pulse = 16L, constant_m = 16L)
  arms <- list(exact = c("full", "positive", "lognormal"),
    finite = c("full", "lognormal"), pulse = c("full", "positive"),
    constant_m = c("full", "positive"))
  run_stage <- function(ids, condition) {
    fit_stage <- if (condition %in% c("pulse", "constant_m")) "mortality" else condition
    lapply(ids, function(id) {
      truth <- pld_simple_truth(id, pulse = condition != "constant_m")
      scheme <- if (condition == "finite") "finite" else if (truth$pulse) "exact" else "null"
      input <- checkpoint(sprintf("observations-%s-%02d.rds", scheme, id),
        function() pld_simple_observations(truth, calibration, finite = condition == "finite"))
      result <- lapply(arms[[condition]], function(arm) {
        record <- checkpoint(sprintf("fit-%s-%02d-%s.rds", condition, id, arm), function() {
          if (!input$ok) stop(input$error)
          pld_simple_fit(input$value, truth, calibration, arm, fit_stage, optimizer)
        })
        list(stage = condition, id = id, arm = arm, converged = isTRUE(record$value$converged),
          numerical = record$value$numerical, gradient = record$value$gradient,
          study_newton = isTRUE(record$value$study_newton),
          error = record$error, warnings = record$warnings, values = record$value$values)
      })
      list(id = id, zeros = input$value$zeros, fits = stats::setNames(result, arms[[condition]]))
    })
  }

  paired_success <- function(histories, reference) {
    sum(vapply(histories, function(record) isTRUE(record$fits$full$converged) &&
      isTRUE(record$fits[[reference]]$converged), logical(1L)))
  }

  histories <- stats::setNames(lapply(stages, function(x) run_stage(1:8, x)), stages)
  gate <- do.call(rbind, lapply(stages, function(x) {
    references <- setdiff(arms[[x]], "full")
    data.frame(stage = x, reference = references,
      pairs = vapply(references, function(arm) paired_success(histories[[x]], arm), integer(1L)))
  }))
  print(gate, row.names = FALSE)
  for (condition in stages) {
    permitted <- all(gate$pairs[gate$stage == condition] >= 7L)
    if (condition %in% c("pulse", "constant_m")) {
      permitted <- all(gate$pairs[gate$stage %in% c("pulse", "constant_m")] >= 7L)
    }
    if (permitted) {
      histories[[condition]] <- c(histories[[condition]], run_stage(9:planned[[condition]], condition))
    }
  }
  result <- list(calibration = calibration, stages = histories, gate = gate, planned = planned,
    provenance = manifest)
  path <- file.path(out, "summary.rds")
  if (!file.exists(path)) saveRDS(result, path)
  cat("Summary:", path, "\n")

}

if (phase != "controlled") {
  source(file.path(study, "pld_validation.R"))
  source(file.path(study, "observation_error_sim.R"))
  source(file.path(study, "pld_fishery_engine.R"))
  sim <- normalizePath(Sys.getenv("TINYAM_PLD_SIMSURVEY"), winslash = "/", mustWork = TRUE)
  revision <- suppressWarnings(system2("git", c("-C", shQuote(sim), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE))
  if (!identical(revision, pld_simsurvey_revision) &&
      !identical(basename(sim), paste0("SimSurvey-", pld_simsurvey_revision))) {
    stop("Set TINYAM_PLD_SIMSURVEY to the SimSurvey source at ", pld_simsurvey_revision)
  }
  pkgload::load_all(sim, quiet = TRUE, helpers = FALSE)
  pld_check_allocate_lengths()
  freeze <- function(file, value) {
    full <- file.path(out, file)
    if (file.exists(full) && !identical(readRDS(full), value)) {
      stop("Study sources differ from ", file, ". Use a new results directory.")
    }
    if (!file.exists(full)) saveRDS(value, full)
  }
  freeze("spatial-data-manifest.rds", list(simsurvey = pld_simsurvey_revision,
    sources = tools::md5sum(sort(list.files(file.path(sim, "R"), "[.]R$",
      full.names = TRUE))), simulation = tools::md5sum(file.path(study,
        c("pld_validation.R", "pld_fishery_engine.R", "observation_error_sim.R"))),
    R = as.character(getRversion()), grid = pld_grid))
  if (phase == "spatial") {
    freeze("spatial-fit-manifest.rds", list(package = tools::md5sum(sort(c(
      file.path(pkg, c("DESCRIPTION", "NAMESPACE")), list.files(file.path(pkg, "R"),
        "[.]R$", full.names = TRUE)))),
      settings = tools::md5sum(file.path(study, "pld_validation_fit.R")),
      RTMB = as.character(packageVersion("RTMB"))))
  }
  arm_names <- c("Lognormal", "PLD")
  run_history <- function(id, plus_ages) {
    draw <- checkpoint(sprintf("spatial-%02d.rds", id), function() pld_simsurvey(id))
    if (!draw$ok) return(list(id = id, error = draw$error, warnings = draw$warnings))
    observed <- do.call(rbind, lapply(names(draw$value$data)[1:2], function(part) {
      d <- draw$value$data[[part]][c("year", "age", "obs")]
      d$component <- part
      d$id <- id
      d
    }))
    fits <- setNames(lapply(plus_ages, function(age) {
      arms <- setNames(lapply(arm_names, function(arm) {
        checkpoint(sprintf("spatial-fit-%02d-age%02d-%s.rds", id, age, arm),
          function() pld_spatial_fit(draw$value, age, arm))
      }), arm_names)
      if (all(vapply(arms, "[[", logical(1L), "ok")) &&
          !identical(arms[[1L]]$value$merged_observations, arms[[2L]]$value$merged_observations)) {
        stop("The paired fits saw different observation rows for history ", id, " at age ", age)
      }
      arms
    }), as.character(plus_ages))
    list(id = id, warnings = draw$warnings, audit = pld_mortality_audit(draw$value),
      reference = draw$value$reference, schedule = draw$value$schedule,
      observed = observed, fits = fits)
  }
  passed_pair <- function(history, age) {
    fits <- history$fits[[as.character(age)]]
    length(fits) == 2L && all(vapply(fits, function(result) isTRUE(result$value$converged), logical(1L)))
  }
  histories <- lapply(1:2, run_history, plus_ages = integer())
  valid_data <- function(history) !is.null(history$audit) && isTRUE(history$audit$adequate)
  gate <- "Two-history accounting and observation check"
  checked <- all(vapply(histories, function(history) !is.null(history$audit), logical(1L)))
  if (phase == "spatial" && checked) {
    histories <- lapply(1:8, run_history, plus_ages = c(20L, 14L))
    n_pairs <- sum(vapply(histories, passed_pair, logical(1L), age = 20L))
    n_adequate <- sum(vapply(histories, valid_data, logical(1L)))
    gate <- sprintf("The pilot had %d of 8 successful age-20+ pairs and %d of 8 histories within the M tolerances.",
      n_pairs, n_adequate)
    if (n_pairs >= 7L && n_adequate == 8L) {
      histories <- c(histories, lapply(9:24, run_history, plus_ages = 20L))
      gate <- paste(gate, "Expanded to 24 histories.")
    } else {
      gate <- paste(gate, "Expansion stopped at eight histories.")
    }
  }
  cat(gate, "\n")
  name <- if (phase == "audit") "spatial-audit.rds" else "spatial-summary.rds"
  saveRDS(list(phase = phase, gate = gate, histories = histories, pulse = pld_pulse,
    M = pld_m, older_window = pld_older_window), file.path(out, name))
}
