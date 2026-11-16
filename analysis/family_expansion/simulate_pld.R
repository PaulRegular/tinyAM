# Run from the tinyAM root. The comparison phase uses the original saved
# observations. Set TINYAM_PLD_SIMSURVEY only when generating spatial histories.
Sys.setenv(OMP_NUM_THREADS = 1L, OPENBLAS_NUM_THREADS = 1L, MKL_NUM_THREADS = 1L)
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1L] else "report"
study <- "analysis/family_expansion"
results <- file.path(study, "results")
archive <- "analysis/formula_expansion/results"
out <- if (length(args) >= 2L) args[2L] else file.path(results, switch(phase,
  controlled = "pld-controlled", audit = "pld-whole-area", spatial = "pld-whole-area",
  "pld-family-comparison"))
if (!phase %in% c("controlled", "audit", "spatial", "comparison", "report")) {
  stop("Use simulate_pld.R {controlled|audit|spatial|comparison|report} [results directory].")
}
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out, winslash = "/")
if (phase == "report") {
  path <- function(name) normalizePath(file.path(out, name), winslash = "/", mustWork = TRUE)
  html <- rmarkdown::render(file.path(study, "pld_validation_report.Rmd"),
    params = list(result = path("summary.rds"), spatial_result = path("spatial-summary.rds")),
    output_dir = out, quiet = TRUE)
  cat("Report:", html, "\n")
  quit(save = "no")
}

pkg <- Sys.getenv("TINYAM_PLD_PACKAGE_SOURCE", ".")
pkgload::load_all(pkg, quiet = TRUE, helpers = FALSE)
source(file.path(study, "pld_simple.R"))
source(file.path(study, "pld_validation_fit.R"))
source(file.path(study, "pld_validation.R"))
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

freeze <- function(name, value) {
  file <- file.path(out, name)
  if (file.exists(file) && !identical(readRDS(file), value)) {
    stop("Study sources or inputs differ from ", name, ". Use a new results directory.")
  }
  if (!file.exists(file)) saveRDS(value, file)
  invisible(value)
}

package_files <- sort(c(file.path(pkg, c("DESCRIPTION", "NAMESPACE")),
  list.files(file.path(pkg, "R"), "[.]R$", full.names = TRUE)))
observations_code <- c(file.path(study,
  c("pld_validation.R", "pld_fishery_engine.R")),
  "analysis/formula_expansion/observation_error_sim.R")
fit_code <- file.path(study, c("pld_simple.R", "pld_validation_fit.R", "simulate_pld.R"))
manifest <- list(package = tools::md5sum(package_files),
  study = tools::md5sum(c(fit_code, observations_code)),
  R = as.character(getRversion()), RTMB = as.character(utils::packageVersion("RTMB")))
arms <- c("lognormal", "full", "gengamma")
planned <- c(exact = 32L, finite = 32L, pulse = 16L, constant_m = 16L)
observation_scheme <- c(exact = "exact", finite = "finite", pulse = "exact", constant_m = "null")

verify_rows <- function(fits) {
  rows <- lapply(Filter(function(fit) fit$ok, fits), function(fit) fit$value$merged_observations)
  if (length(rows) > 1L && !all(vapply(rows[-1L], identical, logical(1L), y = rows[[1L]]))) {
    stop("The fits did not use the same prepared observation values.")
  }
}

controlled <- function(read_observations, calibration, optimizer = "newton") {
  stopifnot(nrow(calibration$values) == 2L)
  print(calibration$values)
  histories <- lapply(names(planned), function(condition) {
    fit_stage <- if (condition %in% c("pulse", "constant_m")) "mortality" else condition
    comparisons <- if (condition == "finite") arms else c(arms, "positive")
    lapply(seq_len(planned[[condition]]), function(id) {
      truth <- pld_simple_truth(id, pulse = condition != "constant_m")
      input <- read_observations(id, truth, observation_scheme[[condition]])
      fits <- setNames(lapply(comparisons, function(arm) {
        checkpoint(sprintf("fit-%s-%02d-%s.rds", condition, id, arm), function() {
          if (!input$ok) stop(input$error)
          pld_simple_fit(input$value, truth, calibration, arm, fit_stage, optimizer)
        })
      }), comparisons)
      verify_rows(fits)
      fits <- lapply(names(fits), function(arm) {
        record <- fits[[arm]]
        list(stage = condition, id = id, arm = arm, converged = isTRUE(record$value$converged),
          numerical = record$value$numerical, gradient = record$value$gradient,
          study_newton = isTRUE(record$value$study_newton),
          error = record$error, warnings = record$warnings, elapsed = record$elapsed,
          values = record$value$values, parameters = record$value$parameters)
      })
      list(id = id, zeros = input$value$zeros, fits = setNames(fits, comparisons))
    })
  })
  summary <- list(calibration = calibration, stages = setNames(histories, names(planned)),
    planned = planned, provenance = c(manifest, list(optimizer = optimizer)))
  saveRDS(summary, file.path(out, "summary.rds"))
}

spatial_history <- function(id, draw, plus_ages) {
  if (!draw$ok) return(list(id = id, error = draw$error, warnings = draw$warnings))
  observed <- do.call(rbind, lapply(names(draw$value$data)[1:2], function(part) {
    d <- draw$value$data[[part]][c("year", "age", "obs")]
    d$component <- part
    d$id <- id
    d
  }))
  fits <- setNames(lapply(plus_ages, function(age) {
    comparison <- setNames(lapply(arms, function(arm) {
      checkpoint(sprintf("spatial-fit-%02d-age%02d-%s.rds", id, age, arm),
        function() pld_spatial_fit(draw$value, age, arm))
    }), arms)
    verify_rows(comparison)
    comparison
  }), as.character(plus_ages))
  list(id = id, warnings = draw$warnings, audit = pld_mortality_audit(draw$value),
    reference = draw$value$reference, schedule = draw$value$schedule, observed = observed, fits = fits)
}

save_spatial <- function(histories, gate, name = "spatial-summary.rds") {
  cat(gate, "\n")
  saveRDS(list(gate = gate, histories = histories, pulse = pld_pulse, M = pld_m,
    older_window = pld_older_window, provenance = manifest), file.path(out, name))
}

if (phase == "comparison") {
  controlled_dir <- if (length(args) >= 3L) args[3L] else file.path(archive, "pld-simple-refined")
  spatial_dir <- if (length(args) >= 4L) args[4L] else file.path(archive, "pld-whole-area")
  calibration_path <- file.path(controlled_dir, "calibration.rds")
  controlled_file <- function(scheme, id) {
    file.path(controlled_dir, sprintf("observations-%s-%02d.rds", scheme, id))
  }
  spatial_files <- file.path(spatial_dir, sprintf("spatial-%02d.rds", 1:8))
  required <- unique(c(calibration_path, unlist(lapply(names(planned), function(stage) {
    controlled_file(observation_scheme[[stage]], seq_len(planned[[stage]]))
  })), spatial_files))
  for (file in required) if (!file.exists(file)) stop("Missing frozen input: ", file)
  freeze("comparison-manifest.rds", c(manifest, list(inputs = tools::md5sum(required))))
  calibration <- readRDS(calibration_path)$value
  controlled(function(id, truth, scheme) readRDS(controlled_file(scheme, id)), calibration)
  histories <- lapply(1:8, function(id) spatial_history(id, readRDS(spatial_files[id]), c(20L, 14L)))
  save_spatial(histories, paste("The original whole-area pilot stopped at eight histories.",
    "The three-family re-fit uses those same histories."))
}

if (phase == "controlled") {
  optimizer <- if (length(args) >= 3L) match.arg(args[3L], c("default", "tight", "newton")) else "newton"
  freeze("controlled-manifest.rds", c(manifest, list(optimizer = optimizer)))
  calibration <- checkpoint("calibration.rds", pld_simple_calibration)$value
  controlled(function(id, truth, scheme) {
    checkpoint(sprintf("observations-%s-%02d.rds", scheme, id),
      function() pld_simple_observations(truth, calibration, finite = scheme == "finite"))
  }, calibration, optimizer)
}

if (phase %in% c("audit", "spatial")) {
  source(observations_code[3L])
  source(observations_code[2L])
  sim <- normalizePath(Sys.getenv("TINYAM_PLD_SIMSURVEY"), winslash = "/", mustWork = TRUE)
  revision <- suppressWarnings(system2("git", c("-C", shQuote(sim), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE))
  if (!identical(revision, pld_simsurvey_revision) &&
      !identical(basename(sim), paste0("SimSurvey-", pld_simsurvey_revision))) {
    stop("Set TINYAM_PLD_SIMSURVEY to the SimSurvey source at ", pld_simsurvey_revision)
  }
  pkgload::load_all(sim, quiet = TRUE, helpers = FALSE)
  pld_check_allocate_lengths()
  freeze("spatial-data-manifest.rds", list(simsurvey = pld_simsurvey_revision,
    sources = tools::md5sum(sort(list.files(file.path(sim, "R"), "[.]R$", full.names = TRUE))),
    simulation = tools::md5sum(observations_code),
    R = as.character(getRversion()), grid = pld_grid))
  if (phase == "spatial") freeze("spatial-fit-manifest.rds", manifest)
  run_history <- function(id, plus_ages) {
    draw <- checkpoint(sprintf("spatial-%02d.rds", id), function() pld_simsurvey(id))
    spatial_history(id, draw, plus_ages)
  }
  histories <- lapply(1:2, run_history, plus_ages = integer())
  checked <- all(vapply(histories, function(history) !is.null(history$audit), logical(1L)))
  if (phase == "audit") {
    save_spatial(histories, "Two-history mortality and observation check", "spatial-audit.rds")
  } else if (checked) {
    histories <- lapply(1:8, run_history, plus_ages = c(20L, 14L))
    save_spatial(histories, "Eight-history whole-area pilot.")
  } else {
    stop("The initial two spatial histories did not both generate observations.")
  }
}
