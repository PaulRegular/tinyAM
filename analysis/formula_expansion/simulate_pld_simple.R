# Run the first eight histories per stage before evaluating the expansion gate.
arguments <- commandArgs(trailingOnly = TRUE)
phase <- if (length(arguments)) arguments[[1L]] else "run"
if (!phase %in% c("run", "report")) stop("Use simulate_pld_simple.R [run|report] [results directory].")
study <- "analysis/formula_expansion"
out <- if (length(arguments) >= 2L) arguments[[2L]] else file.path(study, "results/pld-simple")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out, winslash = "/")
source(file.path(study, "pld_simple.R"))
if (phase == "report") {
  result <- file.path(out, "summary.rds")
  default_result <- if (length(arguments) >= 3L) {
    normalizePath(arguments[[3L]], winslash = "/", mustWork = TRUE)
  } else ""
  tight_result <- if (length(arguments) >= 4L) {
    normalizePath(arguments[[4L]], winslash = "/", mustWork = TRUE)
  } else ""
  stopifnot(file.exists(result))
  html <- rmarkdown::render(file.path(study, "pld_simple_report.Rmd"),
    params = list(result = result, default_result = default_result, tight_result = tight_result),
    output_dir = out, quiet = TRUE)
  cat("Report:", html, "\n")
  quit(save = "no")
}

pkg <- Sys.getenv("TINYAM_PLD_SIMPLE_SOURCE", ".")
pkgload::load_all(pkg, quiet = TRUE, helpers = FALSE)
source(file.path(study, "pld_validation_fit.R"))
Sys.setenv(OMP_NUM_THREADS = 1L, OPENBLAS_NUM_THREADS = 1L, MKL_NUM_THREADS = 1L)
data.table::setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(warn = 1)
optimizers <- c("default", "tight", "newton")
optimizer <- if (length(arguments) >= 3L) match.arg(arguments[[3L]], optimizers) else optimizers[1L]

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

files <- sort(c("DESCRIPTION", "NAMESPACE",
  list.files(file.path(pkg, "R"), "[.]R$", full.names = TRUE)))
files[files %in% c("DESCRIPTION", "NAMESPACE")] <- file.path(pkg, files[files %in% c("DESCRIPTION", "NAMESPACE")])
sources <- c(tools::md5sum(files), tools::md5sum(file.path(study,
  c("pld_simple.R", "pld_validation_fit.R", "simulate_pld_simple.R"))))
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
