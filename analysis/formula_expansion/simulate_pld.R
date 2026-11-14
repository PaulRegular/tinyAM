# Run from the tinyAM repository root. Set TINYAM_PLD_SIMSURVEY to an
# extracted source tree or local checkout when generating observations.
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(OMP_NUM_THREADS = 1L, OPENBLAS_NUM_THREADS = 1L, MKL_NUM_THREADS = 1L)
phase <- if (length(args)) args[1L] else ""
n <- if (length(args) >= 2L) as.integer(args[2L]) else 24L
grid <- if (length(args) >= 3L) as.integer(args[3L]) else 40L
out <- if (length(args) >= 4L) args[4L] else file.path(
  "analysis/formula_expansion/results", paste0("pld-grid", grid, "-n", n))
if (!phase %in% c("controlled", "sampling", "fit", "report", "sampling_fit") ||
    is.na(n) || !n %in% 1:24 || is.na(grid) || grid < 20L) {
  stop("Use simulate_pld.R {controlled|sampling|fit|report|sampling_fit} [replicates] [grid width] [output directory].")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/")
study <- "analysis/formula_expansion"
source(file.path(study, "pld_validation.R"))
data.table::setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(warn = 1)

pld_checkpoint <- function(file, evaluate) {
  if (file.exists(file)) return(invisible(NULL))
  saw <- character()
  start <- proc.time()[["elapsed"]]
  cpu_start <- sum(proc.time()[c("user.self", "sys.self")])
  value <- tryCatch(withCallingHandlers(evaluate(),
    warning = function(w) {
      saw <<- c(saw, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  success <- !inherits(value, "error")
  result <- list(ok = success, value = if (success) value else NULL,
    error = if (!success) conditionMessage(value) else "", warnings = saw,
    elapsed = proc.time()[["elapsed"]] - start,
    cpu_seconds = sum(proc.time()[c("user.self", "sys.self")]) - cpu_start)
  partial <- paste0(file, ".partial")
  saveRDS(result, partial)
  if (!file.rename(partial, file)) stop("Could not save checkpoint ", file)
  cat(sprintf("%s %.1fs (%.1fs CPU) %s\n", basename(file), result$elapsed, result$cpu_seconds,
    if (result$ok) "saved" else paste("error", result$error)))
  flush.console()
  invisible(NULL)
}

pld_record_sources <- function(path, value) {
  if (file.exists(path)) {
    if (!identical(readRDS(path), value)) {
      stop("Study inputs or code changed. Choose a new results directory.")
    }
  } else {
    saveRDS(value, path)
  }
}

if (phase %in% c("controlled", "sampling", "sampling_fit")) {
  sim <- normalizePath(Sys.getenv("TINYAM_PLD_SIMSURVEY"), mustWork = TRUE, winslash = "/")
  head <- suppressWarnings(system2("git", c("-C", shQuote(sim), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE))
  base <- basename(sim)
  document_only <- "bb4a06d8fef8f4eb74ecfa179d6dc831f589e5fe"
  revisions <- c(pld_simsurvey_revision, document_only)
  if (length(head) == 1L && head %in% revisions) {
    revision <- head
  } else if (base %in% paste0("SimSurvey-", revisions)) {
    revision <- sub("^SimSurvey-", "", base)
  } else {
    stop("Use SimSurvey ", pld_simsurvey_revision, " or ", document_only,
      " (identical simulation code). Set TINYAM_PLD_SIMSURVEY to the source tree.")
  }
  pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
  pkgload::load_all(sim, quiet = TRUE, helpers = FALSE)
  generating_sources <- c(file.path(study, "simulate_pld.R"),
    file.path(study, "pld_validation.R"))
  if (phase != "controlled") {
    shared_source <- file.path(study, "observation_error_sim.R")
    source(shared_source)
    generating_sources <- c(generating_sources, shared_source)
    pld_check_allocate_lengths()
  }
  simsurvey_sources <- unname(tools::md5sum(sort(list.files(file.path(sim, "R"),
    "[.]R$", full.names = TRUE))))
  fishery_cache <- Sys.getenv("TINYAM_PLD_FISHERY_CACHE")
  fishery_provenance <- NULL
  if (phase != "controlled" && nzchar(fishery_cache)) {
    fishery_cache <- normalizePath(fishery_cache, winslash = "/", mustWork = TRUE)
    pld_check_fishery_cache(fishery_cache, grid, n, revision, simsurvey_sources)
    fishery_provenance <- tools::md5sum(file.path(fishery_cache, "sampling-manifest.rds"))
  }
  generator <- if (phase == "controlled") pld_controlled_inputs else {
    function(id) pld_simsurvey_inputs(id, grid, file.path(out, sprintf("history-%02d", id)),
      fishery_cache = if (nzchar(fishery_cache)) file.path(fishery_cache, sprintf("history-%02d", id)))
  }
  input_phase <- if (phase == "controlled") phase else "sampling"
  manifest <- file.path(out, paste0(input_phase, "-manifest.rds"))
  generation <- list(n = n, grid = grid, revision = revision,
    generator = tools::md5sum(generating_sources), simsurvey = simsurvey_sources,
    R = as.character(getRversion()), data_table = as.character(packageVersion("data.table")))
  if (!is.null(fishery_provenance)) generation$fishery_cache <- fishery_provenance
  # Starting the fit-and-render workflow after a pure generation run can
  # change this driver while retaining the same simulation implementation.
  if (file.exists(manifest)) {
    existing <- readRDS(manifest)
    driver <- file.path(study, "simulate_pld.R")
    existing$generator <- existing$generator[names(existing$generator) != driver]
    generation$generator <- generation$generator[names(generation$generator) != driver]
    if (!identical(existing, generation)) {
      stop("Study inputs or generating code changed. Choose a new results directory.")
    }
  } else {
    saveRDS(generation, manifest)
  }
  if (phase == "sampling_fit") {
    source(file.path(study, "pld_validation_fit.R"))
    model_files <- sort(c(list.files("R", "[.]R$", full.names = TRUE), "NAMESPACE"))
    pld_record_sources(file.path(out, "fit-sampling-manifest.rds"),
      list(package = tools::md5sum(model_files),
        fitter = unname(tools::md5sum(file.path(study, "pld_validation_fit.R"))),
        R = as.character(getRversion()), RTMB = as.character(packageVersion("RTMB")),
        observation_input = tools::md5sum(file.path(out, "sampling-manifest.rds"))))
  }
  for (id in seq_len(n)) {
    file <- file.path(out, sprintf("%s-%02d.rds", input_phase, id))
    pld_checkpoint(file, function() generator(id))
    if (phase == "sampling_fit") {
      pld_fit_input(file, out, pld_checkpoint)
      # Render in a clean process so knitr never changes simulation RNG or
      # package versions for the next history.
      plot_dir <- Sys.getenv("TINYAM_PLD_REPORT_DIR")
      if (nzchar(plot_dir) && dir.exists(plot_dir)) {
        extra_dir <- if (identical(normalizePath(plot_dir), normalizePath(out))) "" else out
        report <- withr::with_envvar(c(TINYAM_PLD_REPORT_EXTRA = extra_dir),
          suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
            c("--vanilla", shQuote(file.path(study, "simulate_pld.R")), "report", n,
              grid, shQuote(plot_dir)), stdout = TRUE, stderr = TRUE)))
        writeLines(report)
      }
    }
  }
}

if (phase == "fit") {
  pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
  relevant <- list.files(out, "^(controlled|sampling)-[0-9]+[.]rds$", full.names = TRUE)
  if (!length(relevant)) stop("Generate data before running fits.")
  # Keep generation and fitting manifests separate so working on the fitted
  # likelihood does not require repeating identical monitoring simulations.
  model_files <- sort(c(list.files("R", "[.]R$", full.names = TRUE), "NAMESPACE"))
  provenance <- list(package = tools::md5sum(model_files),
    fitter = unname(tools::md5sum(file.path(study, "pld_validation_fit.R"))),
    R = as.character(getRversion()), RTMB = as.character(packageVersion("RTMB")))
  source(file.path(study, "pld_validation_fit.R"))
  for (input in relevant) {
    input_kind <- sub("-.*", "", basename(input))
    manifest <- file.path(out, paste0(input_kind, "-manifest.rds"))
    pld_record_sources(file.path(out, paste0("fit-", input_kind, "-manifest.rds")),
      c(provenance, list(observation_input = tools::md5sum(manifest))))
    pld_fit_input(input, out, pld_checkpoint)
  }
}

if (phase == "report") {
  extra <- Sys.getenv("TINYAM_PLD_REPORT_EXTRA")
  html <- rmarkdown::render(file.path(study, "pld_validation_report.Rmd"),
    params = list(results = out, extra_results = extra), output_dir = out, quiet = TRUE)
  cat("Report:", html, "\n")
}
