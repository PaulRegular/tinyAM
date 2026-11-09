# Use run_observation_error.ps1 for serial execution and resource checks.
study <- file.path("analysis", "formula_expansion")
phase <- Sys.getenv("TINYAM_OE_PHASE")
grid <- as.integer(Sys.getenv("TINYAM_OE_GRID"))
gate <- Sys.getenv("TINYAM_OE_GATE")
population_phases <- c("populationpilot", "populationdiagnostic", "populationextension", "report")
z_phases <- c("zpilot", "zstudy", "zreport")
stopifnot(phase %in% c("pilot", "study", "samplingreport", population_phases, z_phases),
  grid %in% c(20L, 40L), nzchar(gate))
pilot <- phase == "pilot"
population_seed <- 21L
population_pilot <- phase %in% population_phases
results <- file.path(study, "results", if (phase %in% z_phases) {
  "observation_error_z_covariance"
} else if (population_pilot) {
  "observation_error_population_pilot"
} else "observation_error_seed21")
out <- file.path(results, paste0(if (phase %in% c("pilot", "zpilot")) "pilot" else "study", "-grid", grid))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

admit <- function() {
  repeat {
    signal <- tryCatch(suppressWarnings(readLines(gate, n = 1L)), error = function(e) "")
    stamp <- file.info(gate)$mtime
    if (identical(signal, "go") && is.finite(as.numeric(stamp)) &&
        difftime(Sys.time(), stamp, units = "secs") < 45) return(invisible(NULL))
    Sys.sleep(2)
  }
}

save_checkpoint <- function(value, path) {
  partial <- paste0(path, ".partial")
  saveRDS(value, partial, compress = FALSE)
  if (!file.rename(partial, path)) stop("Could not publish checkpoint ", path)
}

verify_manifest <- function(value, path) {
  if (!file.exists(path)) save_checkpoint(value, path)
  if (!identical(value, readRDS(path))) {
    stop("Incompatible checkpoint. Preserve ", dirname(path), " before starting a new run.")
  }
}

checkpoint <- function(path, run) {
  if (file.exists(path)) return(readRDS(path))
  admit()
  began <- Sys.time()
  warnings <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(run(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    NULL
  })
  elapsed <- as.numeric(difftime(Sys.time(), began, units = "secs"))
  attempt <- list(value = value, warnings = unique(warnings), error = error, elapsed = elapsed)
  save_checkpoint(attempt, path)
  cat(format(Sys.time()), basename(path), sprintf("%.2f s", elapsed),
    if (!is.null(error)) paste("FAILED", error) else "", "\n")
  flush.console()
  attempt
}

required <- function(attempt) {
  if (!is.null(attempt$error)) stop(attempt$error)
  attempt$value
}

admit()
options(warn = 1)
pilot_file <- file.path(study, "fit_observation_error_pilot.R")
extension_file <- file.path(study, "fit_observation_error_pilot_extension.R")
if (phase == "report") {
  source(pilot_file, encoding = "UTF-8")
  source(extension_file, encoding = "UTF-8")
  oe_pilot_extension_report(normalizePath(out, winslash = "/", mustWork = TRUE))
  quit(save = "no")
}

if (phase == "zreport") {
  result <- file.path(out, "summary.rds")
  if (!file.exists(result)) {
    stop("No saved Z covariance summary. Run -Phase ZStudy first.")
  }
  result <- normalizePath(result, winslash = "/")
  rmarkdown::render(file.path(study, "observation_error_z_report.Rmd"),
    params = list(result = result), output_dir = dirname(result), quiet = TRUE)
  quit(save = "no")
}

if (phase == "samplingreport") {
  result <- file.path(out, "summary.rds")
  if (!file.exists(result)) {
    stop("No saved sampling summary for this grid. Run -Phase Study first.")
  }
  result <- normalizePath(result, winslash = "/")
  rmarkdown::render(file.path(study, "observation_error_report.Rmd"),
    params = list(result = result), output_dir = dirname(result), quiet = TRUE)
  quit(save = "no")
}

pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
sim_file <- file.path(study, "observation_error_sim.R")
diag_file <- file.path(study, "observation_error_diagnostics.R")
source(sim_file, encoding = "UTF-8")
simsurvey <- normalizePath(Sys.getenv("TINYAM_OE_SIMSURVEY"), winslash = "/", mustWork = TRUE)
stopifnot(identical(basename(simsurvey), paste0("SimSurvey-", .observation_error_revision)))
pkgload::load_all(simsurvey, quiet = TRUE, helpers = FALSE)
Sys.setenv(TINYAM_SIMSURVEY_REVISION = .observation_error_revision)
.observation_error_check_revision()
data.table::setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

if (phase %in% c("zpilot", "zstudy")) {
  source(pilot_file, encoding = "UTF-8")
  source(file.path(study, "observation_error_z_sim.R"), encoding = "UTF-8")
  source(file.path(study, "observation_error_z_summary.R"), encoding = "UTF-8")
  z_pilot <- phase == "zpilot"
  z <- oe_z_simulate(out, if (z_pilot) "pilot" else "study", assessments = !z_pilot)
  z$rows <- subset(z$rows, bin == "modeled")
  z$sampling <- oe_z_summarize(z$rows)
  if (!z_pilot) {
    fit_files <- file.path(study, c("observation_error_fit.R", "observation_error_z_fit.R"))
    for (file in fit_files) source(file, encoding = "UTF-8")
    inputs <- z$assessment_inputs$file
    fit_sources <- c(fit_files, file.path(study, "observation_error_z_summary.R"),
      sort(list.files("R", "[.]R$", full.names = TRUE)))
    verify_manifest(list(sources = tools::md5sum(fit_sources),
      inputs = tools::md5sum(inputs[file.exists(inputs)])),
      file.path(out, "assessment-manifest.rds"))
    z$assessment <- oe_z_assessments(z, out)
  }
  z$rows <- z$overlap_maps <- NULL
  save_checkpoint(z, file.path(out, "summary.rds"))
  cat("Z covariance summary at", file.path(out, "summary.rds"), "\n")
  quit(save = "no")
}

if (phase == "populationextension") {
  source(pilot_file, encoding = "UTF-8")
  source(extension_file, encoding = "UTF-8")
  oe_pilot_extend()
  quit(save = "no")
}
if (phase == "populationpilot") {
  source(pilot_file, encoding = "UTF-8")
  oe_pilot_run()
  quit(save = "no")
}
if (phase == "populationdiagnostic") {
  source(pilot_file, encoding = "UTF-8")
  source(file.path(study, "observation_error_pilot_means.R"), encoding = "UTF-8")
  oe_pilot_mean_diagnostic()
  quit(save = "no")
}

field_seeds <- if (pilot) 17L else c(17L, 438L, 72L)
years <- if (pilot) 10L else c(5L, 10L, 16L)
n_split <- if (pilot) 3L else 30L
calibration <- seq_len(n_split)
validation <- n_split + calibration
oms <- c("coupled", "default", "iid_age")
efforts <- c("half", "standard", "double")
source_version <- list(schema = 2L, population_seed = population_seed, sim = tools::md5sum(sim_file),
  simsurvey_revision = .observation_error_revision,
  simsurvey_sources = tools::md5sum(list.files(file.path(simsurvey, "R"), "[.]R$", full.names = TRUE)),
  versions = vapply(c("tinyAM", "SimSurvey", "RTMB", "data.table", "mgcv"),
    function(p) as.character(utils::packageVersion(p)), character(1)),
  R = as.character(getRversion()))
seed_offset <- 600000000L
source_version$seed_offset <- seed_offset
manifest <- c(list(phase = if (pilot) "pilot" else "study", grid_n = grid, field_seeds = field_seeds,
  years = years, oms = oms, efforts = efforts, calibration = calibration, validation = validation),
  source_version)
manifest_path <- file.path(out, "manifest.rds")
verify_manifest(manifest, manifest_path)
pilot_out <- file.path(results, paste0("pilot-grid", grid))
if (!pilot) {
  previous <- readRDS(file.path(pilot_out, "summary.rds"))
  stopifnot(identical(previous$manifest[names(source_version)], source_version))
}

fixed_checkpoint <- function(path, run) {
  validated_pilot <- file.path(pilot_out, substring(path, nchar(out) + 2L))
  if (!pilot && !file.exists(path) && file.exists(validated_pilot)) {
    save_checkpoint(readRDS(validated_pilot), path)
  }
  required(checkpoint(path, run))
}

seed_for <- function(program, field_seed, year, replicate = 0L) {
  offsets <- c(removals = 100000000L, Gillnet = 200000000L, RV = 300000000L)
  as.integer(offsets[[program]] + field_seed * 100000L + year * 1000L + replicate +
    if (program == "removals") 0L else seed_offset)
}

draw_observations <- function(draw) {
  if (is.list(draw) && !is.data.frame(draw)) draw$rows else draw
}

cat("Starting", phase, "on", grid, "by", grid, "\n")
for (field_seed in field_seeds) {
  field_dir <- file.path(out, paste0("seed", field_seed))
  dir.create(field_dir, showWarnings = FALSE)
  cases <- expand.grid(om = oms, year = years, stringsAsFactors = FALSE)
  complete <- file.path(field_dir, paste0(cases$om, "-", cases$year), "done.rds")
  if (all(file.exists(complete))) next
  pops <- fixed_checkpoint(file.path(field_dir, "populations.rds"),
    function() oe_population(field_seed, grid, population_seed))
  stopifnot(setequal(names(pops), oms))
  for (om in oms) {
    path <- file.path(field_dir, paste0(om, "-overlap.rds"))
    required(checkpoint(path, function() oe_overlap(pops[[om]])))
    for (year in years) {
      case_dir <- file.path(field_dir, paste0(om, "-", year))
      dir.create(case_dir, showWarnings = FALSE)
      done <- file.path(case_dir, "done.rds")
      if (file.exists(done)) next
      admit()
      snapshot <- oe_snapshot(pops[[om]], year)
      truth <- fixed_checkpoint(file.path(case_dir, "removals.rds"),
        function() oe_gillnet_truth(snapshot, seed_for("removals", field_seed, year)))
      logs <- list()
      for (replicate in c(calibration, validation)) {
        for (component in c("Gillnet", "RV")) {
          seed <- seed_for(component, field_seed, year, replicate)
          path <- file.path(case_dir, sprintf("%s-%03d.rds", component, replicate))
          attempt <- checkpoint(path, function() {
            draw <- if (component == "Gillnet") oe_gillnet_draw(truth, seed) else oe_rv_draw(snapshot, seed)
            rows <- data.table::as.data.table(draw_observations(draw))
            stopifnot(nrow(rows) == 102L, !anyDuplicated(rows[, .(age, effort, bin)]),
              all(rows$year == year),
              setequal(rows$effort, efforts), all(rows$component == component))
            if (component == "RV") {
              sets <- unique(rows[, .(effort, n_primary)])
              stopifnot(all(sets$n_primary == c(half = 80, standard = 160, double = 320)[sets$effort]))
            }
            rows[, ':='(om = om, population_id = paste0("seed", field_seed),
              replicate = replicate)]
            rows
          })
          logs[[length(logs) + 1L]] <- data.frame(
            om = om, population_id = paste0("seed", field_seed), year = year,
            component = component, replicate = replicate, seed = seed,
            file = path, elapsed = attempt$elapsed, success = is.null(attempt$error),
            warnings = paste(attempt$warnings, collapse = " | "),
            error = if (is.null(attempt$error)) "" else attempt$error)
          if (pilot && !is.null(attempt$error)) stop(attempt$error)
        }
      }
      save_checkpoint(do.call(rbind, logs), done)
      rm(snapshot, truth, logs, attempt)
      invisible(gc())
    }
  }
  rm(pops)
  invisible(gc())
}

admit()
attempts <- data.table::rbindlist(lapply(list.files(out, "^done[.]rds$",
  recursive = TRUE, full.names = TRUE), readRDS))
rows <- data.table::rbindlist(lapply(attempts$file[attempts$success],
  function(path) readRDS(path)$value), fill = TRUE)
stopifnot(nrow(attempts) == length(field_seeds) * length(years) * length(oms) * 2L * 2L * n_split,
  nrow(rows) == 102L * sum(attempts$success))
overlap <- lapply(field_seeds, function(seed) {
  stats::setNames(lapply(oms, function(om) required(readRDS(file.path(out,
    paste0("seed", seed), paste0(om, "-overlap.rds"))))), oms)
})
names(overlap) <- paste0("seed", field_seeds)

source(diag_file, encoding = "UTF-8")
oe_diagnostic_checks()
diagnostics_key <- unname(tools::md5sum(diag_file))
diagnostics <- required(checkpoint(file.path(out, paste0("diagnostics-", diagnostics_key, ".rds")),
  function() oe_diagnostics(rows, calibration, validation, bootstrap = if (pilot) 30L else 400L)))
summary <- list(manifest = manifest, draw_attempts = attempts, overlap = overlap,
  diagnostics = diagnostics, diagnostics_key = diagnostics_key)
save_checkpoint(summary, file.path(out, "summary.rds"))
cat("Study summary at", file.path(out, "summary.rds"), "\n")
