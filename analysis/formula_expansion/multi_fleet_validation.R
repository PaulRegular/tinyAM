# Run from the repository root: Rscript analysis/formula_expansion/multi_fleet_validation.R core
# The artificial catch split below exercises workflows; it is not an assessment
# trial or a recovery experiment. Scientific trials use separately recorded truth.
pkgload::load_all(".", quiet = TRUE)

fleet_workflow_inputs <- function() {
  d <- tinyAM::cod_obs
  d$catch <- rbind(transform(d$catch, fleet = "Trawl", obs = .6 * obs),
    transform(d$catch, fleet = "Gillnet", obs = .4 * obs))
  d
}

run_fleet_core_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/core") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character(); stage <- "fit"
  results <- withCallingHandlers(tryCatch({
    fit <- tinyAM::fit_tam(fleet_workflow_inputs(), years = 1983:2000, ages = 2:8,
      N_settings = list(process = "off"),
      rec_settings = list(process = "off", mu_form = ~ factor(year)),
      F_settings = list(process = "off", mu_form = ~ fleet + factor(year) + factor(age)),
      M_settings = list(process = "off", mu_supplied = ~ I(.3)), silent = TRUE)
    diagnostics <- tinyAM::check_tam(fit)
    stage <- "projection"
    projected <- update(fit, start_par = fit,
      proj_settings = list(n_proj = 2, n_mean = 2, catch = 1000), silent = TRUE)
    stopifnot(max(abs(projected$rep$total_yield_pred[projected$dat$is_proj] - 1000)) < 1e-6)
    stage <- "simulation"
    set.seed(6515)
    simulations <- tinyAM::sim_tam(fit, n = 2, par_uncertainty = "none")
    stage <- "hindcast"
    hindcast <- tinyAM::fit_hindcast(fit, folds = 1, progress = FALSE)
    stage <- "dashboard"
    tinyAM::vis_tam(fit, open_file = FALSE, output_file = file.path(output, "dashboard.html"),
      render_args = list(quiet = TRUE))
    list(fit = fit, diagnostics = diagnostics, projected = projected,
      simulations = simulations, hindcast = hindcast, error = NULL)
  }, error = function(e) list(error = conditionMessage(e), stage = stage)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings
  results$session <- utils::sessionInfo()
  saveRDS(results, file.path(output, "workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Fleet workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics)
  print(results$hindcast)
  invisible(results)
}

run_fleet_active_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/active") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character(); results <- list()
  withCallingHandlers(tryCatch({
    results$fit <- tinyAM::fit_tam(fleet_workflow_inputs(), years = 1983:2000, ages = 2:8,
      F_settings = list(process = "cor_rw", sd_form = ~fleet, cor_form = ~fleet,
        age_breaks = list(Gillnet = c(2,4,6,9), Trawl = NULL)),
      M_settings = list(process = "off", mu_supplied = ~I(.3)), silent = TRUE)
    results$diagnostics <- tinyAM::check_tam(results$fit)
    tinyAM::vis_tam(results$fit, open_file = FALSE,
      output_file = file.path(output,"dashboard.html"),render_args = list(quiet = TRUE))
  }, error = function(e) { results$error <<- conditionMessage(e) }),
    warning = function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results,file.path(output,"workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Active-F rendering check failed: {results$error}")
  print(results$diagnostics)
  invisible(results)
}

run_fleet_scoped_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/scoped") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  d <- fleet_workflow_inputs()
  d$index <- rbind(transform(d$index, survey = "RV"),
    transform(d$index, survey = "Sentinel", obs = .8 * obs))
  settings <- list(data = d, years = 1983:2000, ages = 2:8,
    N_settings = list(process = "off"),
    rec_settings = list(process = "off", mu_form = ~ factor(year)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3)), silent = TRUE)
  warnings <- character(); results <- list(); stage <- "fit"
  withCallingHandlers(tryCatch({
    results$base <- do.call(tinyAM::fit_tam, c(settings, list(
      F_settings = list(process = "off", mu_form = ~ fleet + factor(year) + factor(age)),
      index_settings = list(q_form = ~ 0 + survey:q_block))))
    results$scoped <- do.call(tinyAM::fit_tam, c(settings, list(
      F_settings = list(process = "off",
        mu_form = ~ 0 + factor(age) + factor(year) + for_fleet("Trawl", ~ 1)),
      index_settings = list(q_form = ~ 0 + for_survey("RV", ~ q_block) +
        for_survey("Sentinel", ~ q_block)))))
    results$diagnostics <- tinyAM::check_tam(results$base, results$scoped)
    results$AIC <- stats::AIC(results$base, results$scoped)
    results$cAIC <- tinyAM::cAIC(results$base, results$scoped, details = TRUE)
    stage <- "hindcast"
    results$hindcast <- tinyAM::fit_hindcast(results$scoped, folds = 1, progress = FALSE)
    stage <- "osa"
    before <- results$scoped$parameter_values
    results$osa <- stats::residuals(results$scoped, type = "osa", subset = 1:3, trace = FALSE)
    stopifnot(identical(before, results$scoped$parameter_values))
    stage <- "dashboard"
    tinyAM::vis_tam(Base = results$base, Scoped = results$scoped, open_file = FALSE,
      output_file = file.path(output, "dashboard.html"), render_args = list(quiet = TRUE))
  }, error = function(e) { results$error <<- conditionMessage(e); results$stage <<- stage }),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results, file.path(output, "workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Scoped workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics); print(results$cAIC); print(results$hindcast)
  invisible(results)
}

if (sys.nframe() == 0L) {
  mode <- commandArgs(TRUE)
  if (!length(mode)) mode <- "core"
  run <- list(core = run_fleet_core_workflow, active = run_fleet_active_workflow,
    scoped = run_fleet_scoped_workflow)
  if (!mode[1L] %in% names(run)) cli::cli_abort("Available validation modes: {paste(names(run), collapse = ', ')}.")
  run[[mode[1L]]]()
}
