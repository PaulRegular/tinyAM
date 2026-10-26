root <- file.path("analysis", "comp_assessments")

# Batch controls. Dashboards are optional and are rendered only after saving fits.
parallel <- TRUE
workers <- 4L
dashboard <- FALSE
timeout <- 20 * 60 # Seconds per assessment; Inf disables the subprocess deadline.

pkgload::load_all(".", quiet = TRUE)
required <- c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "patchwork", "ggtext",
  if (is.finite(timeout)) "callr",
  if (dashboard) c("knitr", "rmarkdown", "flexdashboard"))
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) cli::cli_abort(c(
  "Install the batch's analysis packages before fitting: {paste(missing, collapse = ', ')}.",
  "i" = "These packages are used by the overview figures and optional subprocess/dashboard stages."
))
cli::cli_inform("Updating observation readiness.")
source("analysis/comp_assessments/scripts/database/003_observation_readiness.R")
source(file.path(root, "R", "run_assessment.R"))

cli::cli_inform("Fitting assessments; progress and checkpoints are saved under results/cache/<assessment_id>/.")
runs <- run_assessments(
  database = read_committed_database(),
  parallel = parallel,
  workers = workers,
  save_results = TRUE,
  cache = TRUE,
  dashboard = FALSE, # Render previews below, after saving the complete batch.
  timeout = timeout
)
cli::cli_inform("Saving combined assessment results.")
.assessment_save_rds(runs, file.path(root, "results", "cache", "assessment_runs.rds"))

# Exploratory dashboards do not change aggregate comparison eligibility.
preview_errors <- list()
if (dashboard) for (id in names(runs)) {
  result <- runs[[id]]
  if (!inherits(result$fit, "tam_fit")) next
  output <- file.path(root, "results", "cache", id, "dashboard.html")
  cli::cli_inform("Rendering dashboard: {id}")
  error <- tryCatch({
    view_assessment(result$fit, result, output_file = output, open_file = FALSE)
    runs[[id]]$dashboard_file <- output
    NULL
  }, error = conditionMessage)
  if (!is.null(error)) preview_errors[[id]] <- error
}
if (dashboard) .assessment_save_rds(runs, file.path(root, "results", "cache", "assessment_runs.rds"))
.assessment_save_rds(preview_errors, file.path(root, "results", "cache", "dashboard_errors.rds"))

cli::cli_inform("Building overview figures.")
source(file.path(root, "scripts", "translation", "stats_plots.R"))
source(file.path(root, "scripts", "translation", "trend_plots.R"))
cli::cli_inform("Assessment batch complete. Inspect results/fit_diagnostics.csv for failures and timeouts.")
