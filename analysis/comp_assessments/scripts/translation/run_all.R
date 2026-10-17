root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/scripts/database/003_observation_readiness.R")
source(file.path(root, "R", "run_assessment.R"))

runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE,
  cache = TRUE,
  dashboard = FALSE
)
saveRDS(runs, file.path(root, "results", "cache", "assessment_runs.rds"))

# Exploratory dashboards do not change aggregate comparison eligibility.
preview_errors <- list()
for (id in names(runs)) {
  result <- runs[[id]]
  if (!inherits(result$fit, "tam_fit")) next
  output <- file.path(root, "results", "cache", id, "dashboard.html")
  error <- tryCatch({
    view_assessment(result$fit, result, output_file = output, open_file = FALSE)
    runs[[id]]$dashboard_file <- output
    NULL
  }, error = conditionMessage)
  if (!is.null(error)) preview_errors[[id]] <- error
}
saveRDS(runs, file.path(root, "results", "cache", "assessment_runs.rds"))
saveRDS(preview_errors, file.path(root, "results", "cache", "dashboard_errors.rds"))

source(file.path(root, "scripts", "translation", "stats_plots.R"))
source(file.path(root, "scripts", "translation", "trend_plots.R"))
