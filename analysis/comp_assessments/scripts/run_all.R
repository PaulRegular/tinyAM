# Run from the tinyAM repository root. Refresh the SAIO cache manually first.
root <- file.path("analysis", "comp_assessments")
translation_tier <- "curated"
timeout <- 20 * 60 # Seconds per assessment; Inf removes the hard deadline.

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))
database <- read_database()

# Configure parallelism here, as for fit_retro() and fit_hindcast().
# Use future::plan(future::sequential) when troubleshooting.
future::plan(future::multisession, workers = 1)

## Fit and save first ----
source(file.path(root, "scripts", "observation_readiness.R"))
runs <- run_assessments(database = database, tier = translation_tier,
  cache = TRUE, save_results = TRUE, timeout = timeout)
.assessment_save_rds(runs, file.path(.assessment_cache_root(translation_tier, database),
  "assessment_runs.rds"))
cli::cli_inform("Fitting complete. Inspect results/{translation_tier}/fit_diagnostics.csv for failures and timeouts.")

## Render separately; a render failure never loses a fit ----
dashboards <- render_assessment_dashboards(runs)

## Overview figures can also be rerun without fitting ----
source(file.path(root, "scripts", "stats_plots.R"))
source(file.path(root, "scripts", "trend_plots.R"))
future::plan(future::sequential)
