root <- file.path("analysis", "comp_assessments")

source("analysis/comp_assessments/scripts/database/003_observation_readiness.R")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE,
  cache = TRUE,
  dashboard = TRUE
)
