root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE
)
