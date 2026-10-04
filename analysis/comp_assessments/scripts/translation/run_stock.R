root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "dfo_cod_2j3kl_2025"
x <- run_assessment(assessment_id, database = read_database())

source_data <- x$source
obs <- x$obs
settings <- x$settings
fit <- x$fit
ref <- x$ref
audit <- x$audit
diagnostics <- x$diagnostics
comparison <- x$summary
