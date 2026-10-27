## Select and fit a stock ----
root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "afsc_atka_mackerel_bsai_2024"
database <- read_database()
source <- read_assessment(assessment_id, database)
do_fit <- TRUE
silent <- FALSE

stock <- .assessment_stock(source, do_fit, silent)
if (inherits(stock$error, "error")) stop(stock$error)
assessment <- .assessment_context(source, as.list(stock$env))
list2env(as.list(stock$env), envir = environment())
# Clear optional metadata when changing stocks in the same workspace.
list2env(assessment$comparison, envir = environment())
if (!is.null(fit)) fit <- .assessment_catch_reporting(fit, catch_reporting)
fit

## Preview the current fit ----
if (!is.null(fit)) preview <- view_assessment(fit, assessment)

