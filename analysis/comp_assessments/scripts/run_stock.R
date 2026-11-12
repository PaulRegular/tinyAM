## Select and fit a stock ----
root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

translation_tier <- "curated"
assessment_id <- "noaa_sablefish_pacific_coast_2025"
# Use the same SAIO inputs as the curated sablefish comparison.
database <- read_database(revision = "312bc7210bca3c88e23de5fd882fb4ab0eb03d99")
source <- read_assessment(assessment_id, database)
# Set FALSE to load the inputs for an edited fit_tam() call.
do_fit <- TRUE
silent <- FALSE

stock <- .assessment_stock(source, do_fit, silent, tier = translation_tier)
if (inherits(stock$error, "error")) stop(stock$error)
assessment <- .assessment_context(source, as.list(stock$env))
list2env(as.list(stock$env), envir = environment())
# Clear optional metadata when changing stocks in the same workspace.
list2env(assessment$comparison, envir = environment())
if (!is.null(fit)) fit <- .assessment_catch_reporting(fit, catch_reporting)
fit

## Preview the current fit ----
if (!is.null(fit)) preview <- view_assessment(fit, assessment)

## Try an alternative interactively ----
# trial <- update(fit, N_settings = list(init = "free"))
# check_tam(trial)
# view_assessment(trial, assessment)

