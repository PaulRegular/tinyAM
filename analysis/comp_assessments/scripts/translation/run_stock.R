root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "dfo_herring_4tvn_spring_2024"

## Used for interactive tweaks to the translation script
# source <- read_assessment(assessment_id, database = read_database())

x <- run_assessment(assessment_id, database = read_database(), silent = FALSE)
list2env(x[c("source", "obs", "fit", "background")], envir = environment())
do_fit <- TRUE
silent <- FALSE
if (!is.null(fit)) {
  years <- fit$dat$years
  ages <- fit$dat$ages
}
fit

vis_tam(model_list = list(tinyAM = fit, Accepted = x$ref), background = background)

## For interactive troubleshooting dashboard
# fits <- list(tinyAM = x$fit, Accepted = x$ref); interval <- 0.95
