root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "dfo_cod_4t4vn_2019"

## Used for interactive tweaks to the translation script
# source <- read_assessment(assessment_id, database = read_database())

x <- run_assessment(assessment_id, database = read_database(), silent = FALSE)
x$fit

vis_tam(model_list = list(tinyAM = x$fit, Accepted = x$ref), background = x$background)

## For interactive troubleshooting dashboard
# fits <- list(tinyAM = x$fit, Accepted = x$ref); interval <- 0.95
