root <- file.path("analysis", "comp_assessments")

pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "dfo_herring_4tvn_spring_2024"

database <- read_database()
source <- read_assessment(assessment_id, database)
do_fit <- TRUE
silent <- FALSE
base::source(file.path(root, "scripts", "translation", "stocks",
                       paste0(assessment_id, ".R")), local = TRUE)
fit
