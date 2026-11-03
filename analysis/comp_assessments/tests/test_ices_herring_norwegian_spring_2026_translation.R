source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "ices_herring_norwegian_spring_2026"
database <- read_database(tier = "preliminary")
result <- run_assessment(assessment_id, database = database, fit = FALSE)
obs <- result$obs

stopifnot(
  identical(result$diagnostics$status, "not_fitted"),
  is.null(result$fit),
  tinyAM::check_obs(obs),
  identical(sort(unique(as.integer(obs$catch$year))), 1988:2025),
  identical(sort(unique(as.integer(obs$catch$age))), 2:12),
  nrow(obs$catch) == 38L * 11L,
  nrow(obs$weight) == 38L * 11L,
  nrow(obs$maturity) == 38L * 11L,
  all(obs$pop_covariates$M_assumption[obs$weight$age == 2] == 0.9),
  all(obs$pop_covariates$M_assumption[obs$weight$age > 2] == 0.15),
  setequal(unique(obs$index$survey), c(
    "NASF_1988_2008", "NASF_2015_onward", "IESNS_Barents",
    "IESNS_Norwegian_Sea", "BESS"
  )),
  !"RFID" %in% unique(obs$index$survey)
)

cat("Norwegian spring-spawning herring 2026 translation preparation passed.\n")
