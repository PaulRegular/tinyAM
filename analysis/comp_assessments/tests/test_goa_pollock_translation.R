pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/run_assessment.R")

assessment <- run_assessment("afsc_pollock_goa_2024", fit = FALSE)
comparison_outputs <- assessment$comparison$comparison_outputs
f <- comparison_outputs[comparison_outputs$measure == "fishing_mortality_at_age" &
  comparison_outputs$source_type == "derived_source_output", ]

stopifnot(
  assessment$diagnostics$status == "not_fitted",
  nrow(assessment$obs$catch) > 0,
  nrow(assessment$obs$index) > 0,
  nrow(f) == 371L,
  all(f$year %in% 1970:2023),
  all(f$age %in% 1:8),
  identical(assessment$comparison$comparison_definitions$F$status, "approximate"),
  any(grepl("partial reconstruction", assessment$background, fixed = TRUE))
)

cat("GOA pollock translation retains the approximate, partial F comparison without fitting.\n")
