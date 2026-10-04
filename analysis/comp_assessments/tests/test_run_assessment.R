root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))

result <- run_assessment("dfo_cod_2j3kl_2025", fit = FALSE)
stopifnot(
  result$source$assessment$assessment_id == "dfo_cod_2j3kl_2025",
  result$source$database_source == "working_tree",
  is.data.frame(result$obs$catch),
  is.data.frame(result$obs$index),
  is.list(result$settings),
  is.data.frame(result$audit),
  length(result$background) > 0L,
  is.null(result$fit),
  result$diagnostics$status == "not_fitted",
  result$diagnostics$database_revision == result$source$commit
)

blocked <- run_assessment("afsc_cod_goa_2026", fit = FALSE)
stopifnot(
  blocked$source$assessment$assessment_id == "afsc_cod_goa_2026",
  blocked$diagnostics$status == "blocked_no_stock_spec",
  grepl("No stock translation recipe", blocked$diagnostics$reason)
)

cat("Single-assessment runner tests passed.\n")
