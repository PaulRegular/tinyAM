source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
record <- read_assessment("dfo_petrale_bc_2024", database)
original <- record
for (fit in c(FALSE, TRUE)) {
  stock <- .assessment_stock(record, do_fit = fit, silent = TRUE)
  stopifnot(inherits(stock$error, "error"),
    grepl("original catch/survey age compositions and biological inputs",
      gsub("[[:space:]]+", " ", conditionMessage(stock$error)), fixed = TRUE),
    stock$env$fit_stage == "translation", is.null(stock$env$obs), is.null(stock$env$fit),
    is.null(stock$env$fit_model), is.null(stock$env$retry_start),
    is.null(stock$env$start_par), identical(record, original))
  result <- run_assessment("dfo_petrale_bc_2024", database, fit = fit,
    dashboard = TRUE, cache = FALSE)
  stopifnot(result$diagnostics$status == "translation_failed", is.null(result$fit),
    is.null(result$obs), is.null(result$ref), is.null(result$summary),
    is.null(result$differences), is.null(result$dashboard_file),
    is.na(result$diagnostics$objective), is.na(result$diagnostics$optimizer_code),
    is.na(result$diagnostics$sdreport_success), identical(record, original),
    nrow(result$audit) == nrow(record$assumptions),
    !"dfo_petrale_bc_2024" %in% database$assessments$assessment_id[database$assessments$is_current])
}
summary <- read.csv("analysis/comp_assessments/results/comparison_summary.csv")
series <- read.csv("analysis/comp_assessments/results/population_timeseries.csv")
stopifnot(!"dfo_petrale_bc_2024" %in% summary$assessment_id,
  !"dfo_petrale_bc_2024" %in% series$assessment_id)
cat("BC petrale remains source-blocked in preparation and fitting; no invented model or comparisons are produced.\n")
