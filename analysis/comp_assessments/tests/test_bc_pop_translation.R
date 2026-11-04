source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("dfo_pop_bc_2023")
original <- record
for (fit in c(FALSE, TRUE)) {
  x <- .assessment_stock(record, do_fit = fit, silent = TRUE)
  stopifnot(inherits(x$error, "error"),
    grepl("historical official population outputs",
      gsub("[[:space:]]+", " ", conditionMessage(x$error)), fixed = TRUE),
    x$env$fit_stage == "translation", is.null(x$env$fit), is.null(x$env$obs),
    is.null(x$env$fit_model), is.null(x$env$retry_start), identical(record, original))
  result <- run_assessment("dfo_pop_bc_2023", fit = fit, cache = FALSE)
  stopifnot(result$diagnostics$status == "translation_failed",
    is.null(result$fit), is.null(result$ref), is.null(result$summary),
    is.null(result$dashboard_file), identical(record, original))
}
tables <- .assessment_result_tables(list(result), save_results = FALSE)
stopifnot(nrow(attr(tables, "comparison_summary")) == 0L,
  nrow(attr(tables, "population_timeseries")) == 0L)
cat("BC POP is blocked before fitting/comparison; source records remain intact and no stale comparison rows remain.\n")
