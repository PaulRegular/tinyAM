root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))

result <- run_assessment("dfo_cod_2j3kl_2025", fit = FALSE)
stopifnot(
  result$source$assessment$assessment_id == "dfo_cod_2j3kl_2025",
  result$source$database_source == "working_tree",
  is.data.frame(result$obs$catch),
  is.data.frame(result$obs$index),
  is.data.frame(result$audit),
  length(result$background) > 0L,
  is.null(result$fit),
  result$diagnostics$status == "not_fitted",
  result$diagnostics$database_revision == result$source$commit
)

goa <- run_assessment("afsc_cod_goa_2026", fit = FALSE)
stopifnot(
  goa$source$assessment$assessment_id == "afsc_cod_goa_2026",
  goa$diagnostics$status == "not_fitted",
  nrow(goa$obs$catch) > 0,
  nrow(goa$obs$index) > 0
)

stopifnot(is.null(result$settings), is.null(goa$fit))

cat("Single-assessment runner tests passed.\n")

mock <- list(
  obs_pred = list(catch = data.frame(year = c(2000, 2000, 2001, 2001),
                                    age = c(1, 2, 1, 2), pred = c(10, 20, 30, 40))),
  pop = list(total_yield = data.frame(year = 2000:2001, est = 0),
             total_yield_pred = data.frame(year = 2000:2001, est = 0)),
  opt = list(objective = 123)
)
reporting <- list(
  weights = data.frame(year = c(2000, 2000, 2001, 2001),
                       age = c(1, 2, 1, 2), weight = c(0.1, 0.2, 0.3, 0.4)),
  totals = data.frame(year = 2000:2001, yield = c(6, 26))
)
reported <- .assessment_catch_reporting(mock, reporting)
stopifnot(identical(reported$pop$total_yield$est, c(6, 26)),
          identical(reported$pop$total_yield_pred$est, c(5, 25)),
          identical(reported$rep$total_yield, reported$pop$total_yield$est),
          identical(reported$rep$total_yield_pred, reported$pop$total_yield_pred$est),
          identical(reported$opt, mock$opt),
          all(mock$pop$total_yield$est == 0),
          identical(.assessment_catch_reporting(mock, NULL), mock))
reporting$weights <- reporting$weights[-1, ]
stopifnot(inherits(tryCatch(.assessment_catch_reporting(mock, reporting),
                          error = identity), "error"))
