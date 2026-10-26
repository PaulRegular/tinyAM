source("analysis/comp_assessments/R/run_assessment.R")
x <- read_assessment("ices_anchovy_gulf_cadiz_2025")
original <- x
for (fit in c(FALSE, TRUE)) {
  result <- .assessment_stock(x, do_fit = fit, silent = TRUE)
  stopifnot(inherits(result$error, "error"),
    grepl("current empirical stock weights and native quarterly compositions", gsub("[[:space:]]+", " ", conditionMessage(result$error)), fixed = TRUE),
    result$env$fit_stage == "translation", is.null(result$env$obs), is.null(result$env$fit),
    is.null(result$env$start_par), is.null(result$env$fit_model), identical(x, original))
}
source("analysis/comp_assessments/R/database_to_tam_ref.R")
reference <- database_to_tam_ref(x$assessment$assessment_id, x$outputs,
  years = 1989:2025, ages = 0:3, terminal_year = 2025, age_plus_group = 3,
  assumptions = x$assumptions)
stopifnot(inherits(reference, "tam_ref"), all(!reference$dat$is_proj),
  is.null(reference$pop$N),
  all(is.na(reference$fixed_par$est)), is.null(reference$opt),
  nrow(reference$pop$ssb) == 37L, nrow(reference$pop$recruitment) == 37L,
  reference$pop$ssb$est[reference$pop$ssb$year == 2025] == 18730,
  reference$pop$recruitment$est[reference$pop$recruitment$year == 2025] == 11059600)
cat("Southern anchovy source gate avoids fitting, starting values and invented states; official-only reporting remains available.\n")
