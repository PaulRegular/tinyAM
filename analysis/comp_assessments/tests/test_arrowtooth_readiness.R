source("analysis/comp_assessments/R/run_assessment.R")
x <- read_assessment("afsc_arrowtooth_flounder_bsai_2022")
original <- x
for(fit in c(FALSE,TRUE)) {
  z <- .assessment_stock(x,do_fit=fit,silent=TRUE)
  stopifnot(inherits(z$error,"error"),
    grepl("numerical likelihood compositions",gsub("[[:space:]]+"," ",conditionMessage(z$error)),fixed=TRUE),
    z$env$fit_stage=="translation",is.null(z$env$obs),is.null(z$env$fit),
    is.null(z$env$fit_model),identical(x,original))
}
cat("Arrowtooth readiness trap preserves missing likelihood inputs and runs neither optimizer nor source-start construction.\n")
