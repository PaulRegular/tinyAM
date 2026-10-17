source("analysis/comp_assessments/R/run_assessment.R")
x <- read_assessment("afsc_cod_ai_2024")
original <- x
for(fit in c(FALSE,TRUE)) {
  z <- .assessment_stock(x,do_fit=fit,silent=TRUE)
  stopifnot(inherits(z$error,"error"),
    grepl("matching final production",gsub("[[:space:]]+"," ",conditionMessage(z$error)),fixed=TRUE),
    z$env$fit_stage=="translation",is.null(z$env$obs),is.null(z$env$fit),
    is.null(z$env$fit_model),identical(x,original))
}
cat("AI cod source-run mismatch stops preparation/fitting without inventing final inputs or source starts.\n")
