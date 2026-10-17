source("analysis/comp_assessments/R/run_assessment.R")
x <- read_assessment("dfo_herring_sog_2025")
original <- x
for(fit in c(FALSE,TRUE)) {
  z <- .assessment_stock(x,do_fit=fit,silent=TRUE)
  stopifnot(inherits(z$error,"error"),
    grepl("aggregate spawning biomass",gsub("[[:space:]]+"," ",conditionMessage(z$error)),fixed=TRUE),
    z$env$fit_stage=="translation",is.null(z$env$obs),is.null(z$env$fit),
    is.null(z$env$fit_model),identical(x,original))
}
cat("SoG herring fails clearly before fitting or starting-value construction; source data unchanged.\n")
