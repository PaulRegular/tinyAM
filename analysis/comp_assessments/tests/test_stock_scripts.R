source('analysis/comp_assessments/R/run_assessment.R')
source('analysis/comp_assessments/tests/helper_stock.R')
database <- read_database()
paths <- list.files('analysis/comp_assessments/scripts/translation/stocks', full.names=TRUE)
for (path in paths) {
  id <- sub('[.]R$', '', basename(path))
  record <- read_assessment(id, database)
  original <- record
  stock <- .test_stock(record)
  stopifnot(identical(record, original), is.null(stock$fit),
            is.null(stock$settings$settings), is.call(stock$fit_call),
            identical(stock$fit_call[[1]], quote(tinyAM::fit_tam)))
  # Script variables must stay local when the runner loads another stock.
  result <- run_assessment(id, database, fit=FALSE)
  stopifnot(result$diagnostics$status == 'not_fitted',
            is.null(result$settings), is.null(result$fit),
            isTRUE(all.equal(result$obs, stock$obs)))
}
stopifnot(!exists('warm_fit', inherits=FALSE))

# Run the interactive setup with fitting disabled, in a separate workspace.
interactive_env <- new.env(parent = globalenv())
for (statement in parse('analysis/comp_assessments/scripts/translation/run_stock.R')) {
  if (is.call(statement) && identical(statement[[1]], quote(`<-`)) &&
      identical(statement[[2]], quote(do_fit))) statement[[3]] <- FALSE
  eval(statement, interactive_env)
}
stopifnot(is.null(interactive_env$fit), is.list(interactive_env$obs),
          length(interactive_env$background) > 0L,
          identical(interactive_env$source$assessment$assessment_id, 'dfo_herring_4tvn_spring_2024'))
cat('Top-level stock scripts and fitting-disabled isolation passed.\n')

# Loading another stock must clear absent comparison options, including in the workspace.
interactive_env$comparison_age_groups <- list(N = list(stale = 999))
for (statement in parse('analysis/comp_assessments/scripts/translation/run_stock.R')) {
  if (is.call(statement) && identical(statement[[1]], quote(`<-`))) {
    if (identical(statement[[2]], quote(do_fit))) statement[[3]] <- FALSE
    if (identical(statement[[2]], quote(assessment_id))) statement[[3]] <- 'afsc_cod_goa_2026'
  }
  eval(statement, interactive_env)
}
stopifnot(interactive_env$assessment$assessment_id == 'afsc_cod_goa_2026',
          is.null(interactive_env$assessment$comparison$comparison_age_groups),
          is.null(interactive_env$comparison_age_groups),
          interactive_env$assessment$comparison$age_plus_group == 10,
          is.null(interactive_env$fit))
