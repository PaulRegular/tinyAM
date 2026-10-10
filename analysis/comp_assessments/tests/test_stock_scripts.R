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
cat('Top-level stock scripts and fitting-disabled isolation passed.\n')
