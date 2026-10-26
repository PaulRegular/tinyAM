pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

# Exercise the actual top-level orchestration without fits, exports or rendering.
script <- parse("analysis/comp_assessments/scripts/translation/run_all.R")
for (dashboards in c(FALSE, TRUE)) {
  e <- new.env(parent = globalenv())
  e$source <- function(...) invisible(NULL)
  e$requireNamespace <- function(...) TRUE
  e$read_committed_database <- function() list()
  e$run_assessments <- function(...) {
    e$batch_args <- list(...)
    list(example = list(fit = structure(list(), class = "tam_fit")))
  }
  e$.assessment_save_rds <- function(...) invisible(NULL)
  e$previews <- 0L
  e$view_assessment <- function(...) e$previews <- e$previews + 1L
  for (expr in script) {
    if (is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
        identical(expr[[2L]], quote(dashboard))) expr[[3L]] <- dashboards
    eval(expr, e)
  }
  stopifnot(e$previews == as.integer(dashboards), !e$batch_args$dashboard,
    identical(e$batch_args$timeout, 1200), isTRUE(e$batch_args$cache))
}

missing_packages <- new.env(parent = globalenv())
missing_packages$requireNamespace <- function(...) FALSE
missing_packages$run_assessments <- function(...) stop("Fitting must not begin.")
missing_error <- tryCatch({
  for (expr in script) eval(expr, missing_packages)
  NULL
}, error = identity)
stopifnot(inherits(missing_error, "error"),
  grepl("Install the batch's analysis packages before fitting", conditionMessage(missing_error), fixed = TRUE))

database <- read_database()
id <- "ices_cod_north_sea_2025"
progress_text <- character()
withCallingHandlers(run_assessments(id, database, fit = FALSE, progress = TRUE),
  message = function(m) progress_text <<- c(progress_text, conditionMessage(m)))
stopifnot(any(grepl(paste(id, "started"), progress_text, fixed = TRUE)),
  any(grepl(paste(id, "not_fitted"), progress_text, fixed = TRUE)))
for (bad in list(0, -1, NA_real_, c(1, 2), "10")) {
  stopifnot(inherits(tryCatch(run_assessments(id, database, fit = FALSE,
    timeout = bad, progress = FALSE), error = identity), "error"))
}
stopifnot(inherits(tryCatch(run_assessments(id, database, fit = FALSE,
  progress = NA), error = identity), "error"))
empty_database <- database
empty_database$assessments$is_current <- FALSE
stopifnot(length(run_assessments(database = empty_database, progress = FALSE)) == 0L)

# A hard deadline must become a diagnostic and leave no live child blocking the
# next job. A tiny limit expires during R startup, without running a model.
timed <- run_assessments(id, database, fit = FALSE, timeout = .01, progress = FALSE)
stopifnot(timed[[id]]$diagnostics$status == "timed_out",
  is.null(timed[[id]]$fit), grepl("settings were not changed", timed[[id]]$diagnostics$reason),
  is.finite(timed[[id]]$diagnostics$elapsed_seconds))

ids <- c(id, "ices_herring_north_sea_2026")
for (parallel in c(FALSE, TRUE)) {
  ran <- run_assessments(ids, database, fit = FALSE, parallel = parallel,
    workers = 2L, timeout = 60, progress = FALSE)
  stopifnot(identical(names(ran), ids),
    all(vapply(ran, function(x) x$diagnostics$status == "not_fitted", logical(1))))
}

limited <- run_assessments(ids, database, fit = FALSE, parallel = TRUE,
  workers = 2L, timeout = .01, progress = FALSE)
stopifnot(identical(names(limited), ids),
  all(vapply(limited, function(x) x$diagnostics$status == "timed_out", logical(1))))

# Passing a real fitted RTMB object back from a supervised process must preserve
# the estimates, reporting tables, diagnostics and standalone tidying workflow.
direct <- run_assessment(id, database)
supervised <- run_assessments(id, database, timeout = 60, progress = FALSE)[[id]]
stopifnot(inherits(supervised$fit, "tam_fit"),
  isTRUE(all.equal(direct$fit$opt$objective, supervised$fit$opt$objective)),
  isTRUE(all.equal(direct$fit$parameter_values, supervised$fit$parameter_values)),
  identical(direct$diagnostics$status, supervised$diagnostics$status),
  isTRUE(all.equal(direct$summary, supervised$summary)))
tidied <- tinyAM::tidy_tam(model_list = list(supervised = supervised$fit))
stopifnot(tidied$process_residuals$status$status == "ok")

# Cache each stock before the next serial job is allowed to start.
original_runner <- run_assessment
original_root <- .assessment_root
test_root <- tempfile("assessment_checkpoint_")
dir.create(test_root)
.assessment_root <- test_root
run_assessment <- function(assessment_id, ...) {
  if (assessment_id == ids[[2L]]) {
    stopifnot(file.exists(file.path(test_root, "results", "cache", ids[[1L]], "result.rds")))
  }
  list(assessment_id = assessment_id, fit = NULL, differences = NULL, summary = NULL,
    diagnostics = data.frame(assessment_id = assessment_id, status = "not_fitted"))
}
tryCatch({
  cached <- run_assessments(ids, database, fit = FALSE, cache = TRUE, progress = FALSE)
  stopifnot(all(file.exists(file.path(test_root, "results", "cache", ids, "result.rds"))),
    all(vapply(cached, function(x) x$diagnostics$status == "not_fitted", logical(1))))
}, finally = {
  run_assessment <- original_runner
  .assessment_root <- original_root
  unlink(test_root, recursive = TRUE)
})
cat("Batch controls, subprocess deadlines and immediate checkpoints passed.\n")

checkpoint <- tempfile(fileext = ".rds")
.assessment_save_rds(list(version = 1), checkpoint)
interrupted <- tryCatch(.assessment_save_rds(stop("Interrupted serialization"), checkpoint), error = identity)
stopifnot(inherits(interrupted, "error"), identical(readRDS(checkpoint), list(version = 1)))
.assessment_save_rds(list(version = 2), checkpoint)
stopifnot(identical(readRDS(checkpoint), list(version = 2)))
unlink(checkpoint)
cat("Interrupted checkpoint preservation and replacement passed.\n")
