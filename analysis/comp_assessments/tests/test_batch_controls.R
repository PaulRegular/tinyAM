pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

# The script exposes one caller-owned future plan and a separate render stage.
script <- paste(readLines("analysis/comp_assessments/scripts/run_all.R"), collapse = "\n")
stopifnot(grepl("future::plan(future::multisession, workers = 10)", script, fixed = TRUE),
  grepl("timeout <- 20 * 60", script, fixed = TRUE),
  grepl("render_assessment_dashboards(runs)", script, fixed = TRUE),
  !any(c("parallel", "workers", "dashboard") %in% names(formals(run_assessments))),
  !"dashboard" %in% names(formals(run_assessment)),
  identical(formals(run_assessments)$timeout, quote(Inf)))

database <- read_database()
ids <- c("ices_cod_north_sea_2025", "ices_herring_north_sea_2026")
old_plan <- future::plan()
future::plan(future::sequential)
progress_text <- character()
withCallingHandlers(run_assessments(ids[1], database, fit = FALSE, progress = TRUE),
  message = function(m) progress_text <<- c(progress_text, conditionMessage(m)))
stopifnot(any(grepl(paste(ids[1], "started"), progress_text, fixed = TRUE)),
  any(grepl(paste(ids[1], "not_fitted"), progress_text, fixed = TRUE)))
for (bad in list(0, -1, NA_real_, c(1, 2), "10")) {
  stopifnot(inherits(tryCatch(run_assessments(ids[1], database, fit = FALSE,
    timeout = bad, progress = FALSE), error = identity), "error"))
}
stopifnot(inherits(tryCatch(run_assessments(ids[1], database, fit = FALSE,
  progress = NA), error = identity), "error"))
empty_database <- database
empty_database$assessments <- empty_database$assessments[FALSE, , drop = FALSE]
stopifnot(length(run_assessments(database = empty_database, progress = FALSE)) == 0L)

# Both plans retain all attempts and the same prepared observations.
sequential <- run_assessments(ids, database, fit = FALSE, progress = FALSE)
future::plan(future::multisession, workers = 2)
parallel <- run_assessments(ids, database, fit = FALSE, progress = FALSE)
if (!isTRUE(all.equal(lapply(sequential, `[[`, "obs"), lapply(parallel, `[[`, "obs"), tolerance = 0))) {
  print(lapply(parallel, function(x) x$diagnostics[, c("status", "reason")]))
}
stopifnot(inherits(future::plan(), "multisession"),
  identical(names(sequential), names(parallel)),
  isTRUE(all.equal(lapply(sequential, `[[`, "obs"), lapply(parallel, `[[`, "obs"), tolerance = 0)),
  all(vapply(parallel, function(x) x$diagnostics$status == "not_fitted", logical(1))))
# A tiny deadline expires in startup without running or retuning a model.
limited <- run_assessments(ids, database, fit = FALSE, timeout = .01, progress = FALSE)
stopifnot(identical(names(limited), ids),
  all(vapply(limited, function(x) x$diagnostics$status == "timed_out" &&
    grepl("settings were not changed", x$diagnostics$reason), logical(1))))
future::plan(future::sequential)

# A fitted RTMB result must survive the supervised-process return path.
direct <- run_assessment(ids[1], database)
supervised <- run_assessments(ids[1], database, timeout = 120, progress = FALSE)[[ids[1]]]
stopifnot(inherits(supervised$fit, "tam_fit"),
  isTRUE(all.equal(direct$fit$opt$objective, supervised$fit$opt$objective)),
  isTRUE(all.equal(direct$fit$parameter_values, supervised$fit$parameter_values)),
  identical(direct$diagnostics$status, supervised$diagnostics$status),
  isTRUE(all.equal(direct$summary, supervised$summary)))
tidied <- tinyAM::tidy_tam(model_list = list(supervised = supervised$fit))
stopifnot(tidied$process_residuals$status$status == "ok")

# Every completed attempt is checkpointed before the next sequential job starts.
original_runner <- run_assessment
original_root <- .assessment_root
test_root <- tempfile("assessment_checkpoint_")
dir.create(test_root)
.assessment_root <- test_root
run_assessment <- function(assessment_id, ...) {
  if (assessment_id == ids[[2L]]) {
    stopifnot(file.exists(file.path(.assessment_cache(ids[[1L]], "curated", database), "result.rds")))
  }
  list(assessment_id = assessment_id, fit = NULL, differences = NULL, summary = NULL,
    diagnostics = data.frame(assessment_id = assessment_id, status = "not_fitted",
      database_revision = database$commit, database_key = .assessment_cache_key(database)))
}
tryCatch({
  cached <- run_assessments(ids, database, fit = FALSE, cache = TRUE, progress = FALSE)
  stopifnot(all(file.exists(file.path(.assessment_cache(ids, "curated", database), "result.rds"))))
}, finally = {
  run_assessment <- original_runner
  .assessment_root <- original_root
  unlink(test_root, recursive = TRUE)
})

checkpoint <- tempfile(fileext = ".rds")
.assessment_save_rds(list(version = 1), checkpoint)
interrupted <- tryCatch(.assessment_save_rds(stop("Interrupted serialization"), checkpoint), error = identity)
stopifnot(inherits(interrupted, "error"), identical(readRDS(checkpoint), list(version = 1)))
.assessment_save_rds(list(version = 2), checkpoint)
stopifnot(identical(readRDS(checkpoint), list(version = 2)))
unlink(checkpoint)

# Rendering never refits, rereads the source, changes eligibility or loses failures.
local({
  old_root <- .assessment_root; old_view <- view_assessment
  old_reader <- read_database; old_runner <- run_assessment
  output <- tempfile("assessment_render_"); dir.create(output)
  .assessment_root <<- output
  on.exit({.assessment_root <<- old_root; view_assessment <<- old_view
    read_database <<- old_reader; run_assessment <<- old_runner
    unlink(output, recursive = TRUE)}, add = TRUE)
  read_database <<- function(...) stop("Rendering must not reread source data")
  run_assessment <<- function(...) stop("Rendering must not refit")
  view_assessment <<- function(fit, assessment, output_file, open_file) {
    if (assessment$assessment_id == "failed") stop("Example render error")
    writeLines("<html>Rendered</html>", output_file)
  }
  result <- function(id, convergence) list(assessment_id = id,
    translation_tier = "curated", source = database,
    fit = structure(list(is_converged = convergence), class = "tam_fit"),
    diagnostics = data.frame(database_key = database$database_key),
    summary = data.frame(marker = "unchanged"))
  failed_source <- database; failed_source$translation_tier <- "curated"
  empty <- .failed_assessment_run("absent", failed_source, simpleError("No fitted model"))
  runs <- list(absent = empty, good = result("good", TRUE),
    uncertain = result("uncertain", FALSE), failed = result("failed", TRUE))
  original <- serialize(runs, NULL)
  rendered <- render_assessment_dashboards(runs)
  stopifnot(identical(rendered$status, c("skipped_no_fit", "rendered", "troubleshooting", "render_failed")),
    rendered$error[1] == "No fitted model", rendered$error[4] == "Example render error", identical(serialize(runs, NULL), original),
    file.exists(file.path(output, "results", "curated", "index.html")),
    file.exists(file.path(.assessment_cache_root("curated", database), "dashboard_manifest.csv")))
  all_failed <- render_assessment_dashboards(list(absent = empty))
  stopifnot(identical(all_failed$status, "skipped_no_fit"),
    file.exists(file.path(.assessment_cache_root("curated", failed_source), "dashboard_manifest.csv")))
})
future::plan(old_plan)
cat("Caller-owned future plans, deadlines, checkpoints and separate dashboard status checks passed.\n")
