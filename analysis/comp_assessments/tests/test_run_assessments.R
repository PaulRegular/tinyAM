pkgload::load_all(".", quiet = TRUE)
root <- file.path("analysis", "comp_assessments")
source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
assessment_ids <- c("ices_cod_north_sea_2025", "dfo_cod_2j3kl_2025")
current_ids <- database$assessments$assessment_id[
  !is.na(database$assessments$is_current) & database$assessments$is_current
]
current_batch <- run_assessments(database = database, fit = FALSE)
stopifnot(
  setequal(names(current_batch), current_ids),
  "ices_norway_pout_north_sea_2026_benchmark" %in% names(current_batch),
  current_batch$ices_norway_pout_north_sea_2026_benchmark$diagnostics$status ==
    "not_fitted"
)
comparison_example <- .assessment_comparison_summary(
  data.frame(metric = "ssb", year = 2020:2021, age = NA_integer_,
             source = c(100, 90), tinyAM = c(95, 92),
             percent_difference = c(-5, 100 * 2 / 90)),
  "example_assessment"
)
stopifnot(comparison_example$assessment_id == "example_assessment")
constant_comparison <- .assessment_comparison_summary(
  data.frame(metric = "M", year = 2020:2021, age = 1L,
             source = 0.2, tinyAM = 0.2,
             percent_difference = c(0, 0)),
  "constant_assessment"
)
stopifnot(is.na(constant_comparison$trend_correlation))

serial <- run_assessments(assessment_ids, database = database, fit = FALSE)
stopifnot(
  identical(names(serial), assessment_ids),
  all(vapply(serial, function(x) x$diagnostics$status == "not_fitted", logical(1))),
  identical(attr(serial, "diagnostics")$assessment_id, assessment_ids),
  "assessment_id" %in% names(attr(serial, "comparison_summary"))
)

previous_plan <- future::plan("list")[[1L]]
previous_workers <- future::nbrOfWorkers()
shared_outputs <- file.path(root, "results", c(
  "fit_diagnostics.csv", "comparison_summary.csv"
))
shared_outputs_existed <- file.exists(shared_outputs)
cache_existed <- dir.exists(file.path(root, "results", "cache"))
parallel <- run_assessments(
  assessment_ids, database = database, parallel = TRUE, workers = 2,
  fit = FALSE
)
stopifnot(
  identical(names(parallel), assessment_ids),
  identical(attr(future::plan("list")[[1L]], "call"),
            attr(previous_plan, "call")),
  future::nbrOfWorkers() == previous_workers,
  identical(file.exists(shared_outputs), shared_outputs_existed),
  identical(dir.exists(file.path(root, "results", "cache")), cache_existed),
  all(vapply(parallel, function(x) x$diagnostics$status == "not_fitted", logical(1))),
  identical(attr(parallel, "diagnostics")$assessment_id, assessment_ids)
)

bad_database <- database
bad_database$inputs <- bad_database$inputs[
  !(bad_database$inputs$assessment_id == "dfo_cod_2j3kl_2025" &
      bad_database$inputs$type == "maturity"), , drop = FALSE]
failed <- run_assessments(
  c("dfo_cod_2j3kl_2025", "ices_cod_north_sea_2025"),
  database = bad_database, parallel = TRUE, workers = 2, fit = FALSE
)
goa <- run_assessments("afsc_cod_goa_2026", database = database,
                       fit = FALSE)
stopifnot(
  failed$dfo_cod_2j3kl_2025$diagnostics$status == "translation_failed",
  !is.null(failed$ices_cod_north_sea_2025$obs),
  goa$afsc_cod_goa_2026$diagnostics$status == "not_fitted",
  !is.null(goa$afsc_cod_goa_2026$obs)
)

cat("Multi-assessment runner tests passed.\n")
