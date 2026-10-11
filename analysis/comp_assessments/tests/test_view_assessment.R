pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

assessment <- run_assessment("ices_cod_north_sea_2025")
fit <- assessment$fit
stopifnot(inherits(fit, "tam_fit"), !is.null(assessment$ref))

testthat::test_that("previews use the supplied fit, refresh comparisons and do not export results", {
  before <- fit
  source_before <- assessment$source
  aggregate_files <- file.path(.assessment_root, "results",
    c("fit_diagnostics.csv", "comparison_summary.csv"))
  hashes <- tools::md5sum(aggregate_files)
  captured <- NULL
  testthat::local_mocked_bindings(
    fit_tam = function(...) stop("A preview must never fit a model."),
    vis_tam = function(...) captured <<- list(...), .package = "tinyAM")
  result <- withVisible(view_assessment(fit, assessment, open_file = FALSE))
  testthat::expect_false(result$visible)
  testthat::expect_false(captured$open_file)
  testthat::expect_identical(names(captured$model_list), c("Accepted", "tinyAM"))
  testthat::expect_identical(captured$background, assessment$background)
  testthat::expect_identical(captured$output_file, result$value$dashboard_file)
  testthat::expect_identical(dirname(captured$output_file), normalizePath(tempdir(), winslash = "/"))
  testthat::expect_false(file.exists(captured$output_file))
  order_summary <- function(x) {
    x <- x[order(x$metric), ]
    rownames(x) <- NULL
    x
  }
  testthat::expect_equal(order_summary(result$value$summary), order_summary(assessment$summary))
  view_assessment(fit, assessment)
  testthat::expect_true(captured$open_file)
  changed <- fit
  changed$pop$ssb$est <- fit$pop$ssb$est * 1.5
  result2 <- view_assessment(changed, assessment, output_file = "experiment.html", open_file = FALSE)
  testthat::expect_identical(result2$dashboard_file, "experiment.html")
  testthat::expect_identical(captured$model_list$tinyAM$pop$ssb$est, changed$pop$ssb$est)
  testthat::expect_false(isTRUE(all.equal(result2$summary, result$value$summary)))
  testthat::expect_equal(result2$ref$comparisons$tinyAM[result2$ref$comparisons$metric == "ssb"],
    changed$pop$ssb$est * result2$ref$comparison_scales[["ssb"]])
  testthat::expect_identical(fit, before)
  testthat::expect_identical(assessment$source, source_before)
  testthat::expect_identical(tools::md5sum(aggregate_files), hashes)
})

testthat::test_that("comparison construction reads the current fit grid and biological inputs", {
  captured <- NULL
  original <- database_to_tam_ref
  testthat::local_mocked_bindings(database_to_tam_ref = function(...) {
    captured <<- list(...)
    original(...)
  }, .env = environment(.assessment_comparison))
  changed <- fit
  changed$refit_args$data$maturity$obs <- changed$refit_args$data$maturity$obs * .8
  .assessment_comparison(changed, assessment)
  testthat::expect_identical(captured$obs, changed$refit_args$data)
  testthat::expect_identical(captured$years, changed$dat$years)
  testthat::expect_identical(captured$ages, changed$dat$ages)
  testthat::expect_identical(captured$comparison_scales, assessment$comparison$comparison_scales)
  testthat::expect_identical(captured$assumptions, assessment$source$assumptions)
  incompatible <- assessment
  incompatible$comparison$comparison_age_groups <- list(N = list(`99+` = 99:100))
  testthat::expect_error(.assessment_comparison(fit, incompatible), "distinct fitted ages")
  incompatible$comparison$comparison_age_groups <- list(list(`1+` = 1:2))
  testthat::expect_error(.assessment_comparison(fit, incompatible), "named list of metrics")
  incompatible$comparison$comparison_age_groups <- NULL
  incompatible$comparison$age_plus_group <- max(fit$dat$ages) + 1
  testthat::expect_error(.assessment_comparison(fit, incompatible), "fitted plus age")
})

testthat::test_that("failed numerical checks do not block usable exploratory previews", {
  captured <- NULL
  testthat::local_mocked_bindings(vis_tam = function(...) captured <<- list(...), .package = "tinyAM")
  failed <- fit
  failed$opt$convergence <- 1L
  testthat::expect_warning(result <- view_assessment(failed, assessment, open_file = FALSE),
    "not passed numerical convergence")
  testthat::expect_false(result$diagnostics$is_converged)
  testthat::expect_identical(captured$model_list$tinyAM$opt$convergence, 1L)
  empty <- fit
  empty$obs_pred$catch$pred[] <- NA_real_
  empty$obs_pred$index$pred[] <- NA_real_
  testthat::expect_error(view_assessment(empty, assessment), "no usable reporting predictions")
  testthat::expect_error(view_assessment(NULL, assessment), "fitted tinyAM")
  testthat::expect_error(view_assessment(fit, list()), "assessment context")
  testthat::expect_error(view_assessment(fit, assessment, open_file = NA), "open_file")
  testthat::expect_error(view_assessment(fit, assessment, output_file = NA_character_), "output_file")
})

testthat::test_that("projections rebuild reference grids and summaries use historical years", {
  projected <- suppressWarnings(update(fit, proj_settings = list(n_proj = 2), silent = TRUE))
  result <- .assessment_comparison(projected, assessment)
  testthat::expect_identical(result$ref$dat$years, projected$dat$years)
  testthat::expect_identical(result$ref$dat$is_proj, projected$dat$is_proj)
  testthat::expect_true(all(stats::na.omit(result$summary$terminal_year) <=
    max(projected$dat$years[!projected$dat$is_proj])))
})

testthat::test_that("unavailable metrics retain their status and explanation in summaries", {
  reference <- .assessment_comparison(fit, assessment)$ref
  row <- reference$comparisons[1L, ]
  row$metric <- "unavailable_metric"
  row$year <- NA_integer_
  row$source <- row$tinyAM <- NA_real_
  row$comparison_status <- "unavailable"
  row$reason <- "No numerical accepted output is available."
  reference$comparisons <- rbind(reference$comparisons, row)
  result <- testthat::with_mocked_bindings(.assessment_comparison(fit, assessment),
    database_to_tam_ref = function(...) reference,
    .env = environment(.assessment_comparison))
  summary <- result$summary[result$summary$metric == "unavailable_metric", ]
  testthat::expect_equal(nrow(summary), 1L)
  testthat::expect_identical(summary$comparison_status, "unavailable")
  testthat::expect_identical(summary$reason, row$reason)
  testthat::expect_equal(summary$n, 0)
})

testthat::test_that("missing uncertainty still permits a real temporary comparison dashboard", {
  # Recreate reporting exactly as fit_tam does when sdreport fails.
  no_se <- fit
  no_se$sdrep <- NULL
  parameters <- tinyAM::tidy_par(no_se)
  no_se$fixed_par <- parameters$fixed
  no_se$random_par <- parameters$random
  no_se$obs_pred <- tinyAM::tidy_obs_pred(no_se)
  no_se$pop <- tinyAM::tidy_pop(no_se)
  no_se$rec <- tinyAM::tidy_recruitment(no_se)
  before <- no_se
  testthat::expect_warning(result <- view_assessment(no_se, assessment, open_file = FALSE),
    "not passed numerical convergence")
  testthat::expect_true(file.exists(result$dashboard_file))
  testthat::expect_true(all(is.na(result$fit$pop$ssb$se)))
  testthat::expect_true(all(is.finite(result$fit$pop$ssb$est)))
  testthat::expect_identical(no_se, before)
  connection <- file(result$dashboard_file, "rb")
  html <- paste(readLines(connection, warn = FALSE), collapse = "\n")
  close(connection)
  testthat::expect_match(html, "Numerical checks", fixed = TRUE)
  testthat::expect_match(html, "Accepted", fixed = TRUE)
  unlink(result$dashboard_file)
})

cat("Interactive assessment preview tests passed.\n")
