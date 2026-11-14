source("analysis/comp_assessments/R/translation_registry.R")

testthat::test_that("Step 5 preliminary recipes stop before fitting when source data are not age-ready", {
  ids <- c(
    "noaa_goa_rex_sole_2026",
    "iphc_pacific_halibut_2025",
    "noaa_chilipepper_2025"
  )
  expected <- c(
    noaa_goa_rex_sole_2026 = "final-run age observations",
    iphc_pacific_halibut_2025 = "age-specific observations",
    noaa_chilipepper_2025 = "not age-specific"
  )
  registry <- assessment_translations("preliminary")
  entries <- registry[match(ids, registry$assessment_id), , drop = FALSE]

  testthat::expect_equal(entries$assessment_id, ids)
  testthat::expect_true(all(is.na(entries$source_commit) | !nzchar(entries$source_commit)))
  testthat::expect_true(all(file.exists(entries$script)))

  for (id in ids) {
    env <- new.env(parent = globalenv())
    env$source <- list(assessment = list())
    env$print_sources <- function(assessment) character()
    condition <- tryCatch({
      sys.source(.assessment_script(id, "preliminary"), envir = env)
      NULL
    }, error = identity)

    testthat::expect_s3_class(condition, "error")
    testthat::expect_match(conditionMessage(condition), expected[[id]])
    testthat::expect_identical(env$fit_stage, "source_readiness")
    testthat::expect_true(startsWith(env$background[[1]], "###"))
    testthat::expect_false(exists("obs", envir = env, inherits = FALSE))
    testthat::expect_false(exists("fit", envir = env, inherits = FALSE))
  }
})
