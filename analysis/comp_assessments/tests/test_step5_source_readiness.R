source("analysis/comp_assessments/R/run_assessment.R")

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

  database <- list(
    stocks = data.frame(stock_id = "fixture_stock"),
    assessments = data.frame(assessment_id = ids, stock_id = "fixture_stock",
      assessment_url = "https://example.invalid/assessment"),
    assumptions = data.frame(assessment_id = ids, component = "population",
      setting = "age_range", value = "unknown", source_reference = "fixture", notes = "fixture"),
    inputs = data.frame(assessment_id = character()),
    outputs = data.frame(assessment_id = character())
  )

  for (id in ids) {
    result <- run_assessment(id, database = database, fit = FALSE, tier = "preliminary")
    testthat::expect_identical(result$diagnostics$status, "source_readiness_failed")
    reason <- gsub("[[:space:]]+", " ", result$diagnostics$reason)
    testthat::expect_match(reason, expected[[id]], fixed = TRUE)
    testthat::expect_true(startsWith(result$background[[1L]], "###"))
    testthat::expect_true(any(grepl("https://example.invalid/assessment", result$background, fixed = TRUE)))
    testthat::expect_null(result$obs)
    testthat::expect_null(result$fit)
  }
})
