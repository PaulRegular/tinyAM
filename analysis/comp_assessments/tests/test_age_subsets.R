pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

id <- "ices_herring_western_baltic_2026"
source_data <- read_assessment(id, read_database())
before <- source_data
full <- .test_stock(source_data)
assessment <- .assessment_context(source_data, full)

testthat::test_that("fixed M permits subsets and still requires complete selected ages", {
  fixed <- database_to_tam_M(id, source_data$inputs, years = 1991:1993, ages = 0:7)
  testthat::expect_equal(nrow(fixed$surface), 24L)
  testthat::expect_true(all(fixed$surface$age %in% 0:7))
  source_m <- subset(source_data$inputs, type == "M" & measure == "natural_mortality_at_age")
  testthat::expect_equal(fixed$surface$M_assumption,
    as.numeric(source_m$value[match(fixed$surface$age, source_m$age)]))
  for (ages in list(1:7, 2:5)) {
    selected <- database_to_tam_M(id, source_data$inputs, years = 1991:1993, ages = ages)
    testthat::expect_setequal(selected$surface$age, ages)
  }
  missing <- source_data$inputs[!(source_data$inputs$type == "M" & source_data$inputs$age == 3), ]
  testthat::expect_error(database_to_tam_M(id, missing, years = 1991:1993, ages = 0:7),
    "every requested age")
  duplicate <- rbind(source_data$inputs, source_m[source_m$age == 3, ])
  testthat::expect_error(database_to_tam_M(id, duplicate, years = 1991:1993, ages = 0:7),
    "every requested age")
  testthat::expect_error(database_to_tam_M(id, source_data$inputs, years = 1991:1993, ages = 0:9),
    "every requested age")

  annual <- source_m[rep(1L, 4), ]
  annual$age <- NA_real_
  annual$year <- 1990:1993
  annual$year_basis <- "calendar_year"
  annual$value <- c(.1, .2, .3, .4)
  selected <- database_to_tam_M(id, annual, years = 1991:1992, ages = 0:7)
  testthat::expect_equal(selected$surface$M_assumption, rep(c(.2, .3), 8))
  testthat::expect_error(database_to_tam_M(id, annual[annual$year != 1992, ],
    years = 1991:1992, ages = 0:7), "every requested year")
})

testthat::test_that("observations select biology and retain older catch in the new plus age", {
  obs <- database_to_tam_obs(id, source_data$inputs, years = full$years, ages = 0:7,
    surveys = c("HERAS", "GERAS", "N20", "IBTS/BITSQ1"), assumptions = source_data$assumptions)
  testthat::expect_true(tinyAM::check_obs(obs))
  for (component in names(obs)) testthat::expect_true(all(obs[[component]]$age %in% 0:7))
  testthat::expect_equal(nrow(obs$weight), 35L * 8L)
  old_plus <- subset(full$obs$catch, age >= 7)
  expected <- aggregate(obs ~ year, old_plus, sum)
  new_plus <- subset(obs$catch, age == 7)
  testthat::expect_equal(new_plus$obs, expected$obs[match(new_plus$year, expected$year)])
  testthat::expect_equal(sum(obs$catch$obs), sum(full$obs$catch$obs))
  testthat::expect_equal(obs$index$obs, full$obs$index$obs)
  older <- subset(source_data$inputs, type == "index" & measure == "numbers_at_age" & survey == "HERAS")[1:2, ]
  older$year <- 1991
  older$age <- 7:8
  older$value <- c(10, 20)
  with_older <- database_to_tam_obs(id, rbind(source_data$inputs, older),
    years = full$years, ages = 0:7, surveys = c("HERAS", "GERAS", "N20", "IBTS/BITSQ1"),
    assumptions = source_data$assumptions)
  plus_index <- subset(with_older$index, year == 1991 & age == 7 & survey == "HERAS")
  testthat::expect_equal(plus_index$obs,
    sum(older$value * vapply(older$unit, .translation_index_multiplier, numeric(1))))
  testthat::expect_error(database_to_tam_obs(id, source_data$inputs,
    years = full$years, ages = c(0, 2, 3)), "consecutive")
  testthat::expect_identical(source_data, before)
})

# Evaluate the stock recipe exactly as when its ages line is edited to 0:7.
e <- new.env(parent = globalenv())
e$source <- source_data
e$do_fit <- FALSE
e$silent <- TRUE
script <- file.path(.assessment_root, "scripts/translation/stocks", paste0(id, ".R"))
for (expr in parse(script)) {
  if (is.call(expr) && identical(expr[[1]], quote(`<-`)) && identical(expr[[2]], quote(ages))) {
    expr[[3]] <- quote(0:7)
  }
  eval(expr, e)
}
fit <- do.call(tinyAM::fit_tam, c(list(data = e$obs, years = e$years, ages = e$ages,
                                      start_par = e$start_par, silent = TRUE), full$settings))

testthat::test_that("a reduced-age recipe fits and previews against the original assessment", {
  testthat::expect_identical(e$ages, 0:7)
  testthat::expect_equal(e$age_plus_group, 7)
  testthat::expect_true(isTRUE(fit$is_converged))
  testthat::expect_true(isTRUE(fit$sdrep$pdHess))
  fit_before <- fit
  context_before <- assessment
  result <- .assessment_comparison(fit, assessment)
  testthat::expect_identical(result$ref$dat$ages, 0:7)
  for (metric in c("N", "F", "M")) {
    testthat::expect_equal(nrow(result$ref$pop[[metric]]), 35L * 8L)
    testthat::expect_true(all(result$ref$pop[[metric]]$age %in% 0:7))
  }
  n_source <- subset(source_data$outputs, type == "population" & measure == "numbers_at_age" & age >= 7)
  expected_n <- aggregate(value ~ year, n_source, sum)
  plus <- subset(result$ref$pop$N, age == 7)
  testthat::expect_equal(plus$est, expected_n$value[match(plus$year, expected_n$year)])
  testthat::expect_true(all(is.na(plus$se)))
  testthat::expect_true(all(is.na(plus$lwr)))
  testthat::expect_true(all(is.na(plus$upr)))
  for (metric in c("F", "M")) {
    wanted_measure <- if (metric == "F") "fishing_mortality_at_age" else "natural_mortality_at_age"
    mortality <- subset(source_data$outputs, type == "mortality" & measure == wanted_measure & age >= 7)
    n <- n_source$value[match(paste(mortality$year, mortality$age), paste(n_source$year, n_source$age))]
    expected <- tapply(seq_len(nrow(mortality)), mortality$year,
      function(i) weighted.mean(mortality$value[i], n[i]))
    plus <- subset(result$ref$pop[[metric]], age == 7)
    testthat::expect_equal(plus$est, as.numeric(expected[as.character(plus$year)]))
  }
  testthat::expect_true(all(is.finite(result$ref$pop$abundance$est)))
  testthat::expect_true(all(is.finite(result$ref$pop$ssb$est)))
  testthat::expect_true(any(attr(result$ref, "native_pop")$N$age == 8))
  testthat::expect_identical(fit, fit_before)
  testthat::expect_identical(assessment, context_before)
  testthat::expect_identical(source_data, before)

  html <- tempfile(fileext = ".html")
  preview <- view_assessment(fit, assessment, output_file = html, open_file = FALSE)
  testthat::expect_true(file.exists(preview$dashboard_file))
  testthat::expect_gt(file.info(html)$size, 1000)
  unlink(html)
})

testthat::test_that("a reported group outside the fit is unavailable rather than partly compared", {
  grouped <- assessment
  grouped$comparison$comparison_age_groups <- list(N = list(`7+` = 7:8))
  result <- .assessment_comparison(fit, grouped)
  rows <- subset(result$differences, metric == "N")
  testthat::expect_true(all(rows$comparison_status == "unavailable"))
  testthat::expect_true(all(is.na(rows$percent_difference)))
  testthat::expect_true(all(grepl("outside the fitted age range", rows$reason)))
  duplicate <- grouped
  duplicate$comparison$comparison_age_groups <- list(N = list(`6+` = 6:7, `7+` = 7:8))
  testthat::expect_error(.assessment_comparison(fit, duplicate), "distinct fitted ages")
})

testthat::test_that("a direct reduced-age fit can reuse the original full-age observations", {
  direct <- do.call(tinyAM::fit_tam, c(list(data = full$obs, years = full$years,
    ages = 0:7, start_par = e$start_par, silent = TRUE), full$settings))
  testthat::expect_true(isTRUE(direct$is_converged))
  comparison <- .assessment_comparison(direct, assessment)
  testthat::expect_identical(comparison$ref$dat$ages, 0:7)
  testthat::expect_true(all(comparison$ref$pop$N$age %in% 0:7))
  testthat::expect_true(all(is.finite(comparison$ref$pop$abundance$est)))
  testthat::expect_identical(direct$refit_args$data, full$obs)
  testthat::expect_identical(source_data, before)
})

cat("Age-subset translation and dashboard tests passed.\n")
