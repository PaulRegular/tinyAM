source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
database <- read_database(tier = "preliminary")

testthat::test_that("the first expansion recipes preserve source values and never fit during preparation", {
  ids <- c("ices_ling_faroe_5b_2026", "ices_atlantic_wolffish_iceland_2026")
  testthat::local_mocked_bindings(fit_tam = function(...) {
    stop("Preparation must not fit or construct a source-start retry")
  }, .package = "tinyAM")
  for (id in ids) {
    record <- read_assessment(id, database)
    original <- record
    stock <- .test_stock(record)
    testthat::expect_null(stock$fit)
    testthat::expect_identical(record, original)
    testthat::expect_identical(stock$settings$N_settings$init, "free")
    testthat::expect_identical(stock$settings$F_settings$process, "rw")
    testthat::expect_true(all(abs(stock$obs$pop_covariates$M_assumption - 0.15) < 1e-12))
    testthat::expect_true(all(stock$obs$index$unit == "native survey index"))
    observed <- stock$obs$index
    source_index <- record$inputs[record$inputs$type == "index" &
      record$inputs$measure %in% c("numbers_at_age", "relative_abundance_index") &
      record$inputs$year %in% stock$years, ]
    key <- function(x) paste(x$survey, x$year, x$age)
    testthat::expect_equal(observed$obs, source_index$value[match(key(observed), key(source_index))])
    if (id == "ices_ling_faroe_5b_2026") {
      testthat::expect_equal(stock$years, 1996:2025)
      testthat::expect_equal(stock$ages, 3:12)
      spring <- grepl("spring", observed$survey)
      testthat::expect_equal(observed$samp_time, ifelse(spring, 2 / 12, 7.5 / 12))
    } else {
      testthat::expect_equal(stock$years, 1998:2025)
      testthat::expect_equal(stock$ages, 4:16)
      testthat::expect_false(any(stock$obs$catch$year %in% 1981:1997))
      raw_weight <- record$inputs[record$inputs$type == "weight" &
        record$inputs$year %in% stock$years, ]
      index <- match(paste(stock$obs$weight$year, stock$obs$weight$age),
        paste(raw_weight$year, raw_weight$age))
      testthat::expect_equal(stock$obs$weight$obs, raw_weight$value[index] / 1000)
    }
  }
})

testthat::test_that("whiting retains framework state and q sharing on the production grid", {
  record <- read_assessment("ices_whg_27_7b_ce_k_2026", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 1999:2025)
  testthat::expect_equal(stock$ages, 0:7)
  testthat::expect_equal(stock$settings$F_settings$age_breaks, c(0:6, 8))
  testthat::expect_identical(as.character(stock$dat$F_settings$age_blocks["6"]),
    as.character(stock$dat$F_settings$age_blocks["7"]))
  index <- stock$obs$index
  testthat::expect_true(all(index$samp_time == 0.875))
  testthat::expect_identical(unique(as.character(index$q_group[index$age == 1])),
    unique(as.character(index$q_group[index$age == 2])))
  original <- record$inputs[record$inputs$type == "index" &
    record$inputs$measure == "numbers_at_age", ]
  testthat::expect_equal(index$obs,
    original$value[match(paste(index$year, index$age), paste(original$year, original$age))] / 1000)
  testthat::expect_true(all(index$unit == "thousand fish per km2"))
})

testthat::test_that("redfish reconstructs only sampled ages and keeps biological roles separate", {
  record <- read_assessment("noaa_acadian_redfish_gom_gb_2025", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 1969:2024)
  testthat::expect_identical(stock$settings$N_settings$process, "off")
  testthat::expect_identical(stock$settings$F_settings$process, "off")
  testthat::expect_identical(stock$settings$rec_settings$process, "ar1")
  testthat::expect_true(all(stock$obs$maturity$spawn_time == 0.4))
  testthat::expect_true(all(is.na(stock$obs$catch$obs[stock$obs$catch$year == 2000])))
  proportions <- record$inputs[record$inputs$type == "catch" &
    record$inputs$measure == "proportion_at_age" & record$inputs$year == 2024, ]
  weights <- record$inputs[record$inputs$type == "catch_weight", ]
  weights <- weights$value[match(proportions$age, weights$age)]
  total <- record$inputs$value[record$inputs$type == "catch" &
    record$inputs$measure == "total_biomass" & record$inputs$year == 2024]
  expected <- total * proportions$value / sum(proportions$value * weights)
  observed <- stock$obs$catch[stock$obs$catch$year == 2024, ]
  testthat::expect_equal(observed$obs[match(proportions$age, observed$age)], expected)
  for (survey in unique(stock$obs$index$survey)) {
    index <- stock$obs$index[stock$obs$index$survey == survey, ]
    plateau <- if (grepl("fall", survey)) 4 else 7
    testthat::expect_length(unique(index$q_group[index$age >= plateau]), 1L)
  }
})

testthat::test_that("Irish herring records a source-date blocker before any fitting", {
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Unexpected fit"),
    .package = "tinyAM")
  result <- run_assessment("ices_her_27_nirs_2026", database,
    fit = FALSE, tier = "preliminary")
  testthat::expect_identical(result$diagnostics$status, "source_readiness_failed")
  testthat::expect_match(result$diagnostics$reason, "reference dates")
  testthat::expect_null(result$fit)
})
