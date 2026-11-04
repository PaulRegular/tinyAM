source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
database <- read_database(tier = "preliminary")

testthat::test_that("the expansion register keeps all twenty distinct source targets visible", {
  batch <- read.csv("analysis/comp_assessments/expansion_batch_01.csv", stringsAsFactors = FALSE)
  testthat::expect_equal(nrow(batch), 20L)
  testthat::expect_false(anyDuplicated(batch$stock_id) > 0L)
  testthat::expect_false(anyDuplicated(batch$assessment_id) > 0L)
  testthat::expect_true(all(batch$assessment_id %in% database$assessments$assessment_id))
  # The batch register records its original source review; do not relabel it
  # when current default reads advance to SAIO main.
  testthat::expect_true(all(batch$data_revision == .saio_snapshot_revision("preliminary")))
  testthat::expect_true(all(nzchar(batch$limitation)))
  translations <- assessment_translations("preliminary")
  selected <- translations[match(batch$assessment_id, translations$assessment_id), ]
  testthat::expect_false(anyNA(selected$assessment_id))
  testthat::expect_true(all(file.exists(selected$script)))
})

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

testthat::test_that("remaining NOAA cases record input blockers without fitting", {
  cases <- list(
    noaa_pacific_coast_dover_sole_2021 = "age samples do not represent discards",
    noaa_west_coast_yelloweye_2025 = "age composition",
    noaa_canary_rockfish_2023 = "Conditional ages",
    noaa_black_rockfish_oregon_2023 = "age-specific catches and survey observations"
  )
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Unexpected fit"),
    .package = "tinyAM")
  for (id in names(cases)) {
    result <- run_assessment(id, database, fit = FALSE, tier = "preliminary")
    testthat::expect_identical(result$diagnostics$status, "source_readiness_failed")
    testthat::expect_match(result$diagnostics$reason, cases[[id]])
    testthat::expect_null(result$fit)
    testthat::expect_null(result$obs)
    testthat::expect_true(any(grepl("###", result$background)))
  }
})

testthat::test_that("other source-blocked cases have explicit preliminary recipes", {
  cases <- c(
    dfo_greenland_halibut_2_3klmno_2026_mp = "management procedure",
    dfo_yellowtail_flounder_4t_2025 = "numeric age inputs",
    dfo_western_component_pollock_2025_hcr_update = "HCR update",
    dfo_yellowmouth_2021 = "age-composition arrays",
    ices_cod_west_greenland_inshore_2026 = "named 2026 SAM object",
    ices_cod_west_greenland_offshore_2026 = "named 2026 SAM object",
    ices_seabass_8ab_2026 = "run-matched observations"
  )
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Unexpected fit"),
    .package = "tinyAM")
  for (id in names(cases)) {
    result <- run_assessment(id, database, fit = FALSE, tier = "preliminary")
    testthat::expect_identical(result$diagnostics$status, "source_readiness_failed")
    testthat::expect_match(result$diagnostics$reason, cases[[id]])
    testthat::expect_null(result$fit)
    testthat::expect_null(result$obs)
  }
})

testthat::test_that("Celtic herring retains model-year labels, source conflicts and exclusions", {
  record <- read_assessment("ices_her_27_irls_2026", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 2002:2024)
  testthat::expect_identical(attr(stock$obs, "translation")$year_origin, "1 April")
  testthat::expect_false(any(stock$obs$index$year %in% c(2004, 2017, 2024)))
  testthat::expect_true(all(stock$obs$index$samp_time == 0.5))
  testthat::expect_equal(stock$obs$index$obs[stock$obs$index$year == 2019 & stock$obs$index$age == 2], 98229)
  raw <- record$inputs[record$inputs$type == "catch" & record$inputs$measure == "numbers_at_age", ]
  key <- function(x) paste(x$year, x$age)
  testthat::expect_equal(stock$obs$catch$obs, raw$value[match(key(stock$obs$catch), key(raw))])
  testthat::expect_true(all(stock$obs$catch$unit == "thousand fish"))
  testthat::expect_length(unique(stock$obs$index$q_group[stock$obs$index$age %in% 3:5]), 1L)
})

testthat::test_that("GOA ocean perch uses source biological roles and retains older biology", {
  record <- read_assessment("noaa_goa_pop_2023", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 1990:2023)
  testthat::expect_equal(stock$ages, 2:25)
  testthat::expect_true(all(stock$obs$index$samp_time == 0))
  testthat::expect_equal(max(stock$obs$weight$age), 29)
  testthat::expect_equal(max(stock$obs$maturity$age), 29)
  testthat::expect_identical(stock$settings$M_settings$process, "off")
  testthat::expect_identical(deparse(stock$settings$M_settings$mu_form), "~1")
  testthat::expect_null(stock$settings$M_settings$mu_supplied)
  testthat::expect_null(stock$dat$M_settings$mu_supplied)
  testthat::expect_false(any(is.finite(stock$obs$catch$obs[stock$obs$catch$year == 1991])))
  w <- record$inputs[record$inputs$type == "weight" & is.na(record$inputs$survey), ]
  testthat::expect_equal(stock$obs$weight$obs, w$value[match(stock$obs$weight$age, w$age)] / 1000)
  p <- record$outputs[record$outputs$measure == "maturity_at_age" & record$outputs$year == 1961, ]
  testthat::expect_equal(stock$obs$maturity$obs, p$value[match(stock$obs$maturity$age, p$age)] * 0.5)
})

testthat::test_that("winter flounder uses one count representation and estimates source-shaped M", {
  record <- read_assessment("dfo_sgsl_winter_flounder_2017", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 1973:2016)
  testthat::expect_equal(stock$ages, 2:12)
  testthat::expect_null(stock$dat$M_settings$mu_supplied)
  testthat::expect_equal(stock$dat$rec_settings$phi, plogis(-4))
  testthat::expect_length(unique(stock$obs$pop_covariates$M_group), 2L)
  testthat::expect_length(unique(stock$obs$pop_covariates$M_period), 9L)
  testthat::expect_true(all(stock$obs$index$samp_time == 0.71))
  testthat::expect_equal(exp(stock$dat$log_sd_catch_supplied),
    ifelse(stock$dat$obs$catch$year == 1984, 2, 0.1))
  testthat::expect_equal(ncol(stock$dat$sd_catch_modmat), 0L)
  for (type in c("catch", "index")) {
    original <- record$inputs[record$inputs$type == type & record$inputs$measure == "numbers_at_age", ]
    current <- stock$obs[[type]]
    testthat::expect_equal(current$obs, original$value[match(paste(current$year, current$age),
      paste(original$year, original$age))])
  }
})

testthat::test_that("sablefish reconstruction conserves source biomass and never fills partial removals", {
  record <- read_assessment("noaa_sablefish_pacific_coast_2025", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 2006:2023)
  testthat::expect_equal(stock$ages, 0:30)
  testthat::expect_null(stock$dat$M_settings$mu_supplied)
  testthat::expect_true(all(is.na(stock$obs$catch$obs[stock$obs$catch$year == 2014])))
  testthat::expect_true(all(is.finite(stock$obs$catch$obs[stock$obs$catch$year != 2014])))
  testthat::expect_equal(stock$catch_reconstruction$source_tonnes,
    stock$catch_reconstruction$reconstructed_tonnes)
  testthat::expect_true(all(stock$obs$index$unit == "native survey index"))
  testthat::expect_identical(unique(stock$obs$index$survey), "WCGBTS")
  testthat::expect_true(all(stock$obs$index$samp_time == 0.5))
  biology <- record$inputs[record$inputs$sex == "female" &
    record$inputs$year %in% stock$years, ]
  key <- function(x) paste(x$year, x$age)
  lookup <- function(fleet, measure, target) {
    x <- biology[biology$fleet %in% fleet & biology$measure == measure, ]
    x$value[match(key(target), key(x))]
  }
  testthat::expect_equal(stock$obs$weight$obs,
    lookup("beginning-season population", "start_weight_at_age", stock$obs$weight))
  testthat::expect_equal(stock$obs$maturity$obs,
    0.5 * lookup("mature weight schedule", "mature_weight_at_age", stock$obs$maturity) /
      lookup("mid-season population", "weight_at_age", stock$obs$maturity))
  share <- stock$age_share("catch", "Hook-and-line", 2020)
  testthat::expect_equal(sum(share$share), 1)
  testthat::expect_setequal(share$sex, c("female", "male"))
  testthat::expect_true(all(share$sex_code == "3"))
  # Reported states are for comparisons/optional starts, never for observations.
  changed <- record
  changed$outputs$value[changed$outputs$measure == "numbers_at_age"] <-
    changed$outputs$value[changed$outputs$measure == "numbers_at_age"] * 10
  other <- .test_stock(changed)
  testthat::expect_identical(stock$obs, other$obs)
})
