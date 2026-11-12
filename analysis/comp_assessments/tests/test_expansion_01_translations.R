source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
database <- read_database()

testthat::test_that("the expansion register keeps all twenty distinct source targets visible", {
  batch <- read.csv("analysis/comp_assessments/expansion_batch_01.csv", stringsAsFactors = FALSE)
  testthat::expect_equal(nrow(batch), 20L)
  testthat::expect_false(anyDuplicated(batch$stock_id) > 0L)
  testthat::expect_false(anyDuplicated(batch$assessment_id) > 0L)
  testthat::expect_true(all(batch$assessment_id %in% database$assessments$assessment_id))
  # The batch register records its original source review; do not relabel it
  # when the selected local snapshot is refreshed.
  testthat::expect_true(all(batch$data_revision == "6f08c06c49e4b73591427a1495bfb4568caf49f0"))
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
    nafo_redfish_3m_2024 = "annual age-specific inputs are unavailable",
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

testthat::test_that("sablefish expands six catch streams and four surveys without filling unsampled ages", {
  record <- read_assessment("noaa_sablefish_pacific_coast_2025", database)
  before <- record
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_identical(record, before)
  testthat::expect_equal(stock$years, 1983:2024)
  testthat::expect_equal(stock$ages, 0:30)
  testthat::expect_true(all(abs(stock$dat$log_mu_supplied_m + 2.63061660) < 1e-12))
  testthat::expect_length(stock$dat$M_terms, 1L)
  testthat::expect_equal(stock$dat$M_terms[[1L]]$sd, 0.19036084)
  testthat::expect_equal(pnorm(log(0.04), -2.63061660, stock$dat$M_terms[[1L]]$sd),
    0.001, tolerance = 1e-6)
  testthat::expect_length(stock$dat$M_terms[[1L]]$groups, 1L)

  catch <- stock$obs$catch
  weights <- stock$obs$weight
  fleets <- c("Trawl", "Hook-and-line", "Pot",
    "Trawl discard", "Hook-and-line discard", "Pot discard")
  testthat::expect_setequal(catch$fleet, fleets)
  testthat::expect_equal(nrow(catch), length(stock$years) * length(stock$ages) * length(fleets))
  testthat::expect_equal(stock$dat$F_settings$process, "off")
  temporal <- Filter(function(term) term$type == "ar1", stock$dat$F_terms)
  testthat::expect_length(temporal, 1L)
  testthat::expect_length(temporal[[1L]]$groups, length(fleets))
  testthat::expect_setequal(vapply(stock$dat$F_scopes, `[[`, character(1), "name"), fleets)
  testthat::expect_true(all(catch$F_age <= catch$age))
  testthat::expect_equal(sort(unique(catch$year[catch$retained_period == "late"])), 2002:2024)
  retained_shape <- c("fixed_gear_young", "fixed_gear_rise", "fixed_gear_decline", "fixed_gear_shift")
  testthat::expect_equal(subset(catch, fleet == "Hook-and-line", retained_shape),
    subset(catch, fleet == "Pot", retained_shape), ignore_attr = TRUE)
  testthat::expect_true(all(subset(catch, !fleet %in% c("Hook-and-line", "Pot"),
    retained_shape) == 0))
  groups <- function(gear, yr, ages) {
    subset(catch, fleet == paste(gear, "discard") & year == yr & age %in% ages)$discard_group
  }
  early_hkl <- groups("Hook-and-line", 2010, 0:6)
  testthat::expect_identical(match(early_hkl, unique(early_hkl)), c(1L, 2L, 3L, 3L, 4L, 4L, 4L))
  testthat::expect_identical(groups("Hook-and-line", 2011, c(0, 1, 2, 3, 4, 6, 30)),
    groups("Hook-and-line", 2019, c(0, 2, 3, 4, 5, 6, 30)))
  testthat::expect_identical(groups("Pot", 2010, 0:7),
    c("0", "early 1", rep("tail", 6)))
  testthat::expect_identical(groups("Pot", 2011, 0:7),
    c("0", "late 1", "late 2", rep("tail", 5)))
  testthat::expect_true(any(subset(catch, fleet == "Hook-and-line discard" & age == 0)$obs > 0,
    na.rm = TRUE))
  testthat::expect_true(any(subset(catch, fleet == "Pot discard" & age == 0)$obs > 0,
    na.rm = TRUE))
  catch_weights <- subset(weights, type == "catch")
  key <- function(x) paste(x$year, x$age, x$fleet)
  catch$tonnes <- catch$obs * catch_weights$obs[match(key(catch), key(catch_weights))]
  ages_observed <- aggregate(is.finite(obs) ~ year + fleet, catch, sum)[[3L]]
  testthat::expect_true(all(ages_observed %in% c(0L, length(stock$ages))))
  observed <- aggregate(tonnes ~ year + fleet, catch, sum)
  reconstruction <- stock$catch_reconstruction
  covered <- is.finite(reconstruction$reconstructed_tonnes)
  fleet_year <- function(x) paste(x$year, x$fleet)
  testthat::expect_setequal(fleet_year(observed), fleet_year(reconstruction[covered, ]))
  testthat::expect_equal(as.integer(table(factor(observed$fleet, levels = fleets))),
    c(40L, 38L, 36L, 18L, 18L, 18L))
  testthat::expect_equal(observed$tonnes,
    reconstruction$source_tonnes[match(fleet_year(observed), fleet_year(reconstruction))])
  testthat::expect_true(all(is.na(subset(catch, fleet == "Pot" & year == 2014)$obs)))
  testthat::expect_true(all(is.finite(subset(catch, fleet == "Pot discard" & year == 2014)$obs)))

  index <- stock$obs$index
  survey_years <- list(`Early Triennial survey` = c(1983L, 1989L, 1992L),
    `Late Triennial survey` = c(1995L, 2004L), `NWFSC Slope Survey` = 1998:2002,
    WCGBTS = c(2003:2019, 2021:2024))
  testthat::expect_setequal(index$survey, names(survey_years))
  for (survey in names(survey_years)) {
    testthat::expect_equal(sort(unique(index$year[index$survey == survey])), survey_years[[survey]])
  }
  testthat::expect_true(all(index$unit == "native survey index" & index$samp_time == 0.5))
  testthat::expect_setequal(vapply(stock$dat$q_scope_components, `[[`, character(1), "name"), names(survey_years))
  testthat::expect_gt(ncol(stock$dat$q_modmat), 0L)
  index_mass <- stock$index_reconstruction
  scored <- is.finite(index_mass$reconstructed_biomass_index)
  testthat::expect_equal(sum(scored), 31L)
  testthat::expect_equal(index_mass$reconstructed_biomass_index[scored],
    index_mass$source_biomass_index[scored])
  wcgbts <- subset(index, survey == "WCGBTS")
  testthat::expect_identical(levels(wcgbts$wcgbts_order), c("0", "2+", "1"))
  testthat::expect_true(all(wcgbts$wcgbts_order[wcgbts$age >= 2L] == "2+"))
  # Fit on the survey's original number-equivalent scale, with no official q adjustment.
  source_share <- stock$age_share("index", "WCGBTS", 2024L)
  w <- stock$weight_at(2024L, source_share$age, source_share$sex, "WCGBTS", "weight_at_age")
  total <- subset(index_mass, survey == "WCGBTS" & year == 2024L)$source_biomass_index
  age0 <- subset(wcgbts, year == 2024L & age == 0L)$obs
  testthat::expect_equal(age0, total * sum(source_share$share[source_share$age == 0L]) /
    sum(source_share$share * w))

  stock_weights <- subset(weights, type == "stock")
  biology <- record$inputs
  weight_at <- function(fleet, measure, sex, target) {
    x <- biology[biology$fleet %in% fleet & biology$measure == measure & biology$sex %in% sex, ]
    source_year <- ifelse(target$year < 2003L, -1890L, target$year)
    x$value[match(paste(source_year, target$age), paste(x$year, x$age))]
  }
  start_weight <- function(sex) weight_at("beginning-season population", "start_weight_at_age",
    sex, stock_weights)
  testthat::expect_equal(stock_weights$obs, (start_weight("female") + start_weight("male")) / 2)
  maturity <- stock$obs$maturity
  testthat::expect_true(all(maturity$spawn_time == 0))
  testthat::expect_equal(stock_weights$obs * maturity$obs,
    0.5 * weight_at("mature weight schedule", "mature_weight_at_age", "female", maturity))

  share <- stock$age_share("catch", "Hook-and-line", 2020)
  testthat::expect_equal(sum(share$share), 1)
  testthat::expect_setequal(share$sex, c("female", "male"))
  testthat::expect_true(all(share$sex_code == "3"))
  rates <- subset(stock$comparison_outputs, measure == "fishing_mortality_at_age")
  components <- subset(rates, !is.na(fleet))
  testthat::expect_setequal(components$fleet, fleets)
  summed <- aggregate(value ~ year + age, components, sum)
  total <- subset(rates, is.na(fleet))
  testthat::expect_equal(summed$value,
    total$value[match(paste(summed$year, summed$age), paste(total$year, total$age))])
  # Reported states are for comparisons/optional starts, never for observations.
  changed <- record
  changed$outputs$value[changed$outputs$measure == "numbers_at_age"] <-
    changed$outputs$value[changed$outputs$measure == "numbers_at_age"] * 10
  other <- .test_stock(changed)
  testthat::expect_identical(stock$obs, other$obs)
})
