source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

testthat::test_that("BC Yellowtail uses matched age samples, timing and native biology", {
  database <- read_database()
  record <- read_assessment("dfo_yellowtail_bc_2024", database)
  testthat::local_mocked_bindings(fit_tam = function(...) stop("Do not fit during preparation"),
    .package = "tinyAM")
  stock <- .test_stock(record)
  testthat::expect_null(stock$fit)
  testthat::expect_identical(stock$years, 1979:2018)
  testthat::expect_identical(stock$ages, 0:45)
  testthat::expect_identical(stock$dat$F_settings$mean_ages, 9:20)
  testthat::expect_false(stock$dat$has_fleet)
  testthat::expect_equal(stock$dat$rec$sd, 0.9)
  testthat::expect_equal(exp(stock$dat$log_mu_supplied_m[1]), 0.08 / sqrt(1.09))
  testthat::expect_equal(stock$dat$M_terms[[1]]$sd, sqrt(log(1.09)))

  catch <- stock$obs$catch
  catch_weights <- subset(stock$obs$weight, type == "catch")
  weight_row <- match(with(catch, paste(year, age)), with(catch_weights, paste(year, age)))
  catch$tonnes <- catch$obs * catch_weights$obs[weight_row]
  annual <- aggregate(tonnes ~ year, catch, sum)
  testthat::expect_identical(as.integer(annual$year), setdiff(1979:2018, c(1981L, 1984L)))
  testthat::expect_equal(annual$tonnes, stock$catch_total$value[match(annual$year, stock$catch_total$year)])
  testthat::expect_true(all(is.na(catch$obs[catch$age == 0L])))
  testthat::expect_true(any(catch$obs == 0, na.rm = TRUE))

  index <- stock$obs$index
  index_years <- list(NMFS_TRIENNIAL = c(1983L, 1989L, 1992L, 1995L, 1998L, 2001L),
    QCS_SYNOPTIC = c(2003L, 2004L, 2005L, 2009L, 2011L, 2013L, 2015L, 2017L),
    WCVI_SYNOPTIC = c(2006L, 2008L, 2010L, 2012L, 2016L, 2018L))
  testthat::expect_setequal(unique(index$survey), names(index_years))
  native_month <- c(NMFS_TRIENNIAL = 7L, QCS_SYNOPTIC = 7L, WCVI_SYNOPTIC = 5L)
  for (survey in names(index_years)) {
    rows <- index$survey == survey
    testthat::expect_identical(sort(unique(as.integer(index$year[rows]))), index_years[[survey]])
    testthat::expect_equal(unique(index$samp_time[rows]), (native_month[[survey]] - 1) / 12)
  }
  index_annual <- aggregate(tonnes ~ year + survey, stock$index_numbers, sum)
  testthat::expect_equal(index_annual$tonnes, stock$index_total$value[
    match(stock$survey_year(index_annual), stock$survey_year(stock$index_total))])

  stock_weights <- subset(stock$obs$weight, type == "stock")
  maturity <- stock$obs$maturity
  weight_row <- match(with(maturity, paste(year, age)), with(stock_weights, paste(year, age)))
  female <- stock$female_reproductive_weight
  testthat::expect_true(all(maturity$spawn_time == 0))
  testthat::expect_true(all(maturity$obs[maturity$age < 9] == 0))
  testthat::expect_equal(maturity$obs * stock_weights$obs[weight_row],
    0.5 * female$value[match(maturity$age, female$age)])

  n <- subset(stock$comparison_outputs, type == "population" & measure == "numbers_at_age")
  f <- subset(stock$comparison_outputs, type == "mortality" & measure == "fishing_mortality_at_age")
  testthat::expect_equal(nrow(n), 40L * 46L)
  testthat::expect_equal(nrow(f), 40L * 46L)
  testthat::expect_true(all(n$sex == "combined" & f$sex == "combined"))
  testthat::expect_equal(sum(stock$comparison_outputs$measure == "SSB"), 40L)
  ssb_age <- subset(stock$comparison_outputs, measure == "mature_biomass_at_age")
  ssb <- subset(stock$comparison_outputs, measure == "SSB")
  total <- tapply(ssb_age$value, ssb_age$year, sum)
  testthat::expect_equal(nrow(ssb_age), 40L * 46L)
  testthat::expect_true(all(ssb_age$sex == "female"))
  testthat::expect_true(all(abs(total[as.character(ssb$year)] / ssb$value - 1) < 1e-5))
})
