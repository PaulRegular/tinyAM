source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/bc_canary_obs.R")
source("analysis/comp_assessments/tests/helper_stock.R")
record <- read_assessment("dfo_canary_bc_2022")
original <- record
x <- .test_stock(record)
stopifnot(identical(record, original), is.null(x$fit),
  identical(x$years, 1977:2022), identical(x$ages, 1:60),
  x$settings$N_settings$process == "off", x$settings$N_settings$init == "free",
  x$settings$rec_settings$sd == .9,
  identical(all.vars(x$settings$rec_settings$mu_form), character()),
  x$settings$F_settings$process == "iid", x$settings$M_settings$process == "off",
  setequal(x$obs$index$survey, c("QCS_SYNOPTIC", "WCVI_SYNOPTIC", "NMFS_TRIENNIAL")),
  max(abs(x$obs$index$samp_time - 5/12)) < 1e-12,
  all(x$obs$maturity$obs[x$obs$maturity$age < 10] == 0),
  all(x$obs$maturity$spawn_time == 0),
  all(is.na(x$obs$catch$obs[x$obs$catch$year >= 2018])),
  !any(x$inputs$type == "M"), all(x$comparison_outputs$sex == "female"),
  !any(x$comparison_outputs$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "recruitment")),
  x$comparison_definitions$ssb$status == "approximate",
  identical(x$comparison_aggregates, "ssb"),
  sum(x$dat$is_proj) == 1L, max(x$dat$years) == 2023L,
  is.finite(tinyAM::nll_fun(tinyAM::make_par(x$dat), x$dat)))
tinyAM::check_obs(x$obs)

# The biomass-to-number denominator includes both sexes before selection.
fixture <- record$inputs[rep(1L, 10L), ]
fixture$assessment_id <- "fixture"
fixture$type <- c(rep("catch", 6), rep("catch_weight", 4))
fixture$measure <- c(rep("proportion_at_age", 4), rep("total_biomass", 2), rep("weight_at_age", 4))
fixture$sex <- c("female", "female", "male", "male", "combined", "combined",
                 "female", "female", "male", "male")
fixture$age <- c(1, 2, 1, 2, NA, NA, 1, 2, 1, 2)
fixture$year <- c(rep(2000L, 6), rep(NA_integer_, 4))
fixture$value <- c(.1, .2, .3, .4, .9, .1, 1, 2, 3, 4)
fixture$unit <- c(rep("proportion", 4), "t", "t", rep("kg/fish", 4))
fixture$basis <- c(rep("proportion_numbers", 4), rep("biomass", 2), rep("kg_per_fish", 4))
fixture$fleet <- c(rep("TRAWL_FISHERY", 5), "OTHER_FISHERY", rep("TRAWL_FISHERY", 4))
fixture$survey <- NA_character_
fixture$observation_id <- c(rep("sample:1", 4), rep(NA_character_, 6))
original_fixture <- fixture
v <- .canary_female_numbers(fixture, "catch")
stopifnot(identical(fixture, original_fixture), nrow(v) == 2L,
  max(abs(v$value - c(100/3, 200/3))) < 1e-10,
  all(v$sex == "female"), all(v$fleet == "POOLED_CATCH"),
  all(v$source_type == "derived_tinyAM_observation"))
testthat::expect_error(.canary_female_numbers(fixture[-6, ], "catch"), "matching original")
testthat::expect_error(.canary_female_numbers(fixture[-10, ], "catch"), "positive weights")
testthat::expect_error(.canary_female_numbers(fixture[-c(3, 4), ], "catch"), "joint-sex composition")
bad <- fixture; bad$value[1] <- -1
testthat::expect_error(.canary_female_numbers(bad, "catch"), "joint-sex composition")

# Actual observations reproduce the joint-sex calculation without a female 0.5 factor.
sample <- record$inputs[record$inputs$measure == "proportion_at_age" &
  record$inputs$observation_id == "age:1", ]
sample <- sample[!is.na(sample$observation_id), ]
cw <- record$inputs[record$inputs$type == "catch_weight" &
  record$inputs$fleet == "TRAWL_FISHERY", ]
key <- function(z) paste(z$sex, z$age)
denominator <- sum(sample$value * cw$value[match(key(sample), key(cw))])
total <- sum(record$inputs$value[record$inputs$type == "catch" &
  record$inputs$measure == "total_biomass" & record$inputs$year == 1977]) * 1000
female <- sample[sample$sex == "female", ]
observed <- x$obs$catch[x$obs$catch$year == 1977, ]
stopifnot(max(abs(observed$obs - total/1000*female$value[match(observed$age, female$age)]/denominator)) < 1e-8)
cat("Canary preparation, joint-sex conversion, missing years, comparison limits and no-fit checks passed.\n")
