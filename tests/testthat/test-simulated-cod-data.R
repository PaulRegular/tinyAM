test_that("simulated cod inputs retain complete fleet keys and explicit units", {
  expect_s3_class(sim_cod_obs, "tam_data")
  expect_invisible(check_tam_data(sim_cod_obs, quiet = TRUE))
  catch <- sim_cod_obs$catch
  keys <- catch[c("year", "age", "fleet")]
  expect_false(anyDuplicated(keys) > 0L)
  expect_setequal(unique(as.character(catch$fleet)), c("Gillnet", "Longline", "Trawl"))
  expected <- expand.grid(year = sort(unique(catch$year)),
    age = sort(unique(catch$age)), fleet = unique(as.character(catch$fleet)))
  key <- function(x) paste(x$year, x$age, x$fleet)
  expect_setequal(key(catch), key(expected))
  expect_true(all(is.na(catch$obs) | is.finite(catch$obs) & catch$obs >= 0))
  expect_identical(unique(catch$unit), "thousand fish")
  annual <- aggregate(obs ~ year + fleet, catch, sum)
  expect_true(all(annual$obs > 0))

  index <- sim_cod_obs$index
  expect_length(unique(index$survey), 2L)
  expect_true(all(index$samp_time == 0))
  expect_false(anyDuplicated(index[c("year", "age", "survey")]) > 0L)
  expect_identical(unique(index$unit), "thousand fish")
})

test_that("simulated cod biological inputs distinguish stock and fleet catch weights", {
  weight <- sim_cod_obs$weight
  stock <- weight[weight$type == "stock", ]
  catch <- weight[weight$type == "catch", ]
  expect_true(all(is.na(stock$fleet)))
  expect_setequal(unique(as.character(catch$fleet)),
    unique(as.character(sim_cod_obs$catch$fleet)))
  expect_identical(unique(weight$unit), "kg")
  expect_true(all(is.finite(weight$obs) & weight$obs > 0))
  expect_true(all(sim_cod_obs$maturity$obs >= 0 & sim_cod_obs$maturity$obs <= 1))
  for (fleet in unique(as.character(catch$fleet))) {
    rows <- catch[catch$fleet == fleet, ]
    expect_setequal(paste(rows$year, rows$age), paste(stock$year, stock$age))
  }
})

test_that("simulated cod truth conserves rates, removals and biological units", {
  population <- sim_cod_truth$population
  fleet <- sim_cod_truth$fleet
  key <- function(x) paste(x$year, x$age)
  position <- match(key(fleet), key(population))
  expect_false(anyNA(position))
  expect_equal(fleet$N, population$N[position], tolerance = 1e-12)
  empty <- fleet$N == 0
  expect_true(all(is.na(fleet$F_realized[empty])))
  expect_true(all(fleet$catch_realized[empty] == 0))
  expected_catch <- ifelse(empty, 0,
    fleet$N * fleet$F_realized / fleet$Z * -expm1(-fleet$Z))
  expect_equal(fleet$catch_realized, expected_catch, tolerance = 1e-10)
  expect_equal(fleet$yield, fleet$catch_realized * fleet$weight, tolerance = 1e-10)
  removals <- aggregate(catch_realized ~ year + age, fleet, sum)
  budget <- fleet$death_budget[match(key(removals), key(fleet))]
  expect_true(all(removals$catch_realized <= budget + 1e-10))
  total_f <- aggregate(F_realized ~ year + age, fleet, sum)
  expected_f <- total_f$F_realized[match(key(population), key(total_f))]
  populated <- population$N > 0
  expect_true(all(is.finite(population$F[populated]) & is.finite(population$M[populated])))
  expect_true(all(is.na(population$F[!populated]) & is.na(population$M[!populated])))
  expect_equal(population$Z[populated],
    population$M[populated] + expected_f[populated], tolerance = 1e-12)

  stock <- subset(sim_cod_obs$weight, type == "stock")
  maturity <- sim_cod_obs$maturity
  expected_biomass <- population$N * stock$obs[match(key(population), key(stock))]
  expect_equal(population$B, expected_biomass, tolerance = 1e-12)
  expect_equal(population$SSB,
    expected_biomass * maturity$obs[match(key(population), key(maturity))],
    tolerance = 1e-12)
  observed <- aggregate(obs ~ year + fleet, sim_cod_obs$catch, sum)
  monitoring <- sim_cod_truth$monitoring
  position <- match(paste(observed$year, observed$fleet),
    paste(monitoring$year, monitoring$fleet))
  expect_equal(1000 * observed$obs, monitoring$reported_total[position], tolerance = 1e-10)
})

test_that("survey truth refers to the sampled region without whole-area expansion", {
  survey <- sim_cod_truth$surveys
  population <- sim_cod_truth$population
  position <- match(paste(survey$year, survey$age), paste(population$year, population$age))
  expect_equal(survey$N_total, population$N[position], tolerance = 1e-12)
  expect_true(all(survey$N_region >= 0 & survey$N_region <= survey$N_total))
  populated <- survey$N_total > 0
  expect_true(all(is.na(survey$availability[!populated])))
  expect_true(all(is.finite(survey$availability[populated])))
  expect_equal(survey$availability[populated],
    survey$N_region[populated] / survey$N_total[populated], tolerance = 1e-12)
  expect_equal(survey$index_expected, survey$N_region * survey$selection, tolerance = 1e-12)
  footprint <- sim_cod_truth$footprint
  for (name in unique(survey$survey)) {
    expect_true(any(footprint[[name]]))
    expect_true(any(!footprint[[name]]))
  }
  expect_true(any(footprint$RV & footprint$Sentinel))
  expect_gt(mean(footprint$y[footprint$Sentinel]), mean(footprint$y[footprint$RV]))
})

test_that("active fleet targets respect their accessible vulnerable populations", {
  allocation <- sim_cod_truth$allocations
  expect_false(anyDuplicated(allocation[c("year", "fleet")]) > 0L)
  expect_true(all(is.finite(allocation$base_share) & allocation$base_share > 0))
  expect_true(all(is.finite(allocation$fishery_share) & allocation$fishery_share > 0))
  expect_true(all(allocation$fishery_share <= allocation$base_share))
  expect_true(all(allocation$access_fraction > 0 & allocation$access_fraction <= 1))
  expect_equal(allocation$target,
    allocation$base_share * allocation$accessible_vulnerable_deaths, tolerance = 1e-12)
  expect_equal(allocation$target,
    allocation$fishery_share * allocation$global_vulnerable_deaths, tolerance = 1e-12)
  annual <- sim_cod_truth$fishery_annual
  position <- match(paste(annual$year, annual$fishery),
    paste(allocation$year, allocation$fleet))
  expect_false(anyNA(position))
  expect_equal(annual$reference_target, allocation$target[position], tolerance = 1e-12)
  expect_true(all(annual$n > 0))
})
