test_that("simulated cod inputs retain complete fleet keys and explicit units", {
  expect_equal(sim_cod_truth$provenance$seed, 17)
  expect_s3_class(sim_cod_obs, "tam_data")
  expect_invisible(check_tam_data(sim_cod_obs, quiet = TRUE))
  catch <- sim_cod_obs$catch
  expect_false(anyDuplicated(catch[c("year", "age", "fleet")]) > 0L)
  expect_setequal(as.character(catch$fleet), c("Gillnet", "Longline", "Trawl"))
  expect_setequal(catch$age, 1:20)
  expected <- expand.grid(year = unique(catch$year), age = 1:20,
    fleet = unique(catch$fleet))
  expect_setequal(with(catch, paste(year, age, fleet)),
    with(expected, paste(year, age, fleet)))
  expect_true(all(is.na(catch$obs) | is.finite(catch$obs) & catch$obs >= 0))
  expect_identical(unique(catch$unit), "thousand fish")

  index <- sim_cod_obs$index
  expect_setequal(as.character(index$survey), c("RV", "Sentinel"))
  expect_true(all(index$samp_time == 0))
  expect_false(anyDuplicated(index[c("year", "age", "survey")]) > 0L)
  expect_identical(unique(index$unit), "thousand fish")
  expected <- expand.grid(year = unique(catch$year), age = 1:20,
    survey = unique(index$survey))
  expect_setequal(with(index, paste(year, age, survey)),
    with(expected, paste(year, age, survey)))

  precision <- sim_cod_truth$survey_precision
  position <- match(with(index, paste(year, age, survey)),
    with(precision, paste(year, age, survey)))
  expect_false(anyNA(position))
  status <- precision$estimate_status[position]
  expect_true(all(status[!is.na(index$obs)] %in% c("ok", "partial_age")))
  expect_true(all(is.na(index$obs[status %in% c("no_sets", "unsampled_strata", "missing_lengths")])))
})

test_that("simulated cod weights and maturity reproduce biomass and SSB", {
  weight <- sim_cod_obs$weight
  stock <- subset(weight, type == "stock")
  catch <- subset(weight, type == "catch")
  expect_true(all(is.na(stock$fleet)))
  expect_setequal(as.character(catch$fleet), as.character(sim_cod_obs$catch$fleet))
  expect_identical(unique(weight$unit), "kg")
  expect_true(all(is.finite(weight$obs) & weight$obs > 0))
  expect_true(all(sim_cod_obs$maturity$obs >= 0 & sim_cod_obs$maturity$obs <= 1))
  key <- function(x) paste(x$year, x$age)
  for (f in unique(catch$fleet)) {
    expect_setequal(key(catch[catch$fleet == f, ]), key(stock))
  }

  biology <- sim_cod_truth$catch_biology
  fleet_key <- function(x) paste(x$fleet, x$year, x$age)
  position <- match(fleet_key(catch), fleet_key(biology))
  expect_false(anyNA(position))
  estimated <- biology$mean_weight_kg[position]
  supported <- biology$weight_status[position] == "estimated"
  expect_identical(!is.na(estimated), supported)
  expect_true(any(supported) && any(!supported))
  expect_equal(catch$obs[supported], estimated[supported])
  expect_equal(catch$obs[!supported], stock$obs[match(key(catch[!supported, ]), key(stock))])

  population <- sim_cod_truth$population
  maturity <- sim_cod_obs$maturity
  expected_biomass <- population$N * stock$obs[match(key(population), key(stock))]
  expect_equal(population$B, expected_biomass, tolerance = 1e-12)
  expect_equal(population$SSB,
    expected_biomass * maturity$obs[match(key(population), key(maturity))],
    tolerance = 1e-12)
  expected_ssb <- aggregate(SSB ~ year, population, sum)
  annual <- sim_cod_truth$annual
  expect_equal(annual$SSB, expected_ssb$SSB[match(annual$year, expected_ssb$year)])
  first <- which.min(annual$year)
  last <- which.max(annual$year)
  expect_lt(annual$N[last], annual$N[first])
  expect_lt(annual$SSB[last], annual$SSB[first])
})

test_that("simulated cod fits can pool observations and retain older biology", {
  dat <- prepare_tam(sim_cod_obs, ages = 1:14,
    M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
    index_settings = list(q_form = ~ survey))
  expect_equal(dat$ages, 1:14)
  expect_equal(dat$plus_ages, 14:20)
  for (name in c("catch", "index")) {
    group <- if (name == "catch") "fleet" else "survey"
    raw <- subset(sim_cod_obs[[name]], age >= 14)
    expected <- aggregate(raw["obs"], raw[c("year", group)], function(x) {
      if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
    })
    pooled <- subset(dat$obs[[name]], age == 14)
    keys <- function(x) paste(x$year, x[[group]])
    expect_setequal(keys(pooled), keys(expected))
    expect_equal(pooled$obs, expected$obs[match(keys(pooled), keys(expected))])
  }
  stock <- subset(sim_cod_obs$weight, type == "stock" & age >= 14)
  maturity <- subset(sim_cod_obs$maturity, age >= 14)
  expected_weight <- with(stock, tapply(obs, list(year, age), identity))
  expected_maturity <- with(maturity, tapply(obs, list(year, age), identity))
  expect_equal(as.vector(dat$W_plus_input), as.vector(expected_weight))
  expect_equal(as.vector(dat$P_plus_input), as.vector(expected_maturity))
  for (f in dat$fleets) {
    catch_weight <- subset(sim_cod_obs$weight, type == "catch" & fleet == f & age >= 14)
    expected <- with(catch_weight, tapply(obs, list(year, age), identity))
    expect_equal(as.vector(dat$W_catch_plus_fleet[[f]]), as.vector(expected))
  }
})

test_that("simulated cod truth conserves mortality, removals and yield", {
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
  budget <- fleet$total_deaths[match(key(removals), key(fleet))]
  expect_true(all(removals$catch_realized <= budget + 1e-10))

  total_f <- aggregate(F_realized ~ year + age, fleet, sum)
  expected_f <- total_f$F_realized[match(key(population), key(total_f))]
  populated <- population$N > 0
  expect_true(all(is.finite(population$F[populated]) & is.finite(population$M[populated])))
  expect_true(all(is.na(population$F[!populated]) & is.na(population$M[!populated])))
  expect_equal(population$F, expected_f, tolerance = 1e-12)
  expect_equal(population$Z[populated],
    population$M[populated] + expected_f[populated], tolerance = 1e-12)
  expect_lte(abs(mean(subset(population, age %in% 5:10)$M) - 0.30), 0.03)
  periods <- sim_cod_truth$checks$realized_M_by_period
  expect_setequal(periods$period, c("Before", "After"))
  adult <- subset(population, age %in% 5:10)
  expected_means <- c(Before = mean(adult$M[adult$year < 11]),
    After = mean(adult$M[adult$year >= 11]))
  expect_equal(periods$mean, unname(expected_means[periods$period]))
  expect_true(all(abs(periods$mean - 0.30) <= 0.03))

  terminal <- subset(population, age == 20)
  loss <- sim_cod_truth$checks$terminal_survivors
  i <- match(loss$year, terminal$year)
  expect_equal(loss$N, with(terminal, N * exp(-Z))[i])
  expect_equal(loss$SSB, with(terminal, SSB * exp(-Z))[i])
})

test_that("fleet F walks follow the seeded design and only trawl is reduced", {
  annual <- sim_cod_truth$fishery_annual
  provenance <- sim_cod_truth$provenance
  expect_true(all(annual$target_type == "F_full"))
  expect_true(all(annual$n > 0 & annual$requested_target > 0))
  expect_equal(provenance$fishing_seed, 1701L)
  expect_equal(provenance$fishing_log_step_sd, 0.10)
  expect_equal(provenance$trawl_reduction_year, 11L)
  starting_F <- c(Gillnet = 0.06, Longline = 0.05, Trawl = 0.07)
  increments <- withr::with_seed(1701, matrix(stats::rnorm(57, sd = 0.10), ncol = 3),
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  expected <- sweep(exp(apply(rbind(0, increments), 2, cumsum)), 2, starting_F, "*")
  i <- cbind(annual$year, match(annual$fishery, names(starting_F)))
  expect_equal(annual$F_full_unreduced, expected[i], tolerance = 1e-12)
  reduction <- ifelse(annual$fishery == "Trawl" & annual$year >= 11, 0.2, 1)
  expect_equal(annual$F_full_reference, annual$F_full_unreduced * reduction)
  expect_equal(annual$F_full_reference, annual$requested_target)

  fleet <- sim_cod_truth$fleet
  targets <- annual$requested_target[match(paste(fleet$year, fleet$fleet),
    paste(annual$year, annual$fishery))]
  expect_equal(fleet$F_reference, targets * fleet$selectivity, tolerance = 1e-12)
  ref <- aggregate(F_reference ~ year + age, fleet, sum)
  population <- sim_cod_truth$population
  expected_f <- ref$F_reference[match(paste(population$year, population$age),
    paste(ref$year, ref$age))]
  expect_equal(population$Z, population$M_reference + expected_f, tolerance = 1e-12)
  # The first M draw still uses the population seed, independent of the F walks.
  expect_equal(provenance$reference_M, 0.28)
  first_M <- withr::with_seed(17, exp(log(provenance$reference_M) +
    stats::rnorm(1, sd = 0.04 / sqrt((1 - 0.9^2) * (1 - 0.5^2)))),
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  expect_equal(population$M_reference[population$year == 1 & population$age == 1],
    first_M, tolerance = 1e-12)
})

test_that("simulated catches match reported totals", {
  catch <- sim_cod_obs$catch
  observed <- aggregate(catch["obs"], catch[c("year", "fleet")], function(x) {
    if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  })
  monitoring <- sim_cod_truth$monitoring
  position <- match(paste(observed$year, observed$fleet),
    paste(monitoring$year, monitoring$fleet))
  expect_false(anyNA(position))
  expect_equal(1000 * observed$obs, monitoring$estimated_total[position], tolerance = 1e-10)
  supported <- monitoring$status[position] %in% c("estimated", "zero_reported_catch")
  expect_identical(!is.na(observed$obs), supported)
  expect_equal(1000 * observed$obs[supported],
    monitoring$reported_total[position[supported]], tolerance = 1e-10)
})

test_that("Gillnet and Sentinel share the generating normal curve", {
  provenance <- sim_cod_truth$provenance
  peak <- provenance$gillnet_peak_age
  width <- provenance$gillnet_width
  selectivity <- function(age) {
    stats::dnorm(age, peak, width) / stats::dnorm(peak, peak, width)
  }
  curves <- subset(sim_cod_truth$selectivity, component %in% c("Gillnet", "Sentinel"))
  expect_equal(curves$selectivity, selectivity(curves$age))

  fleet <- subset(sim_cod_truth$fleet, fleet == "Gillnet")
  survey <- subset(sim_cod_truth$surveys, survey == "Sentinel")
  expect_equal(fleet$selectivity, selectivity(fleet$age))
  expect_equal(survey$selection, selectivity(survey$age))
})

test_that("survey truth describes distinct, stable footprints", {
  survey <- sim_cod_truth$surveys
  population <- sim_cod_truth$population
  position <- match(paste(survey$year, survey$age), paste(population$year, population$age))
  expect_equal(survey$N_total, population$N[position], tolerance = 1e-12)
  expect_true(all(survey$N_region >= 0 & survey$N_region <= survey$N_total))
  populated <- survey$N_total > 0
  expect_true(all(is.na(survey$availability[!populated])))
  expect_equal(survey$availability[populated],
    survey$N_region[populated] / survey$N_total[populated], tolerance = 1e-12)
  expect_equal(survey$index_expected, survey$N_region * survey$selection, tolerance = 1e-12)

  footprint <- sim_cod_truth$footprint
  expect_equal(nrow(footprint), 20L * 20L)
  for (name in unique(survey$survey)) {
    expect_true(any(footprint[[name]]) && any(!footprint[[name]]))
  }
  expect_lt(mean(footprint$depth[footprint$Sentinel]), mean(footprint$depth[footprint$RV]))
  area <- vapply(c("RV", "Sentinel"), function(name) {
    sum(footprint$area_km2[footprint[[name]]])
  }, numeric(1))
  monitoring <- sim_cod_truth$survey_monitoring
  sets <- with(monitoring, tapply(n_sets, survey, sum))
  density <- sets[names(area)] / area
  expect_gt(density[["Sentinel"]], density[["RV"]])
  well_populated <- subset(survey, N_total >= 100)
  spread <- vapply(split(well_populated$availability, well_populated$survey),
    function(x) diff(range(x)), numeric(1))
  expect_true(all(is.finite(spread) & spread <= .005))
})
