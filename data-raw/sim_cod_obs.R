# Run from the tinyAM root with the SimSurvey revision documented in README.md.
# Only the compact prepared data are saved. Redirect console output to keep a run log.
library(SimSurvey)
library(data.table)

revision <- "29a1b03350c274f88dac108fb95f9446eac79233"
installed_revision <- utils::packageDescription("SimSurvey", fields = "RemoteSha")
if (is.na(installed_revision)) installed_revision <- Sys.getenv("TINYAM_SIMSURVEY_REVISION")
stopifnot(identical(installed_revision, revision))
output <- Sys.getenv("TINYAM_SIM_COD_OUTPUT",
  "data-raw/cache/sim_cod/seed17-simsurvey-biology/prepared-data.rds")
if (file.exists(output)) stop("Choose a fresh TINYAM_SIM_COD_OUTPUT path.")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Draw the three fishing trajectories separately from M and recruitment.
years <- 1:20
trawl_reduction_year <- 11L
set.seed(1701)
F_full_unreduced <- lapply(c(Gillnet = 0.06, Longline = 0.05, Trawl = 0.07), function(f) {
  f * exp(c(0, cumsum(stats::rnorm(length(years) - 1L, sd = 0.10))))
})
F_full <- F_full_unreduced
F_full$Trawl[years >= trawl_reduction_year] <- 0.2 * F_full$Trawl[years >= trawl_reduction_year]

# Reuse these curves and annual rates when generating Z and fleet targets.
gillnet <- function(age) plogis(1.4 * (age - 5)) * plogis(0.8 * (10.5 - age))
longline <- sim_logistic(k = 1.6, x50 = 3.6)
trawl <- sim_logistic(k = 1.5, x50 = 4.3)
rv <- sim_logistic(k = 1.8, x50 = 2.5)

# Rotate the depth gradient so the shallow edge represents the northern coast.
grid <- make_grid(res = c(14, 14))
footprint <- as.data.frame(grid)
footprint[c("x", "y")] <- list(footprint$y, -footprint$x)
footprint <- transform(footprint,
  Gillnet = y >= 40 & depth >= 25 & depth <= 250,
  Longline = y >= 0 & depth >= 25 & depth <= 350,
  Trawl = y < 70 & depth >= 50 & depth <= 550,
  RV = y < 70 & depth >= 50 & depth <= 550,
  Sentinel = y >= 40 & depth >= 25 & depth <= 250)
grid <- stars::st_as_stars(footprint, crs = sf::st_crs(grid))

# Default recruitment, growth and stock weight. Share one spatial field so
# the fraction available within each survey footprint changes little over time.
# Reference M of 0.28 gives realized M near 0.30 once catches are allocated.
reference_M <- 0.28
set.seed(17)
pop <- sim_abundance(
  years = years,
  Z = function(ages, years) {
    M <- sim_Z(log_mean = log(reference_M), log_sd = 0.04)(ages, years)
    M + outer(gillnet(ages), F_full$Gillnet) +
      outer(longline(ages), F_full$Longline) + outer(trawl(ages), F_full$Trawl)
  },
  maturity = sim_maturity(k = 1, x50 = 5), maturity_by = "age"
) |>
  sim_distribution(grid = grid,
    ays_covar = sim_ays_covar(phi_age = 0, phi_year = 0, group_ages = 1:20, group_years = 1:20)) |>
  round_sim()

# Fisheries draw from the same death budget. Observe catches through reporting
# and random age samples from a subset of operations.
message("Simulating the three fisheries")
fishery <- sim_fisheries(pop,
  Gillnet = sim_fishery(target = F_full[["Gillnet"]], target_type = "F_full",
    selectivity = gillnet, subset_cells = Gillnet,
    catch_model = sim_capture(q = 0.02, capacity_per_gear = 150),
    gear_amount = 60, duration = 24, gear_unit = "net", time_unit = "hour",
    sampling_coverage = 0.30, lengths_cap = 150, ages_cap = 150,
    age_sampling = "random", reporting_cv = 0.05, maturity_sampling = "none"),
  Longline = sim_fishery(target = F_full[["Longline"]], target_type = "F_full",
    selectivity = longline, subset_cells = Longline,
    catch_model = sim_capture(q = 0.001, capacity_per_gear = 1),
    gear_amount = 4000, duration = 24, gear_unit = "hook", time_unit = "hour",
    sampling_coverage = 0.30, lengths_cap = 150, ages_cap = 150,
    age_sampling = "random", reporting_cv = 0.05, maturity_sampling = "none"),
  Trawl = sim_fishery(target = F_full[["Trawl"]], target_type = "F_full",
    selectivity = trawl, subset_cells = Trawl,
    catch_model = sim_capture(q = 2, capacity_per_gear = 10000),
    gear_amount = 1, duration = 4, gear_unit = "trawl", time_unit = "hour",
    sampling_coverage = 0.30, lengths_cap = 150, ages_cap = 150,
    age_sampling = "random", reporting_cv = 0.05, maturity_sampling = "none")
) |>
  run_catch() |>
  catch_error()
fishery <- fishery$fishery

# Both surveys observe the original population at the start of each year.
message("Simulating the offshore RV and inshore Sentinel surveys")
surveys <- sim_surveys(pop,
  RV = sim_survey(q = rv,
    custom_sets = sim_sets(pop, subset_cells = RV, set_den = 2 / 1000),
    lengths_cap = 100, ages_cap = 3, age_length_group = 3, age_space_group = "set",
    maturity_sampling = "none"),
  Sentinel = sim_survey(q = gillnet,
    custom_sets = sim_sets(pop, subset_cells = Sentinel, set_den = 4 / 1000),
    lengths_cap = 100, ages_cap = 3, age_length_group = 3, age_space_group = "set",
    maturity_sampling = "none")
) |>
  run_strat(alk_scale = "set") |>
  strat_error(scope = "survey_frame")

# Convert SimSurvey's estimates from fish to thousand fish.
catch_at_age <- fishery$catch_estimates$catch_at_age
catch <- catch_at_age[,
  .(year, age, fleet = fishery, obs = n_hat / 1000, unit = "thousand fish")]
index <- surveys$age_strat[!is.na(age),
  .(year, age, survey, obs = total / 1000, samp_time = 0, unit = "thousand fish")]
monitoring <- fishery$catch_diagnostics[, .(fleet = fishery, year,
  reported_total = reported_n, n_operations, n_aged, status = age_status)]
monitoring <- merge(monitoring, catch[, .(estimated_total = sum(obs) * 1000), by = .(fleet, year)],
  by = c("fleet", "year"))
survey_precision <- surveys$age_strat[!is.na(age),
  .(year, age, survey, sets, total_se = sd * sampling_units / 1000, estimate_status)]
survey_monitoring <- rbindlist(lapply(surveys$survey_programs, "[[", "samp_totals"), idcol = "survey")

# Use known stock biology. Fill ages lacking a monitored catch weight with
# stock weight so tinyAM has complete weights through age 20.
stock <- stock_biology(pop)$at_age[, .(year, age, weight = mean_weight_kg, maturity = prop_mature)]
catch_biology <- catch_at_age[, .(year, age, fleet = fishery, mean_weight_kg, weight_status)]
weights <- rbind(
  stock[, .(year, age, obs = weight, type = "stock", fleet = NA_character_, unit = "kg")],
  stock[catch_biology, on = .(year, age),
    .(year, age, obs = fcoalesce(mean_weight_kg, weight), type = "catch", fleet, unit = "kg")]
)
maturity <- stock[, .(year, age, obs = maturity, spawn_time = 0)]
sim_cod_obs <- tinyAM::as_tam_data(lapply(list(catch = catch, index = index,
  weight = weights, maturity = maturity), as.data.frame))
tinyAM::check_tam_data(sim_cod_obs, quiet = TRUE)

# Keep realized states and monitoring diagnostics separate from the fitting data.
population <- stock[, .(year, age, weight, maturity, N = as.numeric(pop$N) / 1000,
  Z = as.numeric(pop$Z), B = as.numeric(pop$B) / 1000, SSB = as.numeric(pop$S) / 1000)]
realized <- fishery$catch_error$catch_at_age[, .(year, age, fishery, n_true, mean_weight_true_kg)]
fleet <- merge(fishery$reference, realized, by = c("fishery", "year", "age"))
fleet[stock, on = .(year, age), weight := fcoalesce(mean_weight_true_kg, i.weight)]
fleet[, `:=`(yield = n_true * weight / 1000, n_true = NULL, mean_weight_true_kg = NULL)]
setnames(fleet, "fishery", "fleet")
counts <- c("N", "catch_reference", "catch_realized", "total_deaths", "death_budget")
fleet[, (counts) := lapply(.SD, function(x) x / 1000), .SDcols = counts]
rates <- fleet[, .(F = sum(F_realized), F_reference = sum(F_reference)), by = .(year, age)]
population <- merge(population, rates, by = c("year", "age"))
population[, `:=`(M = Z - F, M_reference = Z - F_reference)]
annual <- merge(population[, .(N = sum(N), B = sum(B), SSB = sum(SSB),
  recruitment = N[age == min(age)]), by = year],
  fleet[, .(yield = sum(yield)), by = year], by = "year")

curves <- list(Gillnet = gillnet, Longline = longline, Trawl = trawl, RV = rv, Sentinel = gillnet)
survey_truth <- surveys$age_strat_error[!is.na(age),
  .(year, age, survey, N_region = N_frame / 1000, N_total = N_stock / 1000,
    index_expected = I / 1000)]
survey_truth[, selection := curves[[.BY$survey]](age), by = survey]
survey_truth[, availability := fifelse(N_total > 0, N_region / N_total, NA_real_)]
selectivity <- rbindlist(lapply(names(curves), function(name) {
  age <- seq(1, 20, by = 0.1)
  data.table(component = name, age, selectivity = curves[[name]](age))
}))
footprint$area_km2 <- abs(prod(stars::st_res(grid)))

tail <- population[age == max(age)]
tail[, N := N * exp(-Z)]
tail[, SSB := N * weight * maturity]
tail <- merge(tail[, .(year, N, SSB)], annual[, .(year, N_total = N, SSB_total = SSB)], by = "year")
tail[, `:=`(proportion = N / N_total, SSB_proportion = SSB / SSB_total)]
M_center <- population[age %in% 5:10 & is.finite(M), M]
fishery_annual <- copy(fishery$annual)
fishery_annual[, F_full_unreduced := F_full_unreduced[[.BY$fishery]][match(year, years)], by = fishery]
checks <- list(terminal_survivors = as.data.frame(tail),
  shortfalls = as.data.frame(fishery$annual[shortfall > 0, .(fishery, year, target, n, shortfall, status)]),
  unsupported_catch = as.data.frame(monitoring[!status %in% c("estimated", "zero_reported_catch")]),
  realized_M = list(ages = 5:10, mean = mean(M_center), sd = sd(M_center),
    min = min(M_center), max = max(M_center)),
  realized_M_by_period = as.data.frame(population[age %in% 5:10 & is.finite(M),
    .(mean = mean(M), sd = sd(M), min = min(M), max = max(M)),
    by = .(period = ifelse(year < trawl_reduction_year, "Before", "After"))]),
  realized_M_at_age = as.data.frame(population[is.finite(M),
    .(mean = mean(M), sd = sd(M), min = min(M), max = max(M)), by = age]))
provenance <- list(simsurvey_revision = revision, seed = 17L, fishing_seed = 1701L,
  fishing_log_step_sd = 0.10, trawl_reduction_year = trawl_reduction_year,
  reference_M = reference_M,
  rng_kind = RNGkind(),
  generator = unname(tools::md5sum("data-raw/sim_cod_obs.R")),
  session = capture.output(utils::sessionInfo()),
  sources = c(population = "https://doi.org/10.1371/journal.pone.0232822",
    data_review = "https://publications.gc.ca/collections/collection_2024/mpo-dfo/Fs97-13-1381-eng.pdf",
    sentinel = "https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-056-eng.pdf"))
sim_cod_truth <- lapply(list(population = population, annual = annual, fleet = fleet,
  surveys = survey_truth, footprint = footprint, selectivity = selectivity,
  fishery_annual = fishery_annual, monitoring = monitoring, catch_biology = catch_biology,
  survey_precision = survey_precision, survey_monitoring = survey_monitoring), as.data.frame)
sim_cod_truth$checks <- checks
sim_cod_truth$provenance <- provenance

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(sim_cod_obs = sim_cod_obs, sim_cod_truth = sim_cod_truth), output, compress = "xz")
message("Staged data at ", output)
message(sprintf("Realized M at ages 5–10: mean %.3f, SD %.3f, range %.3f–%.3f per year",
  mean(M_center), sd(M_center), min(M_center), max(M_center)))
print(checks$realized_M_at_age)
print(checks$realized_M_by_period)
print(checks$shortfalls)
stopifnot(
  abs(mean(M_center) - 0.30) <= 0.03,
  all(abs(checks$realized_M_by_period$mean - 0.30) <= 0.03),
  annual$SSB[which.max(annual$year)] < annual$SSB[which.min(annual$year)],
  all(fishery_annual$F_full_reference == fishery_annual$F_full_unreduced *
    ifelse(fishery_annual$fishery == "Trawl" & fishery_annual$year >= trawl_reduction_year, 0.2, 1))
)
