# Run from the tinyAM root with the SimSurvey revision documented in README.md.
# Only the compact prepared data are saved. Redirect console output to keep a run log.
library(SimSurvey)
library(data.table)

revision <- "42e5def3819a5eec443deaabd5d268c7f126438d"
installed_revision <- utils::packageDescription("SimSurvey", fields = "RemoteSha")
if (is.na(installed_revision)) installed_revision <- Sys.getenv("TINYAM_SIMSURVEY_REVISION")
stopifnot(identical(installed_revision, revision))
output <- Sys.getenv("TINYAM_SIM_COD_OUTPUT",
  "data-raw/cache/sim_cod/f-full-age/prepared-data.rds")
if (file.exists(output)) stop("Choose a fresh TINYAM_SIM_COD_OUTPUT path.")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(438)

# Three fishing rates at full selectivity. Reuse the curves when generating Z.
F_full <- c(Gillnet = 0.06, Longline = 0.05, Trawl = 0.07)
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
pop <- sim_abundance(
  Z = function(ages, years) {
    M <- sim_Z(log_mean = log(0.3), log_sd = 0.04)(ages, years)
    M + F_full[["Gillnet"]] * gillnet(ages) +
      F_full[["Longline"]] * longline(ages) + F_full[["Trawl"]] * trawl(ages)
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
)$fishery[c("annual", "setdet", "samp", "catch_at_age", "reference")]

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
)

# Raise sampled age proportions by reported operation size and inverse operation
# sampling probability. Calibrate each fleet-year to its annual reported total.
aged <- fishery$samp[aged == TRUE, .(n_age = .N), by = .(operation, age)]
aged[, proportion := n_age / sum(n_age), by = operation]
raised <- merge(aged, fishery$setdet[, .(operation, fishery, year, reported_n, sampling_probability)],
  by = "operation")
raised <- raised[, .(expanded = sum(reported_n * proportion / sampling_probability)),
  by = .(fleet = fishery, year, age)]
monitoring <- fishery$annual[, .(fleet = fishery, year, reported_total = reported_n,
  n_operations, n_aged)]
catch <- merge(CJ(fleet = names(F_full), year = pop$years, age = pop$ages),
  raised, by = c("fleet", "year", "age"), all.x = TRUE)
catch <- merge(catch, monitoring[, .(fleet, year, reported_total)], by = c("fleet", "year"))
catch[is.na(expanded), expanded := 0]
catch[, obs := if (reported_total[1] == 0) 0 else if (sum(expanded) > 0) {
  reported_total * expanded / sum(expanded) / 1000
} else NA_real_, by = .(fleet, year)]
monitoring <- merge(monitoring, catch[, .(estimated_total = sum(obs) * 1000), by = .(fleet, year)],
  by = c("fleet", "year"))
monitoring[, status := fifelse(reported_total == 0, "zero_total",
  fifelse(is.na(estimated_total), "no_samples", "ok"))]
catch <- catch[, .(year, age, fleet, obs, unit = "thousand fish")]

# Keep each survey's own footprint in the stratified analysis.
survey_tables <- lapply(names(surveys$survey_programs), function(survey) {
  selected <- select_survey(surveys, survey)
  analysis <- run_strat(selected, alk_scale = "set")
  # SimSurvey reports unassigned abundance at age NA separately from known ages.
  at_age <- analysis$age_strat[!is.na(age),
    .(year, age, obs = total / 1000, sets, total_se = sd * sampling_units / 1000)]
  at_age <- merge(CJ(year = pop$years, age = pop$ages), at_age,
    by = c("year", "age"), all.x = TRUE)
  list(index = at_age[, .(year, age, survey, obs, samp_time = 0, unit = "thousand fish")],
    precision = at_age[, .(year, age, survey, sets, total_se)],
    monitoring = copy(selected$samp_totals)[, survey := survey])
})
index <- rbindlist(lapply(survey_tables, "[[", "index"))

# Biology is supplied for every simulated age, including 14–20. In the fit,
# these ages retain separate biological weights and maturity inside age 14+.
lak <- pop$sim_length(age = pop$ages, length_age_key = TRUE)
weight_at_age <- colSums(lak * pop$sim_weight(as.numeric(rownames(lak))))
population <- as.data.table(expand.grid(age = pop$ages, year = pop$years))
population[, `:=`(N = as.numeric(pop$N) / 1000, Z = as.numeric(pop$Z),
  weight = weight_at_age[match(age, pop$ages)], maturity = pop$sim_maturity(age),
  B = as.numeric(pop$B) / 1000, SSB = as.numeric(pop$S) / 1000)]
fleet <- copy(fishery$reference)
setnames(fleet, "fishery", "fleet")
removals <- fishery$catch_at_age[, .(n = sum(n), kg = sum(kg)), by = .(fleet = fishery, year, age)]
fleet <- merge(fleet, removals, by = c("fleet", "year", "age"), all.x = TRUE)
fleet[, `:=`(weight = fifelse(!is.na(n) & n > 0, kg / n, weight_at_age[match(age, pop$ages)]),
  yield = fcoalesce(kg, 0) / 1000, n = NULL, kg = NULL)]
counts <- c("N", "catch_reference", "catch_realized", "total_deaths", "death_budget")
fleet[, (counts) := lapply(.SD, function(x) x / 1000), .SDcols = counts]

weights <- rbind(
  population[, .(year, age, obs = weight, type = "stock", fleet = NA_character_, unit = "kg")],
  fleet[, .(year, age, obs = weight, type = "catch", fleet, unit = "kg")]
)
maturity <- population[, .(year, age, obs = maturity, spawn_time = 0)]
sim_cod_obs <- tinyAM::as_tam_data(lapply(list(catch = catch, index = index,
  weight = weights, maturity = maturity), as.data.frame))
tinyAM::check_tam_data(sim_cod_obs, quiet = TRUE)

# Keep realized states and monitoring diagnostics separate from the fitting data.
rates <- fleet[, .(F = sum(F_realized), F_reference = sum(F_reference)), by = .(year, age)]
population <- merge(population, rates, by = c("year", "age"))
population[, `:=`(M = Z - F, M_reference = Z - F_reference)]
annual <- merge(population[, .(N = sum(N), B = sum(B), SSB = sum(SSB),
  recruitment = N[age == min(age)]), by = year],
  fleet[, .(yield = sum(yield)), by = year], by = "year")

curves <- list(Gillnet = gillnet, Longline = longline, Trawl = trawl, RV = rv, Sentinel = gillnet)
survey_truth <- rbindlist(lapply(names(surveys$survey_programs), function(survey) {
  cells <- footprint$cell[footprint[[survey]]]
  regional <- pop$sp_N[cell %in% cells, .(N_region = sum(N) / 1000), by = .(year, age)]
  regional <- merge(regional, population[, .(year, age, N_total = N)], by = c("year", "age"))
  regional[, `:=`(survey = survey, selection = curves[[survey]](age))]
  regional[, `:=`(index_expected = N_region * selection,
    availability = fifelse(N_total > 0, N_region / N_total, NA_real_))]
}))
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
checks <- list(terminal_survivors = as.data.frame(tail),
  shortfalls = as.data.frame(fishery$annual[shortfall > 0, .(fishery, year, target, n, shortfall, status)]),
  unsupported_catch = as.data.frame(monitoring[status == "no_samples"]),
  realized_M = list(ages = 5:10, mean = mean(M_center), sd = sd(M_center),
    min = min(M_center), max = max(M_center)),
  realized_M_at_age = as.data.frame(population[is.finite(M),
    .(mean = mean(M), sd = sd(M), min = min(M), max = max(M)), by = age]))
provenance <- list(simsurvey_revision = revision, seed = 438L, rng_kind = RNGkind(),
  generator = unname(tools::md5sum("data-raw/sim_cod_obs.R")),
  session = capture.output(utils::sessionInfo()),
  sources = c(population = "https://doi.org/10.1371/journal.pone.0232822",
    data_review = "https://publications.gc.ca/collections/collection_2024/mpo-dfo/Fs97-13-1381-eng.pdf",
    sentinel = "https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-056-eng.pdf"))
sim_cod_truth <- lapply(list(population = population, annual = annual, fleet = fleet,
  surveys = survey_truth, footprint = footprint, selectivity = selectivity,
  fishery_annual = fishery$annual, monitoring = monitoring,
  survey_precision = rbindlist(lapply(survey_tables, "[[", "precision")),
  survey_monitoring = rbindlist(lapply(survey_tables, "[[", "monitoring"))), as.data.frame)
sim_cod_truth$checks <- checks
sim_cod_truth$provenance <- provenance

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(sim_cod_obs = sim_cod_obs, sim_cod_truth = sim_cod_truth), output, compress = "xz")
message("Staged data at ", output)
message(sprintf("Realized M at ages 5–10: mean %.3f, SD %.3f, range %.3f–%.3f per year",
  mean(M_center), sd(M_center), min(M_center), max(M_center)))
print(checks$realized_M_at_age)
print(checks$shortfalls)
stopifnot(abs(mean(M_center) - 0.30) <= 0.03)
