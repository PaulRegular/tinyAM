# Whole-area monitoring experiment. SimSurvey is a study-only dependency.
pld_simsurvey_revision <- "29a1b03350c274f88dac108fb95f9446eac79233"
pld_years <- 1:30
pld_ages <- 1:20
pld_m <- 0.30
pld_pulse <- 14:16
pld_grid <- 40L
pld_selection <- function(age) exp(-0.5 * ((age - 7.7) / 2.5)^2)
pld_q <- function(age) stats::plogis(2 * (age - 3))
pld_older_window <- 16:23

pld_simsurvey <- function(id) {
  seeds <- c(fishing = 61000L, population = 62000L, spatial = 63000L,
    fishery = 64000L, monitoring = 75000L, survey = 86000L) + 100L * id
  baseline <- withr::with_seed(seeds[["fishing"]],
    0.12 * exp(c(0, cumsum(stats::rnorm(29L, sd = 0.10)))))
  target <- baseline
  target[pld_pulse] <- 3 * target[pld_pulse]
  history <- withr::with_seed(seeds[["population"]], SimSurvey::sim_abundance(
    years = pld_years, Z = function(ages, years) pld_m + outer(pld_selection(ages), target)))
  exact_last <- history$N[20L, 1L]
  for (j in pld_years[-1L]) {
    history$N[20L, j] <- history$N[20L, j] + history$N[20L, j - 1L] *
      exp(-history$Z[20L, j - 1L])
  }
  stopifnot(history$N[20L, 1L] == exact_last, all(history$N[20L, ] > 0))
  expected <- round(history$N)
  grid <- SimSurvey::make_grid(res = rep(280 / pld_grid, 2L))
  pop <- withr::with_seed(seeds[["spatial"]], SimSurvey::sim_distribution(history, grid = grid))
  pop$sp_N[, N := .oe_balance(N, expected[.BY$age, .BY$year]), by = .(year, age)]
  pop <- SimSurvey::round_sim(pop)
  stopifnot(isTRUE(all.equal(pop$N, expected, tolerance = 0)))
  biology <- as.data.frame(SimSurvey::stock_biology(pop)$at_age)
  abundance <- as.numeric(pop$N) / 1000
  population <- data.frame(biology[c("year", "age")], N = abundance,
    Z = as.numeric(pop$Z), ssb = as.numeric(pop$S) / 1000,
    index_expected = abundance * pld_q(biology$age))
  weight <- data.frame(biology[c("year", "age")], obs = biology$mean_weight_kg, unit = "kg")
  # tinyAM multiplies maturity by mean individual weight. Use E[wm]/E[w]
  # so SSB equals SimSurvey's length-integrated mature weight.
  maturity <- data.frame(biology[c("year", "age")],
    obs = biology$mean_mature_weight_kg / biology$mean_weight_kg, spawn_time = 0)
  simulator <- pld_fishery_simulator(pld_fast_fishery_engine)
  fishery <- SimSurvey::sim_fishery(target = target[1L], target_type = "F_full",
    selectivity = pld_selection,
    catch_model = SimSurvey::sim_capture(q = 0.02, capacity_per_gear = 150),
    gear_amount = 1200L, duration = 24, gear_unit = "net", time_unit = "hour",
    sampling_coverage = 0.50, lengths_cap = 150, ages_cap = 150,
    age_sampling = "random", reporting_cv = 0.05, maturity_sampling = "none", light = FALSE)
  draws <- lapply(pld_years, function(year) {
    snapshot <- oe_snapshot(pop, year)
    spec <- fishery
    spec$target <- target[year]
    fished <- withr::with_seed(seeds[["fishery"]] + year,
      simulator(snapshot, Gillnet = spec))
    fished <- withr::with_seed(seeds[["monitoring"]] + year, pld_monitor_fishery(fished, spec))
    catch <- fished$catch_estimates$catch_at_age
    catch <- as.data.frame(catch[, .(year, age, obs = n_hat / 1000, unit = "thousand fish")])
    catch_diagnostics <- as.data.frame(fished$catch_diagnostics)
    annual <- as.data.frame(fished$annual[, .(year, target, n, shortfall, n_operations, status)])
    f <- as.data.frame(fished$reference[, .(year, age, F_reference, F_realized, catch_reference,
      catch_realized, death_budget)])
    survey <- withr::with_seed(seeds[["survey"]] + year, pld_monitor_survey(snapshot))
    index <- survey$index
    stopifnot(identical(as.integer(catch$age), pld_ages),
      identical(as.integer(index$age), pld_ages), all(catch$year == year), all(index$year == year),
      all(f$catch_realized <= f$death_budget),
      abs(sum(f$catch_realized) - annual$n) < 1e-8)
    list(catch = catch, index = index, catch_diagnostics = catch_diagnostics,
      survey_diagnostics = survey$diagnostics, f = f, annual = annual)
  })
  gather <- function(name) as.data.frame(data.table::rbindlist(lapply(draws, "[[", name)))
  reference <- merge(population, gather("f"), by = c("year", "age"), sort = FALSE)
  reference$M_realized <- reference$Z - reference$F_realized
  stopifnot(nrow(reference) == 600L,
    max(abs(reference$F_reference - target[reference$year] * pld_selection(reference$age))) < 1e-12,
    max(abs(reference$Z - reference$F_reference - pld_m)) < 1e-12)
  list(id = id, seeds = seeds, data = tinyAM::as_tam_data(list(catch = gather("catch"),
    index = gather("index"), weight = weight, maturity = maturity)), reference = reference,
    schedule = data.frame(year = pld_years, baseline = baseline, F_full = target),
    fishery_diagnostics = gather("annual"), catch_diagnostics = gather("catch_diagnostics"),
    survey_diagnostics = gather("survey_diagnostics"))
}

pld_monitor_fishery <- function(fished, spec) {
  configure <- getFromNamespace(".fishery_controls", "SimSurvey")
  monitor <- getFromNamespace(".monitor_fisheries", "SimSurvey")
  weights <- setNames(fished$sim_weight(fished$lengths), fished$lengths)
  configuration <- list(Gillnet = configure(spec, fished, "Gillnet"))
  fished$fishery <- monitor(fished$fishery, configuration, fished$years, fished, weights)
  fished$fishery$configurations <- list(Gillnet = spec)
  SimSurvey::run_catch(fished)$fishery
}

pld_monitor_survey <- function(snapshot) {
  trawl_dim <- c(1.5, 0.0018 / 1.5)
  sets <- SimSurvey::sim_sets(snapshot, set_den = 2 / 1000, trawl_dim = trawl_dim)
  stopifnot(max(sets$cell_sets) == 1, all(sets$strat_sets >= 2))
  survey <- SimSurvey::sim_survey(snapshot, custom_sets = sets, q = SimSurvey::sim_logistic(),
    lengths_cap = 100, ages_cap = 3, age_length_group = 3, age_space_group = "set",
    maturity_sampling = "none", light = FALSE)
  index <- as.data.frame(SimSurvey::run_strat(survey, alk_scale = "set")$age_strat)
  index <- index[!is.na(index$age), , drop = FALSE]
  # A positive assigned count survives a partial allocation, while a zero
  # assignment with unresolved age information is missing.
  incomplete <- index$estimate_status %in% "partial_age" & !is.na(index$total) & index$total == 0
  index$total[incomplete] <- NA_real_
  list(index = data.frame(index[c("year", "age")], survey = "RV", samp_time = 0,
    obs = index$total / 1000, unit = "thousand fish"),
    diagnostics = index[c("year", "age", "sets", "estimate_status")])
}

pld_mortality_audit <- function(simulation) {
  use <- subset(simulation$reference, age %in% 5:10 & N > 0)
  departure <- use$M_realized - pld_m
  by_year <- vapply(split(use, use$year), function(d) {
    stats::weighted.mean(d$M_realized - pld_m, d$N)
  }, numeric(1L))
  mean_absolute <- stats::weighted.mean(abs(departure), use$N)
  worst_annual <- max(abs(by_year))
  data.frame(id = simulation$id, mean_absolute = mean_absolute, worst_annual = worst_annual,
    shortfalls = sum(simulation$fishery_diagnostics$shortfall > 0),
    adequate = is.finite(mean_absolute) && mean_absolute <= 0.03 && worst_annual <= 0.06)
}


# GPL-3. SimSurvey 29a1b03 allocation law, batched across cells for speed.
pld_allocate_lengths <- function(size, probability) {
  allocation <- matrix(0, length(size), length(probability))
  for (j in seq_along(probability)[-length(probability)]) {
    p <- probability[j] / sum(probability[j:length(probability)])
    allocation[, j] <- stats::rbinom(length(size), size, p)
    size <- size - allocation[, j]
  }
  allocation[, length(probability)] <- size
  allocation
}

pld_fishery_stock <- function(sim, year, key) {
  rows <- as.data.frame(sim$sp_N[sim$sp_N$year == year, c("cell", "age", "N")])
  by_age <- split(rows, rows$age)
  stock <- data.table::rbindlist(lapply(seq_along(sim$ages), function(i) {
    d <- by_age[[as.character(sim$ages[i])]]
    keep <- which(key[, i] > 0)
    counts <- pld_allocate_lengths(d$N, key[keep, i])
    stopifnot(all(rowSums(counts) == d$N))
    result <- data.table::data.table(cell = rep(d$cell, each = length(keep)),
      age = sim$ages[i], length = rep(as.numeric(rownames(key))[keep], nrow(d)),
      N = as.numeric(t(counts)))
    result[result$N > 0]
  }))
  data.table::setorderv(stock, c("cell", "age", "length"))
  stock
}

pld_fishery_simulator <- function(adapt_engine = identity) {
  simulate <- SimSurvey::sim_fisheries
  engine <- getFromNamespace(".run_fisheries", "SimSurvey")
  if (!".fishery_stock" %in% all.names(body(engine))) {
    stop("The pinned SimSurvey fishery engine changed. Audit the batched length allocation.")
  }
  engine <- adapt_engine(engine)
  private <- new.env(parent = environment(simulate))
  private$.fishery_stock <- pld_fishery_stock
  environment(engine) <- private
  private$.run_fisheries <- engine
  environment(simulate) <- private
  simulate
}

pld_check_allocate_lengths <- function() {
  withr::with_seed(77L, {
    size <- rep.int(80L, 10000L)
    probability <- c(0.05, 0.2, 0.35, 0.4)
    x <- pld_allocate_lengths(size, probability)
    expected <- size[1L] * probability
    theoretical <- size[1L] * (diag(probability) - tcrossprod(probability))
    # Fixed-seed check of exact row totals and approximate first two moments
    # against the multinomial law, including its negative age dependence.
    stopifnot(all(rowSums(x) == size), max(abs(colMeans(x) - expected)) < 0.2,
      max(abs(stats::cov(x) - theoretical)) < 1.2)
  })
}
