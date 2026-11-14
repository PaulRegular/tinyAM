# Simulation and fit helpers for simulate_pld.R. SimSurvey is an analysis dependency.
pld_simsurvey_revision <- "29a1b03350c274f88dac108fb95f9446eac79233"
pld_effort <- data.frame(effort = c("low", "standard", "high"),
  sets_per_1000 = c(1, 2, 4), operation_coverage = c(0.15, 0.30, 0.60))
pld_y50 <- c(catch = 1.5, index = 100)
pld_sd <- c(catch = 0.25, index = 0.35)

pld_m <- function(year) stats::approx(c(1, 5, 9, 12, 16, 20),
  c(0.25, 0.25, 0.60, 0.60, 0.25, 0.25), xout = year)$y

pld_f <- function(age) 0.06 * exp(-0.5 * ((age - 7.7) / 2.5)^2)
pld_q <- function(age) stats::plogis(2 * (age - 3))

pld_seeds <- function(id) {
  c(population = 21000L + id, spatial = 22000L + id, fishery = 23000L + 1000L * id,
    catch_monitoring = 55000L + 1000L * id, rv = 85000L + 1000L * id,
    exact = 26000L + id, multinomial = 27000L + id)
}

pld_demography <- function(seed, ages = 1:20) {
  withr::with_seed(seed, SimSurvey::sim_abundance(ages = ages,
    Z = function(ages, years) outer(rep(1, length(ages)), pld_m(years)) +
      outer(pld_f(ages), rep(1, length(years)))))
}

pld_biology <- function(population) {
  stock <- SimSurvey::stock_biology(population)$at_age
  # SSB = N * weight * maturity in tinyAM. Conditional mature weight and
  # mean maturity are not independent in SimSurvey's default length model.
  stock$effective_maturity <- stock$mean_mature_weight_kg / stock$mean_weight_kg
  list(weight = data.frame(year = stock$year, age = stock$age,
      obs = stock$mean_weight_kg, unit = "kg"),
    maturity = data.frame(year = stock$year, age = stock$age,
      obs = stock$effective_maturity, spawn_time = 0),
    maturity_at_age = as.data.frame(stock))
}

pld_annual_truth <- function(years, ages, N, biology) {
  young <- which.min(ages)
  weights <- biology$maturity_at_age$mean_mature_weight_kg
  annual <- data.frame(year = years, ssb = colSums(N * matrix(weights, nrow(N))) / 1000,
    recruitment = N[young, ] / 1000, abundance = colSums(N) / 1000, M_bar = pld_m(years))
  annual
}

pld_count <- function(rows, value, kind, age = rows$age) {
  cols <- list(year = rows$year, age = age)
  if (kind == "index") cols <- c(cols, list(survey = rows$survey, samp_time = rows$samp_time))
  else if ("fleet" %in% names(rows)) cols$fleet <- rows$fleet
  data.frame(cols, obs = value, unit = "thousand fish")
}

pld_pool <- function(rows, plus_age) {
  # All native constituent values are required for an aggregate observation.
  rows$age <- pmin(rows$age, plus_age)
  columns <- setdiff(names(rows), "obs")
  stats::aggregate(rows["obs"], rows[columns], sum)
}

pld_data <- function(catch, index, biology, plus_age) {
  obs <- list(catch = pld_pool(catch, plus_age), index = pld_pool(index, plus_age),
    weight = biology$weight, maturity = biology$maturity)
  # Catch-grid validation requires existing biological-age rows. The older
  # catches have already been counted in the terminal bin.
  pad <- catch[catch$age > plus_age, , drop = FALSE]
  if (nrow(pad)) {
    pad$obs <- NA_real_
    obs$catch <- rbind(obs$catch, pad)
  }
  tinyAM::as_tam_data(obs)
}

pld_sampling_summary <- function(obs, kind) {
  inputs <- if (kind == "catch") obs$native_catch else obs$native_index
  do.call(rbind, lapply(c(14L, 20L), function(age) {
    tab <- pld_pool(inputs, age)
    data.frame(component = kind, plus_age = age, n_zero = sum(tab$obs == 0, na.rm = TRUE),
      n_positive = sum(tab$obs > 0, na.rm = TRUE), n_missing = sum(is.na(tab$obs)))
  }))
}

pld_effective_mortality <- function(reference) {
  reference <- as.data.frame(reference)
  transform(reference, M_reference = pld_m(year), M_realized = Z - F_realized)
}

# SimSurvey at 29a1b03, R/fishery_engine.R:.fishery_stock and
# R/sim_traits.R:.multinomial_counts, GPL-3. Use the same sequential binomial
# law while batching independent cells by age. This changes the RNG order.
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

pld_fishery_simulator <- function() {
  simulate <- SimSurvey::sim_fisheries
  engine <- getFromNamespace(".run_fisheries", "SimSurvey")
  if (!".fishery_stock" %in% all.names(body(engine))) {
    stop("The pinned SimSurvey fishery engine changed. Audit the batched length allocation.")
  }
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

pld_monitor_rv_year <- function(snapshot, set_density, seed) {
  withr::with_seed(seed, {
    sets <- SimSurvey::sim_sets(snapshot, set_den = set_density / 1000)
    # Sampling is without replacement and survey-wide stratum coverage is
    # needed to interpret a zero as an encounter outcome.
    stopifnot(max(sets$cell_sets) == 1, all(sets$strat_sets <= sets$strat_cells),
      all(sets$strat_sets >= 2))
    survey <- SimSurvey::sim_survey(snapshot, custom_sets = sets, q = SimSurvey::sim_logistic(),
      lengths_cap = 100, ages_cap = 3, age_length_group = 3,
      age_space_group = "set", maturity_sampling = "none", light = FALSE)
    estimated <- SimSurvey::run_strat(survey, alk_scale = "set")
    # Only full_setdet and the public survey roster carry age. setdet is
    # summarized by set, while full_setdet needs set metadata for this year.
    captured <- as.data.frame(survey$full_setdet)
    captured$year <- rep(snapshot$years, nrow(captured))
    list(index = estimated$age_strat, coverage = estimated$survey_coverage,
      caught = stats::aggregate(captured["n"], captured[c("year", "age")], sum),
      sampling = survey$samp_totals)
  })
}

pld_simsurvey_inputs <- function(id, grid_n, checkpoint_dir = NULL) {
  seeds <- pld_seeds(id)
  history <- pld_demography(seeds[["population"]])
  grid <- SimSurvey::make_grid(res = rep(280 / grid_n, 2L))
  cells <- as.data.frame(grid)
  frame_cells <- with(cells, table(strat))
  cell_area <- abs(prod(stars::st_res(grid)))
  high <- pmax(2, round(as.numeric(frame_cells) * cell_area * max(pld_effort$sets_per_1000) / 1000))
  if (any(high > frame_cells)) stop("The grid cannot support high-effort sampling without replacement.")
  # Reuse the prior study's balanced rounding; avoid adding spatially varying
  # changes to our prespecified cohort totals through independent rounding.
  pop <- withr::with_seed(seeds[["spatial"]],
    SimSurvey::sim_distribution(history, grid = grid))
  totals <- round(history$N)
  pop$sp_N[, N := .oe_balance(N, totals[.BY$age, .BY$year]), by = .(year, age)]
  pop <- SimSurvey::round_sim(pop)
  stopifnot(isTRUE(all.equal(pop$N, totals, tolerance = 0)))
  biology <- pld_biology(pop)
  cat("  History", id, "grid and biology prepared. Allocating fishery removals.\n")
  flush.console()

  # A single set of fishery operations supplies all effort comparisons.
  spec <- SimSurvey::sim_fishery(target = 0.06, target_type = "F_full",
    selectivity = function(age) pld_f(age) / 0.06,
    catch_model = SimSurvey::sim_capture(q = 0.02, capacity_per_gear = 150),
    gear_amount = 1200, duration = 24, gear_unit = "net", time_unit = "hour",
    sampling_coverage = 0, reporting_cv = 0.05, lengths_cap = 150,
    ages_cap = 150, age_sampling = "random", maturity_sampling = "none", light = FALSE)
  simulate_fisheries <- pld_fishery_simulator()
  within_draw <- function(stage, year, calculate) {
    if (is.null(checkpoint_dir)) return(calculate())
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(checkpoint_dir, sprintf("%s-year%02d.rds", stage, year))
    if (file.exists(path)) return(readRDS(path))
    value <- calculate()
    temporary <- paste0(path, ".partial")
    saveRDS(value, temporary)
    if (!file.rename(temporary, path)) stop("Could not save checkpoint ", path)
    value
  }
  fishery_years <- lapply(pop$years, function(year) within_draw("fishery", year, function() {
    snapshot <- oe_snapshot(pop, year)
    annual <- withr::with_seed(seeds[["fishery"]] + year,
      simulate_fisheries(snapshot, Gillnet = spec))$fishery
    cat("    Fishery year", year, "operations:", nrow(annual$setdet), "\n")
    flush.console()
    # run_catch() needs the realized catches and monitoring records. Keep
    # only the total deaths by year and age already recorded in reference.
    annual$mortality_ledger <- NULL
    annual
  }))
  # The simulator numbers operations within each call. Offset the IDs to keep
  # a unique join key when monitoring the combined history.
  offset <- 0L
  for (i in seq_along(fishery_years)) {
    fishery <- fishery_years[[i]]
    fishery$setdet$operation <- fishery$setdet$operation + offset
    fishery$joint_catch$operation <- fishery$joint_catch$operation + offset
    offset <- max(c(offset, fishery$setdet$operation))
    fishery_years[[i]] <- fishery
  }
  gather <- function(name) data.table::rbindlist(lapply(fishery_years, `[[`, name))
  fields <- c("setdet", "joint_catch", "annual", "reference")
  truth <- setNames(lapply(fields, gather), fields)
  # The engine returns already monitored summaries. The monitoring adapter
  # expects the allocation ledger's annual columns before those joins.
  monitoring_fields <- c("kg", "reported_n", "reported_kg", "n_measured", "n_aged",
    "n_maturity", "cpue", "reported_cpue", "biomass_cpue", "reported_biomass_cpue")
  truth$annual <- truth$annual[, setdiff(names(truth$annual), monitoring_fields), with = FALSE]
  truth$mortality_ledger <- data.table::data.table()
  truth$has_weight <- TRUE
  rm(fishery_years)
  cat("  History", id, "fishery removals allocated. Starting monitoring.\n")
  flush.console()
  reference <- pld_effective_mortality(truth$reference)
  stopifnot(isTRUE(all.equal(reference$F_reference, pld_f(reference$age))))
  # A pinned internal adapter reapplies SimSurvey monitoring to the same
  # operations. The exported sim_fisheries() would draw new removals.
  fishery_controls <- getFromNamespace(".fishery_controls", "SimSurvey")
  monitor <- getFromNamespace(".monitor_fisheries", "SimSurvey")
  make_cfg <- function(probability) {
    spec$sampling_coverage <- probability
    list(Gillnet = fishery_controls(spec, pop, "Gillnet"))
  }
  stock_lengths <- pop$lengths
  weight_at_length <- setNames(pop$sim_weight(stock_lengths), stock_lengths)
  series <- lapply(seq_len(nrow(pld_effort)), function(i) {
    current <- pld_effort[i, ]
    cat("  History", id, current$effort, "monitoring.\n")
    flush.console()
    config <- make_cfg(current$operation_coverage)
    monitored <- withr::with_seed(seeds[["catch_monitoring"]] + i,
      monitor(data.table::copy(truth), config,
        pop$years, pop, weight_at_length))
    spec$sampling_coverage <- current$operation_coverage
    monitored$configurations <- list(Gillnet = spec)
    fished <- pop
    fished$fishery <- monitored
    fished$replicates <- 1L
    catch_est <- SimSurvey::run_catch(fished)$fishery
    # Do not append NA results to the original capture table and re-expand.
    catch_rows <- catch_est$catch_estimates$catch_at_age
    catch <- pld_count(catch_rows, catch_rows$n_hat / 1000, "catch")
    stopifnot(nrow(catch) == length(pop$ages) * length(pop$years))

    rv <- lapply(pop$years, function(year) within_draw(paste0("rv-", current$effort), year, function() {
      pld_monitor_rv_year(oe_snapshot(pop, year), current$sets_per_1000,
        seeds[["rv"]] + 30L * i + year)
    }))
    rv_table <- function(name) as.data.frame(data.table::rbindlist(lapply(rv, `[[`, name)))
    rv_rows <- rv_table("index")
    rv_rows <- rv_rows[!is.na(rv_rows$age), , drop = FALSE]
    # Partial estimates can retain a positive assigned catch for some ages.
    # Zero assigned catch with unresolved age support is not a known zero.
    uncertain_zero <- rv_rows$estimate_status %in% "partial_age" &
      !is.na(rv_rows$total) & rv_rows$total == 0
    rv_rows$total[uncertain_zero] <- NA_real_
    rv_rows$survey <- "RV"
    rv_rows$samp_time <- 0
    index <- pld_count(rv_rows, rv_rows$total / 1000, "index")
    stopifnot(nrow(index) == length(pop$ages) * length(pop$years))
    expected <- data.frame(year = rep(pop$years, each = length(pop$ages)),
      age = rep(pop$ages, times = length(pop$years)), mu = as.numeric(pop$N * pld_q(pop$ages)) / 1000)
    raw <- list(native_catch = catch, native_index = index)
    list(data = lapply(c(14L, 20L), function(a) pld_data(catch, index, biology, a)),
      observations = raw, effort = current, expected_index = expected,
      zeros = rbind(pld_sampling_summary(raw, "catch"), pld_sampling_summary(raw, "index")),
      catch_diagnostics = as.data.frame(catch_est$catch_diagnostics),
      survey_coverage = rv_table("coverage"), survey_caught = rv_table("caught"),
      survey_sampling = rv_table("sampling"),
      survey_status = rv_rows[c("year", "age", "estimate_status")])
  })
  oldest <- nrow(pop$N)
  survivor <- as.numeric(pop$N[oldest, ] * exp(-pop$Z[oldest, ]))
  list(id = id, seeds = seeds, series = setNames(series, pld_effort$effort),
    annual = pld_annual_truth(pop$years, pop$ages, pop$N, biology),
    reference = reference, annual_fishery = as.data.frame(truth$annual),
    oldest = data.frame(year = pop$years, survivor = survivor,
      proportion = survivor / colSums(pop$N), oldest_age = pop$ages[oldest]))
}

pld_rsample <- function(mu, h50, sd) {
  prob <- -expm1(-log(2) * mu / h50)
  draws <- rep(0, length(mu))
  positive <- stats::runif(length(mu)) < prob
  draws[positive] <- stats::rlnorm(sum(positive),
    log(mu[positive] / prob[positive]) - sd^2 / 2, sd)
  draws
}

pld_rmultinomial <- function(mu, effort, h50) {
  total <- colSums(mu)
  n_sample <- pmax(1, round(log(2) * effort * total / h50))
  p <- sweep(mu, 2, total, "/")
  counts <- vapply(seq_len(ncol(mu)), function(y) {
    drop(stats::rmultinom(1, n_sample[y], p[, y]))
  }, integer(nrow(mu)))
  as.numeric(sweep(counts, 2, total / n_sample, "*"))
}

pld_controlled_inputs <- function(id) {
  seeds <- pld_seeds(id)
  pop <- pld_demography(seeds[["population"]], 1:14)
  biology <- pld_biology(pop)
  # Generate exactly the modeled plus-group transition. No hidden older
  # population or observation pooling is needed in the controlled cases.
  for (y in 2:ncol(pop$N)) {
    old <- nrow(pop$N)
    pop$N[old, y] <- pop$N[old, y] + pop$N[old, y - 1L] * exp(-pop$Z[old, y - 1L])
  }
  catches <- pop$N * pld_f(pop$ages) * (-expm1(-pop$Z)) / pop$Z / 1000
  indexes <- pop$N * pld_q(pop$ages) / 1000
  keys <- expand.grid(age = pop$ages, year = pop$years)
  keys$survey <- "RV"
  keys$samp_time <- 0
  # Effort changes repeat through every mortality phase.
  effort <- rep(c(0.5, 1, 2, 1), length.out = length(pop$years))
  make_data <- function(catch, index) pld_data(
    pld_count(keys, catch, "catch"), pld_count(keys, index, "index"), biology, 14L)
  exact <- withr::with_seed(seeds[["exact"]],
    make_data(pld_rsample(catches, pld_y50[["catch"]], pld_sd[["catch"]]),
      pld_rsample(indexes, pld_y50[["index"]], pld_sd[["index"]])))
  multinomial <- withr::with_seed(seeds[["multinomial"]],
    make_data(pld_rmultinomial(catches, effort, pld_y50[["catch"]]),
      pld_rmultinomial(indexes, effort, pld_y50[["index"]])))
  list(id = id, seeds = seeds, exact = exact, multinomial = multinomial, effort = effort,
    annual = pld_annual_truth(pop$years, pop$ages, pop$N, biology),
    prediction = data.frame(keys, catch = as.numeric(catches), index = as.numeric(indexes)))
}
