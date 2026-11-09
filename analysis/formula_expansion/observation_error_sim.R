# SimSurvey is pinned for the study, independent of the packaged teaching data.
.observation_error_revision <- "29a1b03350c274f88dac108fb95f9446eac79233"
.oe_effort <- c(half = 0.5, standard = 1, double = 2)
.oe_plus_age <- 14L
.oe_tow_area <- 1.5 * 0.02

.observation_error_check_revision <- function() {
  installed <- utils::packageDescription("SimSurvey", fields = "RemoteSha")
  if (is.na(installed)) installed <- Sys.getenv("TINYAM_SIMSURVEY_REVISION")
  if (!identical(installed, .observation_error_revision)) {
    stop("This study requires SimSurvey at ", .observation_error_revision, ".")
  }
}

.observation_error_with_seed <- function(seed, code) {
  withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

.oe_gillnet_q <- function(age) {
  stats::plogis(1.4 * (age - 5)) * stats::plogis(0.8 * (10.5 - age))
}

.oe_rv_q <- function(age) stats::plogis(1.8 * (age - 2.5))

.oe_balance <- function(n, total) {
  count <- floor(n)
  add <- as.integer(total - sum(count))
  stopifnot(add >= 0L, add <= length(n))
  if (add > 0L) {
    i <- head(order(n - count, decreasing = TRUE), add)
    count[i] <- count[i] + 1
  }
  count
}

oe_population <- function(seed, grid_n, population_seed) {
  .observation_error_check_revision()
  stopifnot(grid_n %in% c(20L, 40L))
  grid <- SimSurvey::make_grid(res = rep(280 / grid_n, 2))
  history <- .observation_error_with_seed(population_seed, SimSurvey::sim_abundance(
    Z = function(ages, years) outer(0.3 + 0.06 * .oe_gillnet_q(ages), rep(1, length(years))),
    maturity = SimSurvey::sim_maturity(k = 1, x50 = 5), maturity_by = "age"))
  # The same marginal log-field variance and temporal persistence apply to all ages.
  fields <- list(coupled = SimSurvey::sim_ays_covar(group_ages = history$ages),
    default = SimSurvey::sim_ays_covar(),
    iid_age = SimSurvey::sim_ays_covar(sd = 2.8 / sqrt(1 - 0.5^2),
      phi_age = 0, group_ages = NULL))
  lapply(fields, function(field) {
    pop <- .observation_error_with_seed(seed,
      SimSurvey::sim_distribution(history, grid = grid, ays_covar = field))
    totals <- round(history$N)
    pop$sp_N[, N := .oe_balance(N, totals[.BY$age, .BY$year]), by = .(year, age)]
    pop <- SimSurvey::round_sim(pop)
    stopifnot(isTRUE(all.equal(pop$N, totals, tolerance = 0)))
    pop
  })
}

oe_overlap <- function(pop) {
  cells <- data.table::copy(pop$sp_N)[, .(cell, year, age, N)]
  cells[, share := if (sum(N) > 0) N / sum(N) else NA_real_, by = .(year, age)]
  summary <- cells[, {
    shares <- data.table::dcast(.SD, cell ~ age, value.var = "share")
    ages <- as.integer(names(shares)[-1L])
    weights <- as.matrix(shares[, -1L])
    i <- which(upper.tri(matrix(0, length(ages), length(ages))), arr.ind = TRUE)
    effective <- 1 / colSums(weights^2)
    list(age_a = ages[i[, 1L]], age_b = ages[i[, 2L]],
      overlap = crossprod(sqrt(weights))[i],
      effective_cells_a = effective[i[, 1L]], effective_cells_b = effective[i[, 2L]])
  }, by = year]
  map <- merge(cells[year %in% c(1, 10, 20) & age %in% c(2, 5, 12)],
    data.table::as.data.table(pop$grid_xy)[, .(cell, x, y)], by = "cell")
  list(summary = summary, map = map)
}

oe_snapshot <- function(pop, year) {
  i <- match(year, pop$years)
  stopifnot(length(i) == 1L, !is.na(i))
  keep <- pop$sp_N$year == year
  pop$sp_N <- pop$sp_N[keep]
  stopifnot(nrow(pop$sp_N) == length(pop$ages) * nrow(pop$grid_xy),
    all(pop$sp_N$year == year))
  pop$years <- year
  pop$N <- pop$N[, i, drop = FALSE]
  pop$Z <- pop$Z[, i, drop = FALSE]
  SimSurvey::round_sim(pop)
}

oe_gillnet_truth <- function(snapshot, seed) {
  fished <- .observation_error_with_seed(seed, SimSurvey::sim_fisheries(snapshot,
    Gillnet = SimSurvey::sim_fishery(target = 0.06, target_type = "F_full",
      selectivity = .oe_gillnet_q, catch_model = SimSurvey::sim_capture(q = 0.02,
        capacity_per_gear = 150), gear_amount = 60, duration = 24,
      sampling_coverage = 0, light = FALSE, maturity_sampling = "none")))
  f <- fished$fishery
  operation <- f$setdet[, .(operation, n)]
  caught <- f$joint_catch[, .(count = sum(n)), by = .(operation, age)]
  counts <- matrix(0, nrow(operation), length(snapshot$ages))
  if (nrow(caught)) {
    counts[cbind(match(caught$operation, operation$operation),
      match(caught$age, snapshot$ages))] <- caught$count
  }
  stopifnot(all(rowSums(counts) == operation$n))
  population <- SimSurvey::stock_biology(snapshot)$at_age[, .(year, age,
    weight = mean_weight_kg, maturity = prop_mature, N = as.numeric(snapshot$N) / 1000,
    Z = as.numeric(snapshot$Z), SSB = as.numeric(snapshot$S) / 1000)]
  F_realized <- f$reference$F_realized[match(population$age, f$reference$age)]
  stopifnot(all(is.finite(F_realized[population$N > 0])))
  population[, F := ifelse(N > 0, F_realized, 0)]
  # Integer rounding of the annual death budget can raise F above Z for very
  # rare ages. Retain realized F, and flag the corresponding signed remainder.
  population[, M := Z - F]
  annual <- population[, .(SSB = sum(SSB), recruitment = N[age == min(age)]), by = year]
  list(population = population, annual = annual, operation = operation, counts = counts,
    fishery_annual = f$annual, fishery_reference = f$reference)
}

.oe_random_ages <- function(counts, n) {
  left <- n
  available <- rowSums(counts)
  ages <- matrix(0, nrow(counts), ncol(counts))
  for (a in seq_len(ncol(counts) - 1L)) {
    i <- which(left > 0)
    ages[i, a] <- stats::rhyper(length(i), counts[i, a], available[i] - counts[i, a], left[i])
    left <- left - ages[, a]
    available <- available - counts[, a]
  }
  ages[, ncol(counts)] <- left
  ages
}

oe_gillnet_draw <- function(truth, seed, efforts = names(.oe_effort)) {
  stopifnot(length(efforts) > 0, !anyDuplicated(efforts), all(efforts %in% names(.oe_effort)))
  p <- truth$population
  n <- truth$operation$n
  cap <- pmin(n, 150)
  draw <- .observation_error_with_seed(seed, {
    u <- stats::runif(length(n))
    cv <- 0.05
    reports <- n * stats::rlnorm(length(n), -log1p(cv^2) / 2, sqrt(log1p(cv^2)))
    keep <- which(u < 0.3 * max(.oe_effort[efforts]) & n > 0)
    list(u = u, reports = reports, keep = keep,
      ages = .oe_random_ages(truth$counts[keep, , drop = FALSE], cap[keep]))
  })
  total <- sum(draw$reports)
  target <- colSums(truth$counts)
  data.table::rbindlist(lapply(efforts, function(effort) {
    probability <- 0.3 * .oe_effort[[effort]]
    kept <- draw$u[draw$keep] < probability
    i <- draw$keep[kept]
    weights <- draw$reports[i] / (probability * cap[i])
    expanded <- colSums(draw$ages[kept, , drop = FALSE] * weights)
    obs <- if (total == 0) rep(0, nrow(p)) else total * expanded / sum(expanded)
    status <- if (total == 0) "zero_reported_catch" else
      if (!length(i)) "no_measurement_support" else "estimated"
    rows <- data.table::data.table(year = p$year, age = p$age, component = "Gillnet",
      effort = effort, bin = "native", obs = ifelse(is.nan(obs), NA_real_, obs) / 1000,
      truth = target / 1000, se = NA_real_, n_primary = sum(draw$u < probability),
      n_primary_frame = length(n),
      n_caught = sum(n), n_measured = sum(cap[i]), n_aged = sum(cap[i]),
      n_caught_at_age = target, n_aged_at_age = colSums(draw$ages[kept, , drop = FALSE]),
      status = status, reported_total = total / 1000, N_frame = p$N,
      selectivity = .oe_gillnet_q(p$age),
      expected_encounters = probability * colSums(truth$counts * (cap / pmax(n, 1))))
    modeled <- data.table::copy(rows)
    modeled[, age := pmin(age, .oe_plus_age)]
    modeled <- modeled[, .(obs = sum(obs), truth = sum(truth),
      n_caught_at_age = sum(n_caught_at_age), n_aged_at_age = sum(n_aged_at_age),
      N_frame = sum(N_frame),
      selectivity = if (sum(N_frame) > 0) stats::weighted.mean(selectivity, N_frame) else NA_real_,
      expected_encounters = sum(expected_encounters)), by = age]
    info <- rows[age == 1L, setdiff(names(rows), names(modeled)), with = FALSE]
    rbind(rows, cbind(modeled, info)[, bin := "modeled"], use.names = TRUE)
  }), use.names = TRUE)
}

# Allocate exactly 80, 160 or 320 distinct cells, nested by stratum with at
# least two sets in each. Add the next set where the proportional deficit is greatest.
.oe_rv_allocation <- function(cells) {
  targets <- as.integer(160 * .oe_effort)
  stopifnot(all(cells >= 2L), 2L * length(cells) <= min(targets),
    sum(cells) >= max(targets))
  plan <- matrix(NA_integer_, length(cells), length(targets),
    dimnames = list(NULL, names(.oe_effort)))
  n <- rep(2L, length(cells))
  for (i in seq_along(targets)) {
    while (sum(n) < targets[i]) {
      deficit <- (sum(n) + 1L) * cells / sum(cells) - n
      deficit[n >= cells] <- -Inf
      h <- which.max(deficit)
      n[h] <- n[h] + 1L
    }
    plan[, i] <- n
  }
  plan
}

.oe_rv_design <- function(snapshot) {
  grid <- data.table::copy(data.table::as.data.table(snapshot$grid_xy))
  area <- abs(prod(stars::st_res(snapshot$grid)))
  grid[, rank := sample.int(.N), by = strat]
  strata <- grid[, .(n_cells = .N), by = strat][order(strat)]
  allocation <- .oe_rv_allocation(strata$n_cells)
  plan <- data.table::rbindlist(lapply(colnames(allocation), function(effort) {
    strata[, .(strat, n_cells, n_sets = allocation[, effort],
      n_equivalent = n_cells * area / .oe_tow_area, effort)]
  }))
  frame <- grid[, .(sim = 1L, year = snapshot$years, cell)]
  grid[strata, on = "strat", strat_cells := i.n_cells]
  grid[, c("sim", "year", "tow_area", "cell_area", "cell_sets", "strat_area") :=
    list(1L, snapshot$years, .oe_tow_area, area, 1L, strat_cells * area)]
  grid[, set := seq_len(.N)]
  list(grid = grid, plan = plan, frame = frame, area = area)
}

.oe_rv_truth <- function(snapshot, design, modeled) {
  grid <- design$grid[, .(cell, strat)]
  counts <- merge(snapshot$sp_N[, .(cell, age, N)], grid, by = "cell")
  counts[, p := .oe_rv_q(age) * .oe_tow_area / design$area]
  counts[, c("expected", "capture_var") := list(N * p, N * p * (1 - p))]
  if (modeled) {
    counts[, age := pmin(age, .oe_plus_age)]
    counts <- counts[, lapply(.SD, sum), by = .(strat, cell, age),
      .SDcols = c("N", "expected", "capture_var")]
  }
  counts[, .(N_frame = sum(N) / 1000, truth_tow = mean(expected),
    oracle_var_between = stats::var(expected), oracle_var_within = mean(capture_var)),
    by = .(strat, age)]
}

.oe_rv_bin <- function(survey, inputs, retained, plan, known, snapshot, effort, modeled) {
  strata <- c("sim", "year", "division", "strat", "strat_area", "tow_area")
  set_data <- survey$setdet[set %in% retained]
  expanded <- inputs$af[set %in% retained & !is.na(age), .(set, age, n)]
  captured <- survey$full_setdet[set %in% retained, .(set, age, n)]
  samples <- survey$samp[set %in% retained]
  if (modeled) {
    expanded <- expanded[, .(n = sum(n)), by = .(set, age = pmin(age, .oe_plus_age))]
    captured <- captured[, .(n = sum(n)), by = .(set, age = pmin(age, .oe_plus_age))]
    samples <- data.table::copy(samples)[, age := pmin(age, .oe_plus_age)]
  }
  expanded <- merge(set_data[, c("set", strata), with = FALSE], expanded, by = "set")
  reported <- SimSurvey::strat_means(expanded, "n", c(strata, "age"), c("sim", "year", "age"))
  counts <- samples[, .(n_aged_at_age = sum(aged)), by = age]
  sample_stats <- expanded[, .(mean_tow = mean(n), sd_tow = stats::sd(n)), by = .(strat, age)]
  oracle <- merge(captured, set_data[, .(set, strat)], by = "set")[
    , .(oracle_mean_tow = mean(n), n_caught_at_age = sum(n)), by = .(strat, age)]
  h <- merge(known, plan, by = "strat", allow.cartesian = TRUE)
  h <- merge(h, sample_stats, by = c("strat", "age"))
  h <- merge(h, oracle, by = c("strat", "age"))
  h[, c("area_weight", "fpc") := list(n_equivalent / sum(plan$n_equivalent), 1 - n_sets / n_cells)]
  totals <- h[, .(truth = sum(n_equivalent * truth_tow) / 1000, N_frame = sum(N_frame),
    n_caught_at_age = sum(n_caught_at_age), expected_encounters = sum(n_sets * truth_tow),
    oracle_obs = sum(n_equivalent * oracle_mean_tow) / 1000,
    oracle_se = sqrt(sum(n_equivalent^2 * (fpc * oracle_var_between + oracle_var_within) /
      n_sets)) / 1000, fpc = sum(area_weight * fpc)), by = age]
  totals[counts, on = "age", n_aged_at_age := i.n_aged_at_age]
  totals[is.na(n_aged_at_age), n_aged_at_age := 0L]
  status <- if (inputs$af[set %in% retained & is.na(age), any(n > 0, na.rm = TRUE)]) "partial_age" else "ok"
  rows <- merge(totals, reported[, .(year, age, obs = total / 1000,
    se = sd * sampling_units / 1000)], by = "age")
  if (status == "partial_age") rows[obs == 0, c("obs", "se") := list(NA_real_, NA_real_)]
  rows[, c("component", "effort", "bin", "n_primary", "n_caught", "n_measured", "n_aged",
      "n_equivalent", "n_primary_frame", "status", "selectivity") := list("RV", effort,
      if (modeled) "modeled" else "native", nrow(set_data), nrow(samples),
      sum(samples$measured), sum(samples$aged), sum(plan$n_equivalent), sum(plan$n_cells), status,
      ifelse(N_frame > 0, truth / N_frame, .oe_rv_q(age)))]
  h[, c("year", "bin") := list(snapshot$years, if (modeled) "modeled" else "native")]
  list(rows = rows, strata = h)
}

oe_rv_draw <- function(snapshot, seed, efforts = names(.oe_effort)) {
  stopifnot(length(efforts) > 0, !anyDuplicated(efforts), all(efforts %in% names(.oe_effort)))
  .observation_error_with_seed(seed, {
    design <- .oe_rv_design(snapshot)
    maximum <- efforts[which.max(.oe_effort[efforts])]
    plan <- design$plan[effort == maximum, .(strat, n_sets)]
    sets <- merge(design$grid, plan, by = "strat")[rank <= n_sets]
    sets[, strat_sets := n_sets]
    survey <- SimSurvey::sim_survey(snapshot, q = .oe_rv_q, custom_sets = sets,
      survey_frame = design$frame, lengths_cap = 100L, ages_cap = 3L,
      age_length_group = 3L, age_space_group = "set", maturity_sampling = "none", light = FALSE)
    length_group <- get("length_group", envir = environment(survey$sim_length), inherits = FALSE)
    inputs <- SimSurvey::strat_data(survey, length_group = length_group, alk_scale = "set")
    by_bin <- lapply(c(FALSE, TRUE), function(modeled) {
      known <- .oe_rv_truth(snapshot, design, modeled)
      lapply(efforts, function(level) {
        plan <- design$plan[effort == level]
        retained <- merge(design$grid[, .(set, strat, rank)], plan, by = "strat")[rank <= n_sets, set]
        .oe_rv_bin(survey, inputs, retained, plan, known, snapshot, level, modeled)
      })
    })
    pieces <- unlist(by_bin, recursive = FALSE)
    rows <- data.table::rbindlist(lapply(pieces, "[[", "rows"))
    native <- as.numeric(snapshot$N) * .oe_rv_q(snapshot$ages) / 1000
    expected <- native[match(rows$age, snapshot$ages)]
    expected[rows$bin == "modeled" & rows$age == .oe_plus_age] <-
      sum(native[snapshot$ages >= .oe_plus_age])
    stopifnot(isTRUE(all.equal(rows$truth, expected, tolerance = 1e-12)))
    list(rows = rows,
      strata = data.table::rbindlist(lapply(pieces, "[[", "strata")))
  })
}
