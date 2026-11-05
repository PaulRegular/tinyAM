# Domain calculations used by data-raw/sim_cod_obs.R.
# These helpers neither load SimSurvey nor start a simulation when sourced.

# Conditional expansion of measured/aged fish, followed by operation weighting
# and a ratio calibration to the reported annual fleet total. This deliberately
# accepts no true catch-at-age or unobserved true-age column.
raise_sim_cod_catch <- function(operations, samples, ages, years, fleets,
                                age_length_group = 3) {
  operations <- as.data.frame(operations)
  samples <- as.data.frame(samples)
  required_ops <- c("fishery", "year", "operation", "reported_n", "sampling_probability", "sampled")
  required_samples <- c("operation", "length", "age", "aged")
  stopifnot(all(required_ops %in% names(operations)), all(required_samples %in% names(samples)))
  if (anyDuplicated(operations$operation)) stop("Operation identifiers must be unique.")
  sample_rows <- split(seq_len(nrow(samples)), samples$operation)
  out <- monitors <- list()
  for (fleet in fleets) for (year in years) {
    ops <- operations[operations$fishery == fleet & operations$year == year, , drop = FALSE]
    total <- sum(ops$reported_n)
    if (!is.finite(total) || total < 0) stop("Reported annual catches must be finite and nonnegative.")
    h <- numeric(length(ages))
    status <- if (total == 0) "zero_total" else "ok"
    sampled <- ops[ops$sampled & ops$reported_n > 0, , drop = FALSE]
    if (total > 0 && !nrow(sampled)) status <- "no_samples"
    for (j in seq_len(nrow(sampled))) {
      op <- sampled[j, ]
      s <- samples[sample_rows[[as.character(op$operation)]], , drop = FALSE]
      if ("measured" %in% names(s)) s <- s[s$measured, , drop = FALSE]
      if (!nrow(s) || !is.finite(op$sampling_probability) || op$sampling_probability <= 0) {
        status <- "unsupported_operation"
        next
      }
      if (any(!is.finite(s$length))) stop("Measured lengths must be finite.")
      bin <- floor(s$length / age_length_group)
      age_counts <- numeric(length(ages))
      for (group in unique(bin)) {
        measured <- which(bin == group)
        aged <- measured[s$aged[measured] & !is.na(s$age[measured])]
        if (!length(aged)) {
          status <- "unaged_length_bin"
          next
        }
        if (any(!s$age[aged] %in% ages)) stop("Sample ages lie outside the declared age grid.")
        age_counts <- age_counts + length(measured) *
          tabulate(match(s$age[aged], ages), nbins = length(ages)) / length(aged)
      }
      h <- h + op$reported_n * age_counts / nrow(s) / op$sampling_probability
    }
    estimate <- if (total == 0) h else if (status == "ok" && sum(h) > 0) total * h / sum(h) else rep(NA_real_, length(ages))
    out[[length(out) + 1L]] <- data.frame(year = year, age = ages, fleet = fleet, obs = estimate)
    ss <- samples[samples$operation %in% ops$operation, , drop = FALSE]
    monitors[[length(monitors) + 1L]] <- data.frame(fleet = fleet, year = year,
      reported_total = total, estimated_total = if (all(is.na(estimate))) NA_real_ else sum(estimate),
      n_operations = nrow(ops), n_sampled = sum(ops$sampled), n_measured = nrow(ss),
      n_aged = sum(ss$aged & !is.na(ss$age)), status = status)
  }
  list(catch = do.call(rbind, out), monitoring = do.call(rbind, monitors))
}

.sim_cod_hash <- function(x) {
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(x, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

.sim_cod_git <- function(path) {
  result <- system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout = TRUE)
  if (length(result) != 1L || !grepl("^[[:xdigit:]]{40}$", result)) stop("Cannot identify the SimSurvey source revision.")
  result
}

.sim_cod_grid <- function(settings) {
  g <- settings$grid
  original <- SimSurvey::make_grid(x_range = g$x_range, y_range = g$y_range, res = g$res,
    shelf_depth = g$shelf_depth, shelf_width = g$shelf_width, depth_range = g$depth_range,
    n_div = g$n_div, strat_breaks = g$strat_breaks, strat_splits = g$strat_splits)
  grid <- as.data.frame(original)
  # Rigid rotation: old shallow western edge becomes the northern coast.
  old_x <- grid$x
  grid$x <- grid$y
  grid$y <- -old_x
  grid$Gillnet <- grid$y >= 40 & grid$depth >= 25 & grid$depth <= 250
  grid$Longline <- grid$y >= 0 & grid$depth >= 25 & grid$depth <= 350
  grid$Trawl <- grid$y < 70 & grid$depth >= 50 & grid$depth <= 550
  grid$RV <- grid$y < 70 & grid$depth >= 50 & grid$depth <= 550
  grid$Sentinel <- grid$y >= 40 & grid$depth >= 25 & grid$depth <= 250
  result <- stars::st_as_stars(grid, crs = sf::st_crs(original))
  stopifnot(identical(sort(as.numeric(result$depth)), sort(as.numeric(original$depth))))
  result
}

.sim_cod_selectivity <- function(settings) {
  p <- settings$selectivity
  list(
    Gillnet = function(x) stats::plogis((x - p$Gillnet[["up"]]) / p$Gillnet[["width"]]) *
      stats::plogis((p$Gillnet[["down"]] - x) / p$Gillnet[["down_width"]]),
    Longline = function(x) stats::plogis((x - p$Longline[["up"]]) / p$Longline[["width"]]),
    Trawl = function(x) stats::plogis((x - p$Trawl[["up"]]) / p$Trawl[["width"]]),
    RV = function(x) stats::plogis((x - p$RV[["up"]]) / p$RV[["width"]]),
    Sentinel = function(x) stats::plogis((x - p$Sentinel[["up"]]) / p$Sentinel[["width"]]) *
      stats::plogis((p$Sentinel[["down"]] - x) / p$Sentinel[["down_width"]]))
}

# Match sim_sets(min_sets = 2, resample_cells = FALSE), without drawing sets.
.sim_cod_check_surveys <- function(grid, settings, cell_area) {
  for (survey in names(settings$survey_density)) {
    cells <- table(as.character(grid$strat[grid[[survey]]]))
    if (!length(cells)) stop(survey, " has no eligible survey cells.")
    requested <- pmax(2, round(as.numeric(cells) * cell_area * settings$survey_density[[survey]]))
    bad <- which(requested > cells)
    if (length(bad)) {
      i <- bad[[1L]]
      stop(survey, " requests ", requested[[i]], " sets from ", cells[[i]],
        " cells in stratum ", names(cells)[[i]],
        ". Reduce survey_density or use smaller grid cells; cells are sampled without replacement.")
    }
  }
  invisible(TRUE)
}

.sim_cod_available <- function(population, settings, selectivity) {
  grid <- as.data.frame(population$grid)
  stock <- as.data.frame(population$sp_N)
  key <- population$sim_length(age = settings$ages, length_age_key = TRUE)
  lengths <- as.numeric(rownames(key))
  answer <- list()
  for (label in names(selectivity)) {
    sa <- colSums(key * selectivity[[label]](lengths))
    selected <- stock$cell %in% grid$cell[grid[[label]]]
    regional <- stats::aggregate(stock$N[selected], stock[selected, c("year", "age")], sum)
    names(regional)[3] <- "N_region"
    regional$program <- label
    regional$N_total <- population$N[cbind(match(regional$age, settings$ages), match(regional$year, settings$years))]
    regional$selection <- sa[match(regional$age, settings$ages)]
    regional$index_expected <- regional$N_region * regional$selection
    regional$availability <- regional$N_region / regional$N_total
    answer[[label]] <- regional
  }
  do.call(rbind, answer)
}

.sim_cod_allocations <- function(population, settings, available) {
  result <- list()
  for (fleet in names(settings$shares)) for (year in settings$years) {
    rows <- available[available$program == fleet & available$year == year, ]
    z <- population$Z[cbind(match(rows$age, settings$ages), match(rows$year, settings$years))]
    global <- sum(rows$N_total * rows$selection * -expm1(-z))
    accessible <- sum(rows$N_region * rows$selection * -expm1(-z))
    if (!is.finite(global) || !is.finite(accessible) || global <= 0 || accessible <= 0 || accessible > global) {
      stop("Every active fleet/year requires positive accessible and global vulnerable deaths: ", fleet, " / ", year)
    }
    fraction <- accessible / global
    base <- unname(settings$shares[[fleet]])
    result[[length(result) + 1L]] <- data.frame(fleet = fleet, year = year,
      base_share = base, fishery_share = base * fraction, access_fraction = fraction,
      global_vulnerable_deaths = global, accessible_vulnerable_deaths = accessible,
      target = base * accessible)
  }
  do.call(rbind, result)
}

# Optional scenario design: one common fishing-allocation multiplier sets the
# geometric centre of the reference M surface. Use the baseline Z, not realized
# Z states, observations, fitted SSB, or a search over seeds. Realized M remains
# an outcome of fixed-Z fishing and need not equal this reference surface.
.sim_cod_reference_M <- function(settings, available, allocations) {
  target <- settings$reference_M
  baseline_Z <- exp(settings$mortality$log_mean)
  if (length(target) != 1L || !is.finite(target) || target <= 0 ||
      length(baseline_Z) != 1L || !is.finite(baseline_Z) || target >= baseline_Z) {
    stop("reference_M must be positive and below the scalar baseline Z; positive fleet allocations cannot attain the baseline itself.")
  }
  base <- settings$shares
  if (any(!is.finite(base)) || any(base <= 0) || any(base > 1)) {
    stop("Base fleet shares must be positive fractions no greater than one.")
  }
  grid <- expand.grid(age = settings$ages, year = settings$years)
  keys <- paste(grid$year, grid$age)
  fraction <- numeric(nrow(grid))
  for (fleet in names(base)) {
    rows <- available[available$program == fleet, ]
    age_selection <- rows$selection[match(keys, paste(rows$year, rows$age))]
    allocation <- allocations[allocations$fleet == fleet, ]
    p <- allocation$fishery_share[match(grid$year, allocation$year)]
    fraction <- fraction + p * age_selection
  }
  if (any(!is.finite(fraction)) || any(fraction < 0) || !any(fraction > 0)) {
    stop("Reference-M calibration requires finite positive fishing support.")
  }
  # Native shares and local base shares remain <=1. Stay inside the strictly
  # positive joint remaining-mortality boundary, at representable precision.
  upper <- min(1 / max(base), 1 / max(allocations$fishery_share),
    (1 - 16 * .Machine$double.eps) / max(fraction))
  mean_log_M <- function(scale) mean(log(baseline_Z * (1 - scale * fraction)))
  minimum_M <- exp(mean_log_M(upper))
  if (!is.finite(minimum_M) || target < minimum_M) {
    stop("reference_M = ", format(target), " is not attainable with the current relative fleet shares and selectivity. ",
      "The lowest feasible geometric reference M is ", format(minimum_M),
      "; change the scenario explicitly rather than relaxing allocation bounds.")
  }
  multiplier <- stats::uniroot(function(scale) mean_log_M(scale) - log(target),
    interval = c(0, upper), tol = 1e-10)$root
  scaled <- base * multiplier
  native <- allocations$fishery_share * multiplier
  remaining <- 1 - multiplier * fraction
  if (multiplier <= 0 || any(scaled <= 0 | scaled > 1) ||
      any(native <= 0 | native > 1) || any(remaining <= 0)) {
    stop("Calibrated allocations violate positive remaining mortality or a fleet-share bound.")
  }
  calibrated <- allocations
  calibrated$base_share <- scaled[match(calibrated$fleet, names(base))]
  calibrated$fishery_share <- native
  calibrated$target <- calibrated$base_share * calibrated$accessible_vulnerable_deaths
  grid$M0 <- baseline_Z * remaining
  list(allocations = calibrated, calibration = list(reference_M = target,
    basis = "Z process median/log-scale centre; target is the geometric centre of the M baseline, not arithmetic realized M",
    baseline_Z = baseline_Z, initial_shares = base, multiplier = multiplier,
    calibrated_shares = scaled, M0 = grid,
    feasibility = list(max_multiplier = upper, minimum_reference_M = minimum_M,
      achieved_reference_M = exp(mean(log(grid$M0))), minimum_remaining_fraction = min(remaining),
      maximum_base_share = max(scaled), maximum_native_share = max(native),
      maximum_joint_fishing_fraction = max(1 - remaining))))
}

# Only the fishing copy omits exact-zero spatial rows. One zero anchor per
# year-age keeps even entirely empty groups explicit for round_sim(). The full
# grid remains the location frame; surveys independently use the dense original.
.sim_cod_fishing_population <- function(population) {
  spatial <- population$sp_N
  keep <- spatial$N > 0 | !duplicated(spatial[, c("year", "age")])
  result <- population
  result$sp_N <- data.table::copy(spatial[keep, ])
  totals <- tapply(result$sp_N$N,
    list(age = result$sp_N$age, year = result$sp_N$year), sum)
  if (!isTRUE(all.equal(totals, population$N, tolerance = 0, check.attributes = FALSE)) ||
      !identical(result$grid, population$grid)) {
    stop("Sparse fishing representation must preserve every population total and the complete grid.")
  }
  result
}

.sim_cod_tables <- function(stock, fishery, surveys, settings, selectivity, available, provenance, allocations,
                            calibration = NULL) {
  ages <- settings$ages
  years <- settings$years
  fleets <- names(settings$shares)
  coordinates <- expand.grid(age = ages, year = years)
  raised <- raise_sim_cod_catch(fishery$setdet, fishery$samp,
    ages, years, fleets, settings$monitoring$age_length_group)
  catch <- raised$catch
  catch$obs <- catch$obs / 1000
  catch$unit <- "thousand fish"
  index <- survey_precision <- survey_monitoring <- list()
  for (survey in names(settings$survey_density)) {
    selected <- SimSurvey::select_survey(surveys, survey)
    stratified <- SimSurvey::run_strat(selected, alk_scale = "set")
    table <- as.data.frame(stratified$age_strat)
    # Keep each survey's own footprint area; no expansion to the whole grid.
    table <- merge(coordinates, table, by = c("year", "age"), all.x = TRUE, sort = FALSE)
    index[[survey]] <- data.frame(year = table$year, age = table$age, obs = table$total / 1000,
      survey = survey, samp_time = 0, unit = "thousand fish")
    survey_precision[[survey]] <- data.frame(year = table$year, age = table$age,
      survey = survey, sets = table$sets, total_se = table$sd * table$sampling_units / 1000)
    survey_monitoring[[survey]] <- transform(as.data.frame(selected$samp_totals), survey = survey)
  }
  index <- do.call(rbind, index)
  key <- stock$sim_length(age = ages, length_age_key = TRUE)
  lengths <- as.numeric(rownames(key))
  weights <- stock$sim_weight(lengths)
  mean_weight <- as.numeric(crossprod(weights, key))
  mean_length <- as.numeric(crossprod(lengths, key))
  maturity_at_age <- stock$sim_maturity(ages)
  stock_weight <- data.frame(year = coordinates$year, age = coordinates$age,
    obs = mean_weight[match(coordinates$age, ages)], type = "stock", fleet = NA_character_, unit = "kg")
  catches <- as.data.frame(fishery$catch_at_age)
  catch_weight <- list()
  fleet_truth <- as.data.frame(fishery$reference)
  names(fleet_truth)[names(fleet_truth) == "fishery"] <- "fleet"
  for (fleet in fleets) {
    rows <- catches[catches$fishery == fleet, , drop = FALSE]
    summed <- stats::aggregate(rows[c("n", "kg")], rows[c("year", "age")], sum)
    tab <- merge(coordinates, summed, by = c("year", "age"), all.x = TRUE, sort = FALSE)
    tab$n[is.na(tab$n)] <- 0
    tab$kg[is.na(tab$kg)] <- 0
    selected_key <- key * selectivity[[fleet]](lengths)
    expected_weight <- as.numeric(crossprod(weights, selected_key)) / colSums(selected_key)
    tab$weight <- ifelse(tab$n > 0, tab$kg / tab$n, expected_weight[match(tab$age, ages)])
    catch_weight[[fleet]] <- data.frame(year = tab$year, age = tab$age, obs = tab$weight,
      type = "catch", fleet = fleet, unit = "kg")
    position <- match(paste(fleet_truth$year, fleet_truth$age), paste(tab$year, tab$age))
    use <- fleet_truth$fleet == fleet
    fleet_truth$weight[use] <- tab$weight[position[use]]
    fleet_truth$yield[use] <- tab$kg[position[use]] / 1000
  }
  maturity <- data.frame(year = coordinates$year, age = coordinates$age,
    obs = maturity_at_age[match(coordinates$age, ages)], spawn_time = 0)
  inputs <- list(catch = catch, index = index,
    weight = rbind(stock_weight, do.call(rbind, catch_weight)), maturity = maturity)
  if (!is.null(calibration)) {
    # A known growth-derived covariate, not a hidden abundance or mortality state.
    inputs$catch$mean_length <- mean_length[match(inputs$catch$age, ages)]
    inputs$index$mean_length <- mean_length[match(inputs$index$age, ages)]
    inputs$pop_covariates <- calibration$M0[c("year", "age", "M0")]
  }
  dat <- tinyAM::as_tam_data(inputs)
  tinyAM::check_tam_data(dat, quiet = TRUE)
  total_f <- stats::aggregate(fleet_truth$F_realized, fleet_truth[c("year", "age")], sum)
  names(total_f)[3] <- "F"
  population <- coordinates
  population$N <- as.numeric(stock$N) / 1000
  population$Z <- as.numeric(stock$Z)
  population$weight <- mean_weight[match(population$age, ages)]
  population$maturity <- maturity_at_age[match(population$age, ages)]
  population$B <- population$N * population$weight
  population$SSB <- population$B * population$maturity
  population$F <- total_f$F[match(paste(population$year, population$age), paste(total_f$year, total_f$age))]
  population$M <- population$Z - population$F
  if (!is.null(calibration)) {
    realized <- population$M[is.finite(population$M)]
    provenance$reference_M$realized_M_summary <- list(
      arithmetic_mean = if (length(realized)) mean(realized) else NA_real_,
      geometric_mean = if (length(realized) && all(realized > 0)) exp(mean(log(realized))) else NA_real_,
      minimum = if (length(realized)) min(realized) else NA_real_,
      maximum = if (length(realized)) max(realized) else NA_real_,
      nonpositive_cells = sum(realized <= 0), missing_cells = sum(!is.finite(population$M)))
  }
  annual <- stats::aggregate(population[c("N", "B", "SSB")], population["year"], sum)
  annual$recruitment <- population$N[match(paste(annual$year, min(ages)), paste(population$year, population$age))]
  realized_yield <- stats::aggregate(fleet_truth$yield, fleet_truth["year"], sum)
  annual$yield <- realized_yield$x[match(annual$year, realized_yield$year)]
  available$N_region <- available$N_region / 1000
  available$N_total <- available$N_total / 1000
  available$index_expected <- available$index_expected / 1000
  regional <- available[available$program %in% names(settings$survey_density), ]
  names(regional)[names(regional) == "program"] <- "survey"
  footprint <- as.data.frame(stock$grid)
  footprint$area_km2 <- abs(prod(stars::st_res(stock$grid)))
  curves <- do.call(rbind, lapply(names(selectivity), function(label) {
    data.frame(component = label, length = seq(0, 140), selectivity = selectivity[[label]](seq(0, 140)))
  }))
  for (column in intersect(c("N", "catch_reference", "catch_realized", "total_deaths", "death_budget"), names(fleet_truth))) {
    fleet_truth[[column]] <- fleet_truth[[column]] / 1000
  }
  loss <- population$N[population$age == max(ages)] * exp(-population$Z[population$age == max(ages)])
  last_weight <- population$weight[population$age == max(ages)]
  last_maturity <- population$maturity[population$age == max(ages)]
  checks <- list(terminal_survivors = data.frame(year = years, N = loss,
      proportion = loss / annual$N, B = loss * last_weight,
      SSB = loss * last_weight * last_maturity,
      B_proportion = loss * last_weight / annual$B,
      SSB_proportion = loss * last_weight * last_maturity / annual$SSB),
    negative_remaining_M = population[is.finite(population$M) & population$M < 0, c("year", "age", "M")],
    shortfalls = as.data.frame(fishery$annual)[c("fishery", "year", "target", "n", "shortfall", "status")],
    unsupported_catch = raised$monitoring[raised$monitoring$status != "ok" & raised$monitoring$status != "zero_total", ])
  accessibility <- lapply(fleets, function(fleet) {
    regional <- available[available$program == fleet, ]
    z <- population$Z[match(paste(regional$year, regional$age), paste(population$year, population$age))]
    allocation <- allocations[allocations$fleet == fleet, ]
    # Separate-fleet necessary bound, not a guarantee under competition. Use
    # ceil(expected deaths) as an upper bound on balanced integer age budgets.
    regional$target <- regional$N_total * allocation$fishery_share[match(regional$year, allocation$year)] *
      regional$selection * -expm1(-z)
    regional$accessible_expected_deaths <- regional$N_region * regional$selection * -expm1(-z)
    regional$accessible_bound <- pmin(regional$N_region,
      ceiling(regional$N_total * 1000 * -expm1(-z)) / 1000)
    table <- stats::aggregate(regional[c("target", "accessible_bound", "accessible_expected_deaths")], regional["year"], sum)
    table$fleet <- fleet
    table$target_exceeds_bound <- table$target > table$accessible_bound
    table$local_reference_fraction <- table$target / table$accessible_expected_deaths
    table
  })
  checks$accessibility <- do.call(rbind, accessibility)
  truth <- list(population = population, annual = annual, fleet = fleet_truth,
    surveys = regional, footprint = footprint, selectivity = curves,
    fishery_annual = as.data.frame(fishery$annual),
    monitoring = raised$monitoring, survey_precision = do.call(rbind, survey_precision),
    survey_monitoring = do.call(rbind, survey_monitoring),
    allocations = allocations,
    checks = checks, provenance = provenance)
  list(sim_cod_obs = dat, sim_cod_truth = truth)
}
