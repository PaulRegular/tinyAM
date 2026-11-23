# Default SimSurvey mortality histories split between fishing and other deaths.
delta_matched_sets <- c(low = 80L, standard = 160L, high = 320L)
delta_matched_coverage <- c(low = .25, standard = .50, high = 1)
delta_matched_shares <- c(high_f = .90, low_f = .10)
delta_matched_flat <- function(age) 1 + 0 * age

delta_matched_survey_sets <- function(snapshot, sizes) {
  cells <- data.table::as.data.table(as.data.frame(snapshot$grid))
  area <- prod(stars::st_res(snapshot$grid))
  groups <- split(seq_len(nrow(cells)), cells$strat, drop = TRUE)
  available <- lengths(groups)
  if (any(available < 2L) || min(sizes) < 2L * length(groups) ||
      max(sizes) > nrow(cells) || is.unsorted(sizes)) {
    stop("The requested survey design cannot sample every stratum without replacement.")
  }

  # Add sets in proportion to stratum size. Each effort uses a prefix of the
  # same random permutation within a stratum.
  allocated <- rep.int(2L, length(groups))
  shuffled <- lapply(groups, function(rows) rows[sample.int(length(rows))])
  selections <- lapply(sizes, function(size) {
    for (set in seq_len(size - sum(allocated))) {
      priority <- ifelse(allocated < available, available / (allocated + 1L), -Inf)
      stratum <- which.max(priority)
      allocated[stratum] <<- allocated[stratum] + 1L
    }
    unlist(Map(function(rows, n) rows[seq_len(n)], shuffled, allocated), use.names = FALSE)
  })
  selected <- selections[[length(selections)]]
  sets <- cells[selected, ]
  sets$sim <- 1L
  sets$year <- snapshot$years
  sets$strat_cells <- unname(available[as.character(sets$strat)])
  sets$tow_area <- 0.0018
  sets$cell_area <- area
  sets$strat_area <- sets$strat_cells * area
  sets[, strat_sets := .N, by = "strat"]
  sets$cell_sets <- 1L
  sets$set <- seq_len(nrow(sets))
  frame <- data.table::data.table(sim = 1L, year = snapshot$years, cell = cells$cell)
  data.table::setattr(sets, "survey_frame", frame)
  list(sets = sets, select = lapply(selections, match, table = selected))
}

delta_matched_monitor_survey <- function(snapshot, sizes) {
  design <- delta_matched_survey_sets(snapshot, sizes)
  sampled <- SimSurvey::sim_survey(snapshot, custom_sets = design$sets,
    q = delta_matched_flat, lengths_cap = 100, ages_cap = 3, age_length_group = 3,
    age_space_group = "set", maturity_sampling = "none", light = FALSE)
  lapply(design$select, function(selected) {
    survey <- sampled
    survey$setdet <- sampled$setdet[sampled$setdet$set %in% selected, ]
    survey$setdet[, strat_sets := .N, by = "strat"]
    survey$samp <- sampled$samp[sampled$samp$set %in% selected, ]
    survey$full_setdet <- sampled$full_setdet[sampled$full_setdet$set %in% selected, ]
    survey$samp_totals <- survey$setdet[, list(n_sets = .N, n_caught = sum(n),
      kg_caught = sum(kg), n_measured = sum(n_measured), n_aged = sum(n_aged)),
      by = c("sim", "year")]
    index <- as.data.frame(SimSurvey::run_strat(survey, alk_scale = "set")$age_strat)
    index <- index[!is.na(index$age), , drop = FALSE]
    unresolved_zero <- index$estimate_status %in% "partial_age" &
      !is.na(index$total) & index$total == 0
    index$total[unresolved_zero] <- NA_real_
    stopifnot(nrow(survey$setdet) == length(selected),
      identical(as.integer(index$age), as.integer(snapshot$ages)))
    list(data = data.frame(index[c("year", "age")], survey = "RV", samp_time = 0,
      obs = index$total / 1000, unit = "thousand fish"),
      sampling = data.frame(component = "index", year = snapshot$years,
        n_aged = sum(survey$setdet$n_aged), n_measured = sum(survey$setdet$n_measured),
        n_sets = nrow(survey$setdet), n_operations = NA_integer_,
        n_sampled_operations = NA_integer_))
  })
}

delta_matched_monitor_catch <- function(raw, spec, coverage) {
  spec$sampling_coverage <- 1
  sampled <- delta_monitor_fishery(raw, spec)
  operations <- sampled$setdet$operation
  order <- sample.int(length(operations))
  lapply(coverage, function(fraction) {
    fishery <- sampled
    if (fraction < 1 && length(operations)) {
      # Fixed-size simple random sampling of operations. Reporting multipliers
      # and biological samples remain shared across nested selections.
      count <- max(1L, as.integer(round(fraction * length(operations))))
      selected <- operations[order[seq_len(count)]]
      fishery$setdet <- data.table::copy(sampled$setdet)
      fishery$setdet$sampled <- operations %in% selected
      fishery$setdet$sampling_probability <- rep(count / length(operations), length(operations))
      fishery$samp <- sampled$samp[sampled$samp$operation %in% selected, ]
      raw$fishery <- fishery
      fishery <- SimSurvey::run_catch(raw)$fishery
    }
    catch <- as.data.frame(fishery$catch_estimates$catch_at_age)
    diagnostics <- fishery$catch_diagnostics
    stopifnot(identical(as.integer(catch$age), as.integer(raw$ages)))
    list(data = data.frame(catch[c("year", "age")], obs = catch$n_hat / 1000,
      unit = "thousand fish"),
      sampling = data.frame(component = "catch", year = raw$years,
        n_aged = diagnostics$n_aged, n_measured = diagnostics$n_measured, n_sets = NA_integer_,
        n_operations = diagnostics$n_operations, n_sampled_operations = diagnostics$n_sampled_operations))
  })
}

delta_matched_spatial <- function(id, years = delta_years, grid_n = delta_grid) {
  seeds <- c(population = 104000L, spatial = 204000L, fishery = 304000L,
    catch_train = 404000L, catch_holdout = 504000L,
    index_train = 604000L, index_holdout = 704000L) + 100L * id
  history <- withr::with_seed(seeds[["population"]], SimSurvey::sim_abundance(years = years))
  terminal <- nrow(history$N)
  for (j in seq_along(years)[-1L]) {
    history$N[terminal, j] <- history$N[terminal, j] +
      history$N[terminal, j - 1L] * exp(-history$Z[terminal, j - 1L])
  }
  expected <- round(history$N)
  grid <- SimSurvey::make_grid(res = rep(280 / grid_n, 2L))
  pop <- withr::with_seed(seeds[["spatial"]], SimSurvey::sim_distribution(history, grid = grid))
  pop$sp_N[, N := .oe_balance(N, expected[match(.BY$age, pop$ages), match(.BY$year, years)]),
    by = .(year, age)]
  pop <- SimSurvey::round_sim(pop)
  stopifnot(isTRUE(all.equal(pop$N, expected, tolerance = 0)))
  biology <- as.data.frame(SimSurvey::stock_biology(pop)$at_age)
  population <- data.frame(biology[c("year", "age")], N = as.numeric(pop$N) / 1000,
    Z = as.numeric(pop$Z), ssb = as.numeric(pop$S) / 1000)
  population$index_expected <- population$N
  weight <- data.frame(biology[c("year", "age")], obs = biology$mean_weight_kg, unit = "kg")
  maturity <- data.frame(biology[c("year", "age")],
    obs = biology$mean_mature_weight_kg / biology$mean_weight_kg, spawn_time = 0)
  simulator <- delta_fishery_simulator(delta_fast_fishery_engine)
  spec <- SimSurvey::sim_fishery(target = delta_matched_shares[[1L]], target_type = "share",
    selectivity = delta_matched_flat,
    catch_model = SimSurvey::sim_capture(q = .02, capacity_per_gear = 150),
    gear_amount = 1200L, duration = 24, gear_unit = "net", time_unit = "hour",
    sampling_coverage = 1, lengths_cap = 150, ages_cap = 150, age_sampling = "random",
    reporting_cv = .05, maturity_sampling = "none", light = FALSE)
  efforts <- names(delta_matched_sets)
  replicates <- c("train", "holdout")
  gather <- function(x) as.data.frame(data.table::rbindlist(x))
  batches <- lapply(seq_along(years), function(j) {
    year <- years[j]
    snapshot <- oe_snapshot(pop, year)
    survey <- lapply(replicates, function(rep) withr::with_seed(seeds[[paste0("index_", rep)]] + year,
      delta_matched_monitor_survey(snapshot, delta_matched_sets)))
    names(survey) <- replicates
    run_fishery <- function(target) {
      spec$target <- target
      raw <- withr::with_seed(seeds[["fishery"]] + year, simulator(snapshot, Gillnet = spec))
      ledger <- as.data.frame(raw$fishery$reference)
      ledger <- ledger[match(pop$ages, ledger$age), , drop = FALSE]
      stopifnot(all(ledger$catch_realized <= ledger$death_budget),
        abs(sum(ledger$catch_realized) - raw$fishery$annual$n) < 1e-8)
      truth <- data.frame(ledger[c("year", "age", "F_reference")],
        F_realized = ifelse(ledger$N > 0, ledger$F_realized, 0),
        catch_expected = ledger$catch_realized / 1000, death_budget = ledger$death_budget / 1000)
      catcher <- lapply(replicates, function(rep) withr::with_seed(seeds[[paste0("catch_", rep)]] + year,
        delta_matched_monitor_catch(raw, spec, delta_matched_coverage)))
      names(catcher) <- replicates
      annual <- as.data.frame(raw$fishery$annual[, c("year", "target", "achieved", "shortfall")])
      list(truth = truth, annual = annual, monitoring = lapply(replicates, function(rep) {
        stats::setNames(lapply(efforts, function(effort) {
          s <- survey[[rep]][[effort]]
          c <- catcher[[rep]][[effort]]
          list(catch = c$data, index = s$data,
            sampling = transform(rbind(c$sampling, s$sampling), replicate = rep))
        }), efforts)
      }))
    }
    lapply(delta_matched_shares, run_fishery)
  })

  bundles <- list()
  for (driver in names(delta_matched_shares)) {
    draws <- lapply(batches, `[[`, driver)
    reference <- merge(population, gather(lapply(draws, `[[`, "truth")), by = c("year", "age"))
    reference$M_reference <- reference$Z - reference$F_reference
    reference$M_realized <- reference$Z - reference$F_realized
    reference$F_share <- with(reference, ifelse(Z > 0, F_realized / Z, NA_real_))
    reference$deaths_expected <- with(reference, N * -expm1(-Z))
    # The rounded age budgets can exceed continuous expected deaths by less
    # than one fish. Preserve the signed M estimates for the audit.
    excess <- reference$catch_expected - reference$deaths_expected
    severe <- any(excess > .001 + 1e-10)
    target_share <- delta_matched_shares[[driver]]
    stopifnot(isTRUE(all.equal(reference$F_reference, target_share * reference$Z,
      tolerance = 1e-12)), !severe, is.finite(mean(reference$M_realized)),
      mean(reference$M_realized) > 0)
    reconstructed <- with(reference, N * F_realized / Z * -expm1(-Z))
    annual <- gather(lapply(draws, `[[`, "annual"))
    annual$target_share <- target_share
    annual$realized_share <- with(reference, unname(tapply(catch_expected, year, sum) /
      tapply(deaths_expected, year, sum)))[match(annual$year, sort(unique(reference$year)))]
    annual$F_mean <- with(reference, unname(tapply(F_realized, year, mean)))[match(annual$year, sort(unique(reference$year)))]
    annual$M_mean <- with(reference, unname(tapply(M_realized, year, mean)))[match(annual$year, sort(unique(reference$year)))]
    annual$Z_mean <- with(reference, unname(tapply(Z, year, mean)))[match(annual$year, sort(unique(reference$year)))]
    audit <- data.frame(id = id, driver = driver, valid =
      isTRUE(all.equal(reconstructed, reference$catch_expected, tolerance = 1e-12)),
      target_share = target_share,
      realized_share = sum(reference$catch_expected) / sum(reference$deaths_expected),
      F_mean = mean(reference$F_realized), M_mean = mean(reference$M_realized),
      Z_mean = mean(reference$Z), nonpositive_m = sum(reference$M_realized <= 0),
      minimum_m = min(reference$M_realized),
      max_abs_balance = max(abs(reference$Z - reference$M_realized - reference$F_realized)),
      max_abs_removal_error = max(abs(reconstructed - reference$catch_expected)),
      shortfalls = sum(annual$shortfall > 0), shortfall_n = sum(annual$shortfall))
    for (effort in efforts) {
      monitor <- lapply(seq_along(replicates), function(i) {
        choices <- lapply(draws, function(d) d$monitoring[[i]][[effort]])
        data <- tinyAM::as_tam_data(list(catch = gather(lapply(choices, `[[`, "catch")),
          index = gather(lapply(choices, `[[`, "index")), weight = weight, maturity = maturity))
        list(data = data, sampling = gather(lapply(choices, `[[`, "sampling")))
      })
      name <- paste(driver, effort, sep = "-")
      bundles[[name]] <- list(id = id, driver = driver, effort = effort, seeds = seeds,
        reference = reference, data = monitor[[1L]]$data, holdout = monitor[[2L]]$data,
        sampling = rbind(monitor[[1L]]$sampling, monitor[[2L]]$sampling), audit = audit,
        annual_mortality = annual)
    }
  }
  bundles
}

# Preserve survivors and total removals when reducing the assessment age grid.
delta_matched_plus_reference <- function(reference, plus_age = 14L) {
  fields <- c("year", "age", "N", "Z", "F_realized", "M_realized", "ssb", "catch_expected", "index_expected")
  terminal <- subset(reference, age >= plus_age)
  pooled <- do.call(rbind, lapply(split(terminal, terminal$year), function(d) {
    total_n <- sum(d$N)
    total_c <- sum(d$catch_expected)
    total_s <- sum(d$N * exp(-d$Z))
    Z <- if (total_n > 0) -log(total_s / total_n) else mean(d$Z)
    F <- if (total_n > 0 && Z > 0) total_c / total_n * Z / -expm1(-Z) else 0
    data.frame(year = d$year[1L], age = plus_age, N = total_n, Z = Z,
      F_realized = F, M_realized = Z - F, ssb = sum(d$ssb), catch_expected = total_c,
      index_expected = sum(d$index_expected))
  }))
  out <- rbind(reference[reference$age < plus_age, fields], pooled)
  out[order(out$year, out$age), , drop = FALSE]
}
