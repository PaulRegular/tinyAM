# Generate the compact, 3Ps-inspired teaching data; never run during installation.
# Rscript data-raw/sim_cod_obs.R prepare /path/to/SimSurvey
# Rscript data-raw/validate_sim_cod.R
# SimSurvey source revision is deliberately fixed. No automatic install/update.
#
# Population settings reproduce Regular et al. (2020), S1/S2 supplements. The
# rigid grid rotation changes orientation only. DFO Data Report 1381 (2024)
# documents predominantly inshore fixed-gear and offshore mobile-gear landings,
# later coastal RV strata, and RV exclusion of strata deeper than 300 fathoms.
# DFO Res. Doc. 2026/056 describes coastal Sentinel sampling and GN catches mainly
# at 50-80 cm, versus LT catches extending to roughly 30 cm. These sources guide
# footprints and contrasting shapes, not exact selectivity midpoint estimates.
#
# Shares, shape parameters, geographic/depth cutoffs, encounter coefficients,
# capacities and sample caps below are illustrative design choices. Operations
# represent aggregate fishing days/hauls, not calibrated literal survey effort.
# No contemporary critical-zone TAC or closures are imposed. All fleets have
# positive allocations in all years; access can still cause recorded shortfalls.
# For fleet f/year t, A=sum(N_region*s*(1-exp(-Z))) and
# G=sum(N_total*s*(1-exp(-Z))). Effective fishery_share=base_share*A/G makes
# the reference annual catch base_share*A. This allocates existing deaths using
# operating-model access; it is neither an assessment predictor nor F feedback
# into future populations. Targeting can still concentrate local/age removals.
# Stock/catch weights and age maturity are assumed known teaching inputs. The
# weight relationship is the 3Ps reference allometry in Data Report 1381 Table10;
# maturity(age) is deliberately illustrative, not the default female length curve.

sim_cod_settings <- function() {
  list(
    simsurvey_revision = "282e3c14a2fbb5d312a2fde79ee3e8d8511b4513",
    seed = 438L, ages = 1:20, years = 1:20,
    recruitment = list(log_mean = log(30e6), log_sd = 0.5, random_walk = TRUE),
    mortality = list(log_mean = log(0.5), log_sd = 0.2, phi_age = 0.9, phi_year = 0.5),
    growth = list(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1, length_group = 3),
    grid = list(x_range = c(-140, 140), y_range = c(-140, 140), res = c(3.5, 3.5),
      shelf_depth = 200, shelf_width = 100, depth_range = c(0, 1000), n_div = 1,
      strat_breaks = seq(0, 1000, by = 40), strat_splits = 2),
    spatial = list(sd = 2.8, range = 300, lambda = 1, model = "matern",
      phi_age = 0.5, phi_year = 0.9, group_ages = 5:20),
    depth = list(mu = 200, sigma = 70),
    weight = list(alpha = 10^-5.2106, beta = 3.0879),
    maturity = list(k = 1, x50 = 5),
    shares = c(Gillnet = 0.30, Longline = 0.15, Trawl = 0.15),
    monitoring = list(sampling_coverage = 0.30, lengths_cap = 150, ages_cap = 3,
      age_sampling = "stratified", age_length_group = 3, sampling_targeting = 0,
      reporting_factor = 1, reporting_cv = 0.05),
    survey_density = c(RV = 2 / 1000, Sentinel = 6 / 1000),
    survey_sampling = list(lengths_cap = 100, ages_cap = 3,
      age_sampling = "stratified", age_length_group = 3, age_space_group = "set"),
    max_operations = 5000L,
    sources = c(
      population = "https://doi.org/10.1371/journal.pone.0232822",
      data_review = "https://publications.gc.ca/collections/collection_2024/mpo-dfo/Fs97-13-1381-eng.pdf",
      management = "https://www.dfo-mpo.gc.ca/fisheries-peches/ifmp-gmp/groundfish-poisson-fond/2025/groundfish-poisson-fond-div3ps-nafo-opano-eng.html",
      sentinel = "https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-056-eng.pdf"))
}

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
      s <- samples[samples$operation == op$operation, , drop = FALSE]
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
  original <- do.call(SimSurvey::make_grid, settings$grid)
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

.sim_cod_selectivity <- function() {
  list(
    Gillnet = function(length) stats::plogis((length - 50) / 5) * stats::plogis((80 - length) / 5),
    Longline = function(length) stats::plogis((length - 40) / 5),
    Trawl = function(length) stats::plogis((length - 45) / 5),
    RV = function(length) stats::plogis((length - 30) / 5),
    Sentinel = function(length) stats::plogis((length - 50) / 5) * stats::plogis((80 - length) / 5))
}

.sim_cod_fisheries <- function(settings, selectivity, allocations) {
  # Operations are aggregate fishing-day/haul units. Capacity and encounter
  # coefficients are illustrative; no claim of calibrated net/hook saturation.
  common <- c(list(selectivity_by = "length", targeting = 1, exploration = 0.05,
    knowledge_sd = 0.2, max_operations = settings$max_operations, light = FALSE,
    maturity_sampling = "none"), settings$monitoring)
  shares <- lapply(names(settings$shares), function(fleet) {
    rows <- allocations[allocations$fleet == fleet, ]
    stats::setNames(rows$fishery_share[match(settings$years, rows$year)], settings$years)
  })
  names(shares) <- names(settings$shares)
  gillnet <- do.call(SimSurvey::sim_fishery, c(list(fishery_share = shares$Gillnet,
    catch_model = SimSurvey::sim_capture(q = 0.02, capacity_per_gear = 150),
    selectivity = selectivity$Gillnet, gear_amount = function(operation) {
      if (operation$fishery_operation == 1L) message("Gillnet: year ", operation$year)
      60
    }, duration = 24,
    gear_unit = "net", time_unit = "hour", subset_cells = quote(Gillnet)), common))
  longline <- do.call(SimSurvey::sim_fishery, c(list(fishery_share = shares$Longline,
    catch_model = SimSurvey::sim_capture(q = 0.001, capacity_per_gear = 1),
    selectivity = selectivity$Longline, gear_amount = function(operation) {
      if (operation$fishery_operation == 1L) message("Longline: year ", operation$year)
      4000
    }, duration = 24,
    gear_unit = "hook", time_unit = "hour", subset_cells = quote(Longline)), common))
  trawl <- do.call(SimSurvey::sim_fishery, c(list(fishery_share = shares$Trawl,
    catch_model = SimSurvey::sim_capture(q = 2, capacity_per_gear = 10000),
    selectivity = selectivity$Trawl, gear_amount = function(operation) {
      if (operation$fishery_operation == 1L) message("Trawl: year ", operation$year)
      1
    }, duration = 4,
    gear_unit = "trawl", time_unit = "hour", subset_cells = quote(Trawl)), common))
  list(Gillnet = gillnet, Longline = longline, Trawl = trawl)
}

.sim_cod_population <- function(settings) {
  set.seed(settings$seed)
  abundance <- SimSurvey::sim_abundance(ages = settings$ages, years = settings$years,
    R = do.call(SimSurvey::sim_R, settings$recruitment),
    Z = do.call(SimSurvey::sim_Z, settings$mortality), N0 = SimSurvey::sim_N0("exp"),
    growth = do.call(SimSurvey::sim_growth, settings$growth),
    weight = do.call(SimSurvey::sim_weight, settings$weight),
    maturity = do.call(SimSurvey::sim_maturity, settings$maturity), maturity_by = "age")
  distributed <- SimSurvey::sim_distribution(abundance, grid = .sim_cod_grid(settings),
    ays_covar = do.call(SimSurvey::sim_ays_covar, settings$spatial),
    depth_par = do.call(SimSurvey::sim_parabola, settings$depth))
  list(value = SimSurvey::round_sim(distributed), rng = .Random.seed)
}

.sim_cod_surveys <- function(simulation, settings, selectivity) {
  configs <- lapply(names(settings$survey_density), function(label) {
    # Pass the literal grid-column name: sim_sets evaluates subset_cells inside
    # its grid table, without retaining this calling frame's local variables.
    sets <- do.call(SimSurvey::sim_sets, list(sim = simulation,
      subset_cells = as.name(label), n_sims = 1,
      set_den = settings$survey_density[[label]], min_sets = 2,
      trawl_dim = c(1.5, 0.02)))
    do.call(SimSurvey::sim_survey, c(list(q = selectivity[[label]], q_by = "length",
      custom_sets = sets, maturity_sampling = "none", light = FALSE), settings$survey_sampling))
  })
  names(configs) <- names(settings$survey_density)
  do.call(SimSurvey::sim_surveys, c(list(sim = simulation, n_sims = 1), configs))
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

.sim_cod_tables <- function(simulation, settings, selectivity, available, provenance, allocations) {
  ages <- settings$ages
  years <- settings$years
  fleets <- names(settings$shares)
  coordinates <- expand.grid(age = ages, year = years)
  raised <- raise_sim_cod_catch(simulation$fishery$setdet, simulation$fishery$samp,
    ages, years, fleets, settings$monitoring$age_length_group)
  catch <- raised$catch
  catch$obs <- catch$obs / 1000
  catch$unit <- "thousand fish"
  index <- survey_precision <- survey_monitoring <- list()
  for (survey in names(settings$survey_density)) {
    selected <- SimSurvey::select_survey(simulation, survey)
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
  key <- simulation$sim_length(age = ages, length_age_key = TRUE)
  lengths <- as.numeric(rownames(key))
  weights <- simulation$sim_weight(lengths)
  mean_weight <- as.numeric(crossprod(weights, key))
  maturity_at_age <- simulation$sim_maturity(ages)
  stock_weight <- data.frame(year = coordinates$year, age = coordinates$age,
    obs = mean_weight[match(coordinates$age, ages)], type = "stock", fleet = NA_character_, unit = "kg")
  catches <- as.data.frame(simulation$fishery$catch_at_age)
  catch_weight <- list()
  fleet_truth <- as.data.frame(simulation$fishery$reference)
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
  dat <- tinyAM::as_tam_data(list(catch = catch, index = index,
    weight = rbind(stock_weight, do.call(rbind, catch_weight)), maturity = maturity))
  tinyAM::check_tam_data(dat, quiet = TRUE)
  total_f <- stats::aggregate(fleet_truth$F_realized, fleet_truth[c("year", "age")], sum)
  names(total_f)[3] <- "F"
  population <- coordinates
  population$N <- as.numeric(simulation$N) / 1000
  population$Z <- as.numeric(simulation$Z)
  population$weight <- mean_weight[match(population$age, ages)]
  population$maturity <- maturity_at_age[match(population$age, ages)]
  population$B <- population$N * population$weight
  population$SSB <- population$B * population$maturity
  population$F <- total_f$F[match(paste(population$year, population$age), paste(total_f$year, total_f$age))]
  population$M <- population$Z - population$F
  annual <- stats::aggregate(population[c("N", "B", "SSB")], population["year"], sum)
  annual$recruitment <- population$N[match(paste(annual$year, min(ages)), paste(population$year, population$age))]
  realized_yield <- stats::aggregate(fleet_truth$yield, fleet_truth["year"], sum)
  annual$yield <- realized_yield$x[match(annual$year, realized_yield$year)]
  available$N_region <- available$N_region / 1000
  available$N_total <- available$N_total / 1000
  available$index_expected <- available$index_expected / 1000
  regional <- available[available$program %in% names(settings$survey_density), ]
  names(regional)[names(regional) == "program"] <- "survey"
  footprint <- as.data.frame(simulation$grid)
  footprint$area_km2 <- abs(prod(stars::st_res(simulation$grid)))
  curves <- do.call(rbind, lapply(names(selectivity), function(label) {
    data.frame(component = label, length = seq(0, 140), selectivity = selectivity[[label]](seq(0, 140)))
  }))
  for (column in intersect(c("N", "catch_reference", "catch_realized", "total_deaths", "death_budget"), names(fleet_truth))) {
    fleet_truth[[column]] <- fleet_truth[[column]] / 1000
  }
  loss <- population$N[population$age == max(ages)] * exp(-population$Z[population$age == max(ages)])
  checks <- list(terminal_survivors = data.frame(year = years, N = loss,
      proportion = loss / annual$N),
    negative_remaining_M = population[is.finite(population$M) & population$M < 0, c("year", "age", "M")],
    shortfalls = as.data.frame(simulation$fishery$annual)[c("fishery", "year", "target", "n", "shortfall", "status")],
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
    fishery_annual = as.data.frame(simulation$fishery$annual),
    monitoring = raised$monitoring, survey_precision = do.call(rbind, survey_precision),
    survey_monitoring = do.call(rbind, survey_monitoring),
    allocations = allocations,
    checks = checks, provenance = provenance)
  list(sim_cod_obs = dat, sim_cod_truth = truth)
}

build_sim_cod_obs <- function(simsurvey_repo, repo = ".", cache_dir = file.path(repo, "data-raw", "cache", "sim_cod"),
                              overwrite = FALSE) {
  settings <- sim_cod_settings()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  if (missing(simsurvey_repo) || !nzchar(simsurvey_repo)) stop("Supply a clean SimSurvey checkout at the recorded revision; see data-raw/sim_cod_obs.R.")
  revision <- .sim_cod_git(simsurvey_repo)
  if (!identical(revision, settings$simsurvey_revision)) stop("SimSurvey revision differs from the recorded generation source.")
  if (length(system2("git", c("-C", shQuote(simsurvey_repo), "status", "--porcelain", "--untracked-files=no"), stdout = TRUE))) stop("Use a clean SimSurvey source checkout.")
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload to load the recorded SimSurvey source.")
  pkgload::load_all(repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  pkgload::load_all(simsurvey_repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  records <- list()
  cached <- function(name, key, code) {
    path <- file.path(cache_dir, paste0(name, ".rds"))
    if (!overwrite && file.exists(path)) {
      old <- readRDS(path)
      if (identical(old$key, key)) {
        records[[name]] <<- old[c("warnings", "error", "elapsed_seconds")]
        if (!is.null(old$error)) stop(name, ": previous generation failed; inspect ", path,
          " and rerun with overwrite = TRUE after resolving the cause. ", old$error)
        message(name, ": reused cached result")
        return(old$value)
      }
    }
    warnings <- character()
    message(name, ": generating (warnings and failures retained)")
    start <- proc.time()[["elapsed"]]
    error <- NULL
    result <- tryCatch(withCallingHandlers(code(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    record <- list(key = key, value = result, warnings = warnings, error = error,
      elapsed_seconds = proc.time()[["elapsed"]] - start)
    saveRDS(record, path)
    records[[name]] <<- record[c("warnings", "error", "elapsed_seconds")]
    if (!is.null(error)) stop(name, ": ", error)
    if (length(warnings)) message(name, ": retained ", length(warnings), " warning(s)")
    result
  }
  runtime <- vapply(c("SimSurvey", "sf", "stars", "data.table"), function(x) as.character(utils::packageVersion(x)), character(1))
  source_files <- list.files(file.path(simsurvey_repo, "R"), "[.]R$", full.names = TRUE)
  provenance <- list(simsurvey_revision = revision, seed = settings$seed, settings = settings,
    rng_kind = RNGkind(),
    sources = stats::setNames(unname(tools::md5sum(source_files)), basename(source_files)),
    generator_md5 = unname(tools::md5sum(file.path(repo, "data-raw", "sim_cod_obs.R"))),
    runtime = c(R = as.character(getRversion()), runtime),
    units = c(abundance = "thousand fish", weight = "kg", biomass = "tonnes", mortality = "per year"),
    generated = Sys.time())
  population_key <- .sim_cod_hash(list(settings[c("seed", "ages", "years", "recruitment",
    "mortality", "growth", "grid", "spatial", "depth", "weight", "maturity")],
    body(.sim_cod_population), body(.sim_cod_grid), provenance$sources, provenance$runtime))
  population <- cached("population", population_key, function() .sim_cod_population(settings))
  selectivity <- .sim_cod_selectivity()
  available <- .sim_cod_available(population$value, settings, selectivity)
  allocations <- .sim_cod_allocations(population$value, settings, available)
  provenance$allocation <- list(method = "base_share * accessible/global vulnerable expected deaths",
    table_hash = .sim_cod_hash(allocations))
  fishery_key <- .sim_cod_hash(list(population_key, allocations, settings$monitoring,
    settings$max_operations, body(.sim_cod_allocations), body(.sim_cod_fisheries), body(.sim_cod_selectivity)))
  survey_key <- .sim_cod_hash(list(fishery_key, settings$survey_density, settings$survey_sampling,
    body(.sim_cod_surveys)))
  fisheries <- cached("fisheries", fishery_key, function() {
    assign(".Random.seed", population$rng, envir = .GlobalEnv)
    configs <- .sim_cod_fisheries(settings, selectivity, allocations)
    result <- do.call(SimSurvey::sim_fisheries, c(list(sim = population$value, n_sims = 1), configs))
    list(value = result, rng = .Random.seed)
  })
  surveys <- cached("surveys", survey_key, function() {
    assign(".Random.seed", fisheries$rng, envir = .GlobalEnv)
    .sim_cod_surveys(fisheries$value, settings, selectivity)
  })
  provenance$records <- records
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  error <- NULL
  result <- tryCatch(withCallingHandlers(
    .sim_cod_tables(surveys, settings, selectivity, available, provenance, allocations),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
  provenance$records$tables <- list(warnings = warnings, error = error,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  if (!is.null(error)) {
    saveRDS(list(provenance = provenance), file.path(cache_dir, "generation-report.rds"))
    stop("Table preparation failed; see the retained generation report. ", error)
  }
  result$sim_cod_truth$provenance <- provenance
  dir.create(file.path(repo, "data"), showWarnings = FALSE)
  sim_cod_obs <- result$sim_cod_obs
  sim_cod_truth <- result$sim_cod_truth
  save(sim_cod_obs, file = file.path(repo, "data", "sim_cod_obs.rda"), compress = "xz", version = 2)
  save(sim_cod_truth, file = file.path(repo, "data", "sim_cod_truth.rda"), compress = "xz", version = 2)
  saveRDS(list(provenance = provenance, checks = sim_cod_truth$checks,
    monitoring = sim_cod_truth$monitoring), file.path(cache_dir, "generation-report.rds"))
  message("Generated sim_cod_obs and sim_cod_truth; inspect retained warnings, shortfalls and sampling coverage before fitting.")
  invisible(result)
}

sim_cod_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  mode <- if (length(args)) args[[1]] else "prepare"
  if (!identical(mode, "prepare")) stop("Use prepare /path/to/SimSurvey; source the script for pure helpers.")
  source <- if (length(args) > 1L) args[[2]] else Sys.getenv("SIMSURVEY_REPO")
  build_sim_cod_obs(source)
}

if (sys.nframe() == 0L) sim_cod_main()
