# Sourced by simulate_observation_error.R after the existing simulation helpers.
# The earlier sampling and population-pilot checkpoints are independent.

oe_z_manifest <- function(phase = c("study", "pilot"), assessments = FALSE) {
  phase <- match.arg(phase)
  stopifnot(is.logical(assessments), length(assessments) == 1L, !is.na(assessments))
  scenarios <- expand.grid(spatial = c("coupled", "default", "iid_age"),
    mortality = c("iid", "ar1", "shared"), stringsAsFactors = FALSE)
  scenarios <- rbind(scenarios, data.frame(spatial = "default", mortality = "temporal_only"))
  scenarios$scenario <- with(scenarios, paste(spatial, mortality, sep = "_"))
  scenarios$phi_age <- unname(c(iid = 0, ar1 = .9, shared = 1, temporal_only = 0)[scenarios$mortality])
  scenarios$phi_year <- ifelse(scenarios$mortality == "temporal_only", .5, 0)
  scenarios$assessment <- assessments & phase == "study" &
    scenarios$scenario %in% c("default_iid", "default_shared", "iid_age_iid")

  # Recruitment streams come from four designated histories in the SD pilot.
  # The classifications below refer to those histories under constant Z.
  histories <- data.frame(
    history_id = sprintf("history%02d", 1:4),
    pilot_population = c("pop07", "pop02", "pop03", "pop04"),
    pilot_regime = c("stable", "decline", "rise_fall", "stable"),
    recruitment_seed = c(10024L, 10002L, 10004L, 10015L),
    mortality_seed = 970000100L + 1:4,
    spatial_seed = 970000200L + 1:4,
    monitoring_seed = 971000000L + 100000L * (1:4),
    assessment = assessments & phase == "study" & 1:4 <= 3L
  )
  if (phase == "pilot") histories <- histories[1L, , drop = FALSE]
  list(scenarios = scenarios, histories = histories,
    settings = list(phase = phase, grid_n = 40L, effort = "standard",
      years = if (phase == "pilot") 5L else c(5L, 16L),
      repeats = if (phase == "pilot") 2L else 20L,
      log_z_marginal_sd = .20, simsurvey_revision = .observation_error_revision))
}

oe_z_demography <- function(history, scenario, log_z_marginal_sd = .20) {
  ages <- years <- seq_len(20L)
  baseline <- .3 + .06 * .oe_gillnet_q(ages)
  recruitment <- .observation_error_with_seed(history$recruitment_seed,
    SimSurvey::sim_R()(years))
  initial <- SimSurvey::sim_N0()(recruitment[1L], baseline, ages)

  # Hold arithmetic E[Z] and marginal log-Z SD fixed across all designs.
  # For a shared age anomaly use one series, avoiding an AR1 endpoint.
  phi <- c(scenario$phi_age, scenario$phi_year)
  log_mean <- log(baseline) - log_z_marginal_sd^2 / 2
  if (scenario$mortality == "shared") {
    common <- .observation_error_with_seed(history$mortality_seed,
      SimSurvey::sim_Z(log_mean = 0, log_sd = log_z_marginal_sd * sqrt(1 - phi[2L]^2),
        phi_age = 0, phi_year = phi[2L])(1L, years))
    Z <- outer(exp(log_mean), drop(common))
  } else {
    innovation_sd <- log_z_marginal_sd * sqrt(prod(1 - phi^2))
    Z <- .observation_error_with_seed(history$mortality_seed,
      SimSurvey::sim_Z(log_mean = outer(log_mean, rep(1, length(years))),
        log_sd = innovation_sd, phi_age = phi[1L], phi_year = phi[2L])(ages, years))
  }
  dimnames(Z) <- list(age = ages, year = years)
  stopifnot(all(is.finite(Z)), all(Z > .06 * .oe_gillnet_q(ages)))
  out <- SimSurvey::sim_abundance(ages, years, Z = function(ages, years) Z,
    R = function(years) recruitment, N0 = function(R0, Z0, ages) initial,
    maturity = SimSurvey::sim_maturity(k = 1, x50 = 5), maturity_by = "age")
  stopifnot(identical(out$Z, Z), all(out$N[, 1L] == initial),
    all(out$R == recruitment),
    isTRUE(all.equal(unname(out$N[-1L, -1L]),
      unname(out$N[-length(ages), -length(years)] * exp(-Z[-length(ages), -length(years)])),
      tolerance = 1e-12)))
  out
}

.oe_z_field <- function(demography, spatial, seed, grid_n) {
  grid <- SimSurvey::make_grid(res = rep(280 / grid_n, 2L))
  field <- switch(spatial,
    coupled = SimSurvey::sim_ays_covar(group_ages = demography$ages),
    default = SimSurvey::sim_ays_covar(),
    iid_age = SimSurvey::sim_ays_covar(sd = 2.8 / sqrt(1 - .5^2),
      phi_age = 0, group_ages = NULL))
  pop <- .observation_error_with_seed(seed,
    SimSurvey::sim_distribution(demography, grid = grid, ays_covar = field))
  prob <- pop$sp_N[, .(cell, year, age,
    probability = N / demography$N[cbind(match(age, demography$ages),
      match(year, demography$years))])]
  stopifnot(max(abs(prob[, sum(probability), by = .(year, age)]$V1 - 1)) < 1e-12)
  list(grid = pop$grid, grid_xy = pop$grid_xy, probability = prob)
}

.oe_z_population <- function(demography, field) {
  counts <- data.table::copy(field$probability)
  N <- demography$N
  counts[, N := .oe_balance(probability * N[.BY$age, .BY$year],
    round(N[.BY$age, .BY$year])), by = .(year, age)]
  counts[, probability := NULL]
  pop <- SimSurvey::round_sim(c(demography,
    list(grid = field$grid, grid_xy = field$grid_xy, sp_N = counts)))
  stopifnot(identical(pop$Z, demography$Z),
    isTRUE(all.equal(pop$N, round(demography$N), tolerance = 0)))
  pop
}

.oe_z_keys <- function(x, scenario, history) {
  if (is.null(x)) return(NULL)
  d <- data.table::copy(data.table::as.data.table(x))
  d[, ':='(scenario = scenario$scenario, spatial = scenario$spatial,
    mortality = scenario$mortality, phi_age = scenario$phi_age,
    phi_year = scenario$phi_year, history_id = history$history_id)]
  d
}

.oe_z_seed <- function(history, program, year, replicate = 0L) {
  offset <- c(Removals = 0L, Gillnet = 30000L, RV = 60000L)
  as.integer(history$monitoring_seed + offset[[program]] + year * 100L + replicate)
}

.oe_z_attempt <- function(attempt, scenario, history, component, year, replicate, seed, path) {
  .oe_z_keys(data.frame(component = component, year = year, replicate = replicate,
    seed = seed, file = path, elapsed = attempt$elapsed, success = is.null(attempt$error),
    warnings = paste(attempt$warnings, collapse = " | "),
    error = if (is.null(attempt$error)) "" else attempt$error), scenario, history)
}

.oe_z_draw <- function(snapshot, removals, component, seed) {
  if (component == "Gillnet") {
    rows <- oe_gillnet_draw(removals, seed, efforts = "standard")
  } else {
    rows <- oe_rv_draw(snapshot, seed, efforts = "standard")$rows
    stopifnot(all(rows$n_primary == 160L))
  }
  stopifnot(nrow(rows) == 34L, !anyDuplicated(rows[, .(age, bin)]),
    all(rows$component == component), all(rows$year == snapshot$years),
    all(rows$effort == "standard"), all(rows$truth >= 0),
    all(is.na(rows$obs) | (is.finite(rows$obs) & rows$obs >= 0)),
    setequal(rows$age[rows$bin == "modeled"], 1:14))
  rows
}

.oe_z_truth_tables <- function(pop) {
  p <- SimSurvey::stock_biology(pop)$at_age[, .(year, age, weight = mean_weight_kg,
    maturity = prop_mature, N = as.numeric(pop$N) / 1000, Z = as.numeric(pop$Z),
    SSB = as.numeric(pop$S) / 1000)]
  list(population = p, annual = p[, .(SSB = sum(SSB), recruitment = N[age == 1L]), by = year])
}

.oe_z_case <- function(pop, scenario, history, settings, cache, assessment_path) {
  assessment <- !is.null(assessment_path)
  years <- if (assessment) pop$years else settings$years
  keep <- seq_along(years)
  streams <- c(Gillnet = "catch", RV = "index")
  draws <- attempts <- removals_log <- shortfalls <- realized <- vector("list", length(years))
  training <- validation <- stats::setNames(vector("list", length(streams)), streams)
  overlap <- oe_overlap(pop)
  truth <- .oe_z_truth_tables(pop)
  for (i in keep) {
    year <- years[i]
    admit()
    snapshot <- oe_snapshot(pop, year)
    seed <- .oe_z_seed(history, "Removals", year)
    path <- file.path(cache, sprintf("removals-%02d.rds", year))
    removal <- checkpoint(path, function() .oe_pilot_removals(snapshot, seed))
    removals_log[[i]] <- .oe_z_attempt(removal, scenario, history, "Removals", year, 0L, seed, path)
    if (settings$phase == "pilot") required(removal)
    if (is.null(removal$error)) {
      shortfalls[[i]] <- removal$value$fishery_annual
      realized[[i]] <- removal$value$population[, .(year, age, F, M)]
    }
    sampling <- year %in% settings$years
    repeats <- if (sampling) seq_len(settings$repeats) else 1:2
    records <- list()
    for (replicate in repeats) {
      for (component in names(streams)) {
        seed <- .oe_z_seed(history, component, year, replicate)
        path <- file.path(cache, sprintf("%s-%02d-%02d.rds", component, year, replicate))
        result <- checkpoint(path, function() {
          if (component == "Gillnet") required(removal)
          .oe_z_draw(snapshot, removal$value, component, seed)
        })
        records[[length(records) + 1L]] <- .oe_z_attempt(result, scenario,
          history, component, year, replicate, seed, path)
        if (settings$phase == "pilot") required(result)
        if (!is.null(result$error)) next
        if (sampling) {
          rows <- .oe_z_keys(result$value, scenario, history)
          rows[, replicate := replicate]
          draws[[i]][[length(draws[[i]]) + 1L]] <- rows
        }
        if (assessment && replicate %in% 1:2) {
          series <- streams[[component]]
          rows <- as.data.frame(result$value[result$value$bin == "modeled", ])
          if (replicate == 1L) training[[series]][[i]] <- rows else validation[[series]][[i]] <- rows
        }
      }
    }
    attempts[[i]] <- data.table::rbindlist(records)
    rm(snapshot, removal)
    invisible(gc())
  }
  realized <- data.table::rbindlist(realized)
  p <- truth$population
  p[, c("F", "M") := list(NA_real_, NA_real_)]
  if (nrow(realized)) p[realized, on = .(year, age), c("F", "M") := list(i.F, i.M)]
  assessment_input <- NULL
  if (assessment) {
    bind <- function(x) as.data.frame(data.table::rbindlist(x))
    complete <- function(x) length(x) == length(years) && all(lengths(x) > 0L)
    expected <- all(vapply(c(training, validation), complete, logical(1))) &&
      nrow(realized) == nrow(truth$population)
    if (expected && !file.exists(assessment_path)) {
      dir.create(dirname(assessment_path), recursive = TRUE, showWarnings = FALSE)
      input <- list(scenario = scenario, history = history, training = lapply(training, bind),
        validation = lapply(validation, bind), truth = lapply(truth, as.data.frame))
      save_checkpoint(input, assessment_path)
    }
    assessment_input <- .oe_z_keys(data.frame(file = assessment_path,
      complete = expected, stringsAsFactors = FALSE), scenario, history)
  }
  list(rows = data.table::rbindlist(unlist(draws, recursive = FALSE), fill = TRUE),
    draw_attempts = data.table::rbindlist(attempts),
    removal_attempts = data.table::rbindlist(removals_log),
    population = .oe_z_keys(p, scenario, history),
    annual = .oe_z_keys(truth$annual, scenario, history),
    overlap = .oe_z_keys(overlap$summary, scenario, history),
    overlap_maps = .oe_z_keys(overlap$map, scenario, history),
    target_shortfalls = .oe_z_keys(data.table::rbindlist(shortfalls), scenario, history),
    assessment_inputs = assessment_input)
}

oe_z_simulate <- function(output, phase = c("pilot", "study"), assessments = FALSE) {
  phase <- match.arg(phase)
  .observation_error_check_revision()
  manifest <- oe_z_manifest(phase, assessments)
  source_dir <- file.path("analysis", "formula_expansion")
  manifest$sources <- tools::md5sum(file.path(source_dir, c("observation_error_z_sim.R",
    "observation_error_sim.R", "fit_observation_error_pilot.R")))
  simsurvey_source <- Sys.getenv("TINYAM_OE_SIMSURVEY")
  if (nzchar(simsurvey_source)) {
    stopifnot(identical(basename(simsurvey_source), paste0("SimSurvey-", .observation_error_revision)))
    simsurvey_files <- list.files(file.path(simsurvey_source, "R"), "[.]R$", full.names = TRUE)
    stopifnot(length(simsurvey_files) > 0L)
    manifest$simsurvey_sources <- stats::setNames(unname(tools::md5sum(simsurvey_files)), basename(simsurvey_files))
  }
  manifest$versions <- c(R = as.character(getRversion()),
    vapply(c("SimSurvey", "data.table"), function(pkg) {
      as.character(utils::packageVersion(pkg))
    }, character(1)))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  verify_manifest(manifest, file.path(output, "manifest.rds"))
  settings <- manifest$settings
  cases <- list()
  for (h in seq_len(nrow(manifest$histories))) {
    history <- manifest$histories[h, , drop = FALSE]
    for (spatial in unique(manifest$scenarios$spatial)) {
      scenarios <- manifest$scenarios[manifest$scenarios$spatial == spatial, , drop = FALSE]
      field <- NULL
      for (s in seq_len(nrow(scenarios))) {
        scenario <- scenarios[s, , drop = FALSE]
        cache <- file.path(output, "checkpoints", scenario$scenario, history$history_id)
        dir.create(cache, recursive = TRUE, showWarnings = FALSE)
        path <- file.path(cache, "case.rds")
        assessment_path <- if (scenario$assessment && history$assessment) {
          file.path(output, "assessment", scenario$scenario, history$history_id, "input.rds")
        }
        if (file.exists(path)) {
          result <- readRDS(path)
        } else {
          admit()
          demography <- oe_z_demography(history, scenario, settings$log_z_marginal_sd)
          if (is.null(field)) {
            field <- .oe_z_field(demography, spatial, history$spatial_seed, settings$grid_n)
          }
          pop <- .oe_z_population(demography, field)
          result <- .oe_z_case(pop, scenario, history, settings, cache, assessment_path)
          save_checkpoint(result, path)
          rm(pop, demography)
          invisible(gc())
        }
        cases[[length(cases) + 1L]] <- result
      }
    }
  }
  keys <- names(cases[[1L]])
  summary <- stats::setNames(lapply(keys, function(key) {
    data.table::rbindlist(lapply(cases, "[[", key), fill = TRUE)
  }), keys)
  summary$mortality <- summary$population[, .(scenario, spatial, mortality, history_id,
    phi_age, phi_year, year, age, Z)]
  summary$manifest <- manifest
  summary
}
