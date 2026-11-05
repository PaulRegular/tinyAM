# Reproduce or explore the 3Ps-inspired teaching scenario. SimSurvey is used here
# only, not by tinyAM installation or documentation builds.
#
# source("data-raw/sim_cod_obs.R")
# result <- build_sim_cod_obs("../SimSurvey")  # selected 20 x 20 teaching case
# Trial results and three stage caches use configuration-specific folders.
# write_data = TRUE is an explicit decision to replace the packaged .rda files.
# Existing packaged data and earlier simulation caches are never changed by default.
#
# Population centres follow Regular et al. (2020), S1/S2; smaller mortality and
# length variation simplify this teaching case. Coastal fixed gear,
# offshore trawling/RV and Sentinel size coverage are inspired by DFO Data Report
# 1381 (2024) and Res. Doc. 2026/056. Footprints, mortality shares, capacities and
# sample caps are illustrative choices, not a reproduction of the assessment.
# All fleets operate; no contemporary rebuilding TAC or closures are imposed.
# The RV curve is logistic and the gillnet/Sentinel curves are double logistic.
# Operation-level simulation can require several GB of RAM; generate it on its own.
# Optional reference_M calibrates the geometric centre of an M baseline using
# exp(log_mean_Z), the Z process median/log-scale centre. It is not E[Z] or a
# target for arithmetic realized M. NULL skips this calibration in explorations.

sim_cod_settings <- function(res = c(14, 14)) {
  list(
    simsurvey_revision = "282e3c14a2fbb5d312a2fde79ee3e8d8511b4513",
    seed = 438L, ages = 1:20, years = 1:20,
    recruitment = list(log_mean = log(30e6), log_sd = 0.5, random_walk = TRUE),
    mortality = list(log_mean = log(0.5), log_sd = 0.1, phi_age = 0.9, phi_year = 0.5),
    growth = list(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.03, length_group = 3),
    grid = list(x_range = c(-140, 140), y_range = c(-140, 140), res = res,
      shelf_depth = 200, shelf_width = 100, depth_range = c(0, 1000), n_div = 1,
      strat_breaks = seq(0, 1000, by = 40), strat_splits = 2),
    spatial = list(sd = 2.8, range = 300, lambda = 1, model = "matern",
      phi_age = 0, phi_year = 0, group_ages = 1:20, group_years = 1:20),
    depth = list(mu = 200, sigma = 70),
    weight = list(alpha = 10^-5.2106, beta = 3.0879),
    maturity = list(k = 1, x50 = 5),
    selectivity = list(Gillnet = c(up = 50, width = 5, down = 80, down_width = 5),
      Longline = c(up = 40, width = 5), Trawl = c(up = 45, width = 5),
      RV = c(up = 30, width = 5), Sentinel = c(up = 50, width = 5, down = 80, down_width = 5)),
    shares = c(Gillnet = 0.30, Longline = 0.15, Trawl = 0.15),
    reference_M = 0.3,
    monitoring = list(sampling_coverage = 0.30, lengths_cap = 150, ages_cap = 3,
      age_sampling = "stratified", age_length_group = 3, sampling_targeting = 0,
      reporting_factor = 1, reporting_cv = 0.05),
    survey_density = c(RV = 2 / 1000, Sentinel = 4 / 1000),
    survey_sampling = list(lengths_cap = 100, ages_cap = 3,
      age_sampling = "stratified", age_length_group = 3, age_space_group = "set"),
    max_operations = 5000L,
    sources = c(
      population = "https://doi.org/10.1371/journal.pone.0232822",
      data_review = "https://publications.gc.ca/collections/collection_2024/mpo-dfo/Fs97-13-1381-eng.pdf",
      management = "https://www.dfo-mpo.gc.ca/fisheries-peches/ifmp-gmp/groundfish-poisson-fond/2025/groundfish-poisson-fond-div3ps-nafo-opano-eng.html",
      sentinel = "https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-056-eng.pdf"))
}

build_sim_cod_obs <- function(simsurvey_repo, settings = sim_cod_settings(),
                              repo = ".", cache_dir = NULL, write_data = FALSE,
                              overwrite = FALSE) {
  helpers <- new.env(parent = globalenv())
  sys.source(file.path(repo, "data-raw", "sim_cod_helpers.R"), envir = helpers)
  if (!identical(helpers$.sim_cod_git(simsurvey_repo), settings$simsurvey_revision)) {
    stop("Use the SimSurvey revision recorded in sim_cod_settings().")
  }
  if (length(system2("git", c("-C", shQuote(simsurvey_repo), "status", "--porcelain",
      "--untracked-files=all", "--", "R", "DESCRIPTION", "NAMESPACE"), stdout = TRUE))) stop("Use clean SimSurvey package sources.")
  pkgload::load_all(repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  pkgload::load_all(simsurvey_repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  source_files <- file.path(repo, "data-raw", c("sim_cod_obs.R", "sim_cod_helpers.R"))
  provenance <- list(simsurvey_revision = settings$simsurvey_revision,
    seed = settings$seed, rng_kind = RNGkind(), settings = settings,
    sources = stats::setNames(unname(tools::md5sum(source_files)), basename(source_files)),
    runtime = c(R = as.character(getRversion()), vapply(c("SimSurvey", "sf", "stars", "data.table"),
      function(package) as.character(utils::packageVersion(package)), character(1))),
    units = c(abundance = "thousand fish", weight = "kg", biomass = "tonnes", mortality = "per year"))
  # One readable configuration/source fingerprint, not a dependency graph of
  # reflected function bodies. Different trials cannot share a cache folder.
  fingerprint <- helpers$.sim_cod_hash(provenance)
  if (is.null(cache_dir)) cache_dir <- file.path(repo, "data-raw", "cache", "sim_cod", fingerprint)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  records <- list()
  stage <- function(name, code, cache = TRUE) {
    path <- file.path(cache_dir, paste0(name, ".rds"))
    if (cache && !overwrite && file.exists(path)) {
      saved <- readRDS(path)
      if (!identical(saved$fingerprint, fingerprint)) stop("The selected cache belongs to a different configuration.")
      if (!is.null(saved$error)) stop("The previous ", name, " attempt failed: ", saved$error)
      assign(".Random.seed", saved$rng, envir = .GlobalEnv)
      records[[name]] <<- saved[c("warnings", "error", "elapsed_seconds")]
      message(name, ": reused cache")
      return(saved$value)
    }
    message(name, ": generating")
    warnings <- character(); error <- NULL; started <- proc.time()[["elapsed"]]
    value <- tryCatch(withCallingHandlers(force(code), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    saved <- list(fingerprint = fingerprint, value = value, rng = .Random.seed,
      warnings = warnings, error = error, elapsed_seconds = proc.time()[["elapsed"]] - started)
    if (cache) saveRDS(saved, path)
    records[[name]] <<- saved[c("warnings", "error", "elapsed_seconds")]
    if (!is.null(error)) {
      saveRDS(list(provenance = provenance, records = records), file.path(cache_dir, "generation-report.rds"))
      stop(name, ": ", error)
    }
    value
  }

  # 1. One population; the original paper settings are visible in these calls.
  population <- stage("population", {
    set.seed(settings$seed)
    SimSurvey::sim_abundance(
      ages = settings$ages, years = settings$years,
      R = SimSurvey::sim_R(log_mean = settings$recruitment$log_mean,
        log_sd = settings$recruitment$log_sd, random_walk = settings$recruitment$random_walk),
      Z = SimSurvey::sim_Z(log_mean = settings$mortality$log_mean,
        log_sd = settings$mortality$log_sd, phi_age = settings$mortality$phi_age,
        phi_year = settings$mortality$phi_year),
      N0 = SimSurvey::sim_N0("exp"),
      growth = SimSurvey::sim_growth(Linf = settings$growth$Linf, L0 = settings$growth$L0,
        K = settings$growth$K, log_sd = settings$growth$log_sd, length_group = settings$growth$length_group),
      weight = SimSurvey::sim_weight(alpha = settings$weight$alpha, beta = settings$weight$beta),
      maturity = SimSurvey::sim_maturity(k = settings$maturity$k, x50 = settings$maturity$x50),
      maturity_by = "age") |>
      SimSurvey::sim_distribution(grid = helpers$.sim_cod_grid(settings),
        ays_covar = SimSurvey::sim_ays_covar(sd = settings$spatial$sd, range = settings$spatial$range,
          lambda = settings$spatial$lambda, model = settings$spatial$model,
          phi_age = settings$spatial$phi_age, phi_year = settings$spatial$phi_year,
          group_ages = settings$spatial$group_ages, group_years = settings$spatial$group_years),
        depth_par = SimSurvey::sim_parabola(mu = settings$depth$mu, sigma = settings$depth$sigma)) |>
      SimSurvey::round_sim()
  })
  # Survey cells are sampled without replacement. Reject infeasible set counts
  # before the expensive fishery stage, rather than failing after it completes.
  helpers$.sim_cod_check_surveys(as.data.frame(population$grid), settings,
    abs(prod(stars::st_res(population$grid))))
  selectivity <- helpers$.sim_cod_selectivity(settings)
  available <- helpers$.sim_cod_available(population, settings, selectivity)
  allocations <- helpers$.sim_cod_allocations(population, settings, available)
  calibration <- NULL
  if (!is.null(settings$reference_M)) {
    calibrated <- helpers$.sim_cod_reference_M(settings, available, allocations)
    allocations <- calibrated$allocations
    calibration <- calibrated$calibration
    provenance$reference_M <- calibration
  }
  shares <- split(allocations$fishery_share, allocations$fleet)
  shares <- lapply(shares, stats::setNames, settings$years)
  m <- settings$monitoring

  # 2. Fishery branch: allocate existing deaths, then sample actual catches.
  # Access-adjusted share = base share * accessible/global vulnerable deaths.
  # SimSurvey retains a configuration's calling environment. Keep this small,
  # so configurations cannot retain the whole population or build workspace.
  program_settings <- list2env(list(settings = settings, selectivity = selectivity,
    shares = shares, m = m), parent = globalenv())
  fisheries <- stage("fisheries", {
    programs <- evalq(list(
      Gillnet = SimSurvey::sim_fishery(fishery_share = shares$Gillnet,
        catch_model = SimSurvey::sim_capture(q = 0.02, capacity_per_gear = 150),
        selectivity = selectivity$Gillnet, selectivity_by = "length", subset_cells = Gillnet,
        gear_amount = 60, duration = 24, gear_unit = "net", time_unit = "hour",
        targeting = 1, knowledge_sd = 0.2, exploration = 0.05, max_operations = settings$max_operations,
        sampling_coverage = m$sampling_coverage, sampling_targeting = m$sampling_targeting,
        lengths_cap = m$lengths_cap, ages_cap = m$ages_cap, age_sampling = m$age_sampling,
        age_length_group = m$age_length_group, reporting_factor = m$reporting_factor,
        reporting_cv = m$reporting_cv, maturity_sampling = "none", light = FALSE),
      Longline = SimSurvey::sim_fishery(fishery_share = shares$Longline,
        catch_model = SimSurvey::sim_capture(q = 0.001, capacity_per_gear = 1),
        selectivity = selectivity$Longline, selectivity_by = "length", subset_cells = Longline,
        gear_amount = 4000, duration = 24, gear_unit = "hook", time_unit = "hour",
        targeting = 1, knowledge_sd = 0.2, exploration = 0.05, max_operations = settings$max_operations,
        sampling_coverage = m$sampling_coverage, sampling_targeting = m$sampling_targeting,
        lengths_cap = m$lengths_cap, ages_cap = m$ages_cap, age_sampling = m$age_sampling,
        age_length_group = m$age_length_group, reporting_factor = m$reporting_factor,
        reporting_cv = m$reporting_cv, maturity_sampling = "none", light = FALSE),
      Trawl = SimSurvey::sim_fishery(fishery_share = shares$Trawl,
        catch_model = SimSurvey::sim_capture(q = 2, capacity_per_gear = 10000),
        selectivity = selectivity$Trawl, selectivity_by = "length", subset_cells = Trawl,
        gear_amount = 1, duration = 4, gear_unit = "trawl", time_unit = "hour",
        targeting = 1, knowledge_sd = 0.2, exploration = 0.05, max_operations = settings$max_operations,
        sampling_coverage = m$sampling_coverage, sampling_targeting = m$sampling_targeting,
        lengths_cap = m$lengths_cap, ages_cap = m$ages_cap, age_sampling = m$age_sampling,
        age_length_group = m$age_length_group, reporting_factor = m$reporting_factor,
        reporting_cv = m$reporting_cv, maturity_sampling = "none", light = FALSE)), program_settings)
    SimSurvey::sim_fisheries(helpers$.sim_cod_fishing_population(population), n_sims = 1,
      Gillnet = programs$Gillnet, Longline = programs$Longline, Trawl = programs$Trawl)$fishery
  })

  # 3. Independent survey branch: use the same dense annual population, not the
  # large fishery result. Neither survey depletes the other or observes removals.
  surveys <- stage("surveys", {
    rv_sets <- SimSurvey::sim_sets(population, subset_cells = RV,
      n_sims = 1, set_den = settings$survey_density[["RV"]], min_sets = 2, trawl_dim = c(1.5, 0.02))
    sentinel_sets <- SimSurvey::sim_sets(population, subset_cells = Sentinel,
      n_sims = 1, set_den = settings$survey_density[["Sentinel"]], min_sets = 2, trawl_dim = c(1.5, 0.02))
    a <- settings$survey_sampling
    SimSurvey::sim_surveys(population, n_sims = 1,
      RV = SimSurvey::sim_survey(q = selectivity$RV, q_by = "length", custom_sets = rv_sets,
        lengths_cap = a$lengths_cap, ages_cap = a$ages_cap, age_sampling = a$age_sampling,
        age_length_group = a$age_length_group, age_space_group = a$age_space_group,
        maturity_sampling = "none", light = FALSE),
      Sentinel = SimSurvey::sim_survey(q = selectivity$Sentinel, q_by = "length", custom_sets = sentinel_sets,
        lengths_cap = a$lengths_cap, ages_cap = a$ages_cap, age_sampling = a$age_sampling,
        age_length_group = a$age_length_group, age_space_group = a$age_space_group,
        maturity_sampling = "none", light = FALSE))
  })
  provenance$records <- records
  provenance$generated <- Sys.time()
  provenance$session <- utils::capture.output(utils::sessionInfo())
  result <- stage("tables", helpers$.sim_cod_tables(population, fisheries, surveys, settings, selectivity,
    available, provenance, allocations, calibration), cache = FALSE)
  result$sim_cod_truth$provenance$records <- records
  saveRDS(result, file.path(cache_dir, "prepared-data.rds"))
  if (write_data) {
    sim_cod_obs <- result$sim_cod_obs; sim_cod_truth <- result$sim_cod_truth
    save(sim_cod_obs, file = file.path(repo, "data", "sim_cod_obs.rda"), compress = "xz", version = 2)
    save(sim_cod_truth, file = file.path(repo, "data", "sim_cod_truth.rda"), compress = "xz", version = 2)
  }
  message("Prepared data in ", cache_dir, if (write_data) "; packaged data updated." else "; packaged data unchanged.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply the pinned SimSurvey checkout; source this script to choose settings or write_data = TRUE.")
  build_sim_cod_obs(args[[1]])
}
