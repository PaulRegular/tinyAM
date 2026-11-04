# Manual generator check; SimSurvey is deliberately not a tinyAM dependency.
# Rscript data-raw/validate_sim_cod_representation.R /path/to/pinned/SimSurvey

validate_sim_cod_representation <- function(simsurvey_repo, repo = ".") {
  code <- new.env(parent = globalenv())
  sys.source(file.path(repo, "data-raw", "sim_cod_obs.R"), envir = code)
  if (!identical(code$.sim_cod_git(simsurvey_repo), code$sim_cod_settings()$simsurvey_revision)) {
    stop("Use the SimSurvey revision recorded by the data generator.")
  }
  pkgload::load_all(simsurvey_repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(9)
  population <- SimSurvey::sim_abundance(ages = 1:3, years = 1:2,
    R = function(years) rep(100, length(years)),
    Z = function(ages, years) matrix(0.5, length(ages), length(years)),
    growth = SimSurvey::sim_growth(length_group = 10))
  population <- SimSurvey::sim_distribution(population,
    grid = SimSurvey::make_grid(res = c(70, 70)),
    ays_covar = function(x, ages, years, cells) array(0,
      c(length(ages), length(years), length(cells)),
      dimnames = list(age = ages, year = years, cell = cells)),
    depth_par = function(d) ifelse(d$x < 0, 0, -1000))
  population <- SimSurvey::round_sim(population)
  configurations <- list(
    FleetA = SimSurvey::sim_fishery(fishery_share = 0.2,
      catch_model = SimSurvey::sim_capture(q = 5000, capacity_per_gear = 10),
      selectivity = function(age) rep(1, length(age)), sampling_coverage = 0.5,
      lengths_cap = 10, ages_cap = 3, max_operations = 20,
      maturity_sampling = "none", light = FALSE),
    FleetB = SimSurvey::sim_fishery(fishery_share = 0.1,
      catch_model = SimSurvey::sim_capture(q = 3000, capacity_per_gear = 8),
      selectivity = function(age) stats::plogis(age - 2), sampling_coverage = 0.6,
      lengths_cap = 8, ages_cap = 2, max_operations = 20,
      maturity_sampling = "none", light = FALSE))
  results <- list()
  for (empty_age_year in c(FALSE, TRUE)) {
    dense <- data.table::copy(population)
    if (empty_age_year) {
      dense$sp_N$N[dense$sp_N$age == 3 & dense$sp_N$year == 2] <- 0
      dense <- SimSurvey::round_sim(dense)
    }
    original <- data.table::copy(dense)
    sparse <- code$.sim_cod_fishing_population(dense)
    stopifnot(identical(dense, original), identical(sparse$grid, dense$grid),
      identical(SimSurvey::round_sim(sparse)$N, dense$N))
    set.seed(802)
    dense_result <- do.call(SimSurvey::sim_fisheries,
      c(list(sim = dense, n_sims = 1), configurations))
    dense_rng <- .Random.seed
    set.seed(802)
    sparse_result <- code$.sim_cod_run_fisheries(dense, configurations)
    sparse_rng <- .Random.seed
    stopifnot(identical(dense_rng, sparse_rng), identical(dense, original),
      identical(sparse_result$sp_N, dense$sp_N),
      identical(sparse_result$N, dense_result$N),
      identical(sparse_result$N_at_length, dense_result$N_at_length),
      identical(sparse_result$B, dense_result$B), identical(sparse_result$S, dense_result$S))
    dense_fishery <- dense_result$fishery
    sparse_fishery <- sparse_result$fishery
    dense_ledger <- dense_fishery$mortality_ledger
    sparse_ledger <- sparse_fishery$mortality_ledger
    dense_ledger <- as.data.frame(dense_ledger[dense_ledger$starting > 0, ])
    sparse_ledger <- as.data.frame(sparse_ledger[sparse_ledger$starting > 0, ])
    dense_fishery$mortality_ledger <- sparse_fishery$mortality_ledger <- NULL
    stopifnot(isTRUE(all.equal(dense_fishery, sparse_fishery, tolerance = 0)),
      isTRUE(all.equal(dense_ledger, sparse_ledger, tolerance = 0, check.attributes = FALSE)))
    results[[length(results) + 1L]] <- data.frame(empty_age_year = empty_age_year,
      dense_rows = nrow(dense$sp_N), sparse_rows = nrow(sparse$sp_N),
      outputs_exact = TRUE, rng_exact = TRUE, dense_spatial_restored = TRUE)
  }
  result <- do.call(rbind, results)
  print(result, row.names = FALSE)
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Supply the pinned SimSurvey source checkout.")
  validate_sim_cod_representation(args[[1]])
}
