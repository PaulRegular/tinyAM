# Extend the saved pilot on the same training and holdout histories.
.oe_pilot_extra_models <- c("gillnet_smooth3", "rv_smooth3", "gillnet_iid",
  "rv_iid", "rv_design_annual", "rv_cadigan")
.oe_pilot_sd_years <- c(5L, 10L, 16L)
.oe_pilot_sd_draws <- 20L

.oe_pilot_extra_fit <- function(base, train, holdout, truth, model, baseline, calibration) {
  spec <- function(input, model, training_rv) {
    result <- .oe_pilot_spec(input, "constant", training_rv)
    switch(model,
      gillnet_smooth3 = result$settings$catch_settings$sd_form <- ~ s(age, k = 3),
      rv_smooth3 = result$settings$index_settings$sd_form <- ~ s(age, k = 3),
      gillnet_iid = result$settings$catch_settings$sd_form <- ~ iid(age),
      rv_iid = result$settings$index_settings$sd_form <- ~ iid(age),
      rv_design_annual = {
        offsets <- oe_pilot_rv_offset(training_rv, model)
        result <- oe_pilot_apply_rv_offset(result$input, result$settings, offsets)
      },
      rv_cadigan = {
        if (!isTRUE(baseline$usable) || !NROW(baseline$holdout)) {
          stop("The matched constant-SD fit has no usable training predictions.")
        }
        # These columns are calculated from the training fit. The holdout
        # observation and its residual are not used to construct the SD.
        index <- subset(baseline$holdout, component == "RV", c("year", "age", "pred", "sd"))
        offsets <- oe_pilot_rv_offset(training_rv, model, index, calibration)
        result <- oe_pilot_apply_rv_offset(result$input, result$settings, offsets)
      },
      stop("Unknown SD candidate ", model))
    result
  }
  .oe_pilot_fit(base, train, holdout, truth, model, specification = spec)
}

.oe_pilot_extra_records <- function(stock, cache, old_cache, calibration) {
  data <- readRDS(file.path(old_cache, "inputs.rds"))
  seeds <- .oe_pilot_seed_pairs(stock)
  fit <- function(base, train, holdout, truth, model, draw) {
    baseline <- if (model == "rv_cadigan") {
      readRDS(file.path(old_cache, sprintf("fit-%s-constant.rds", draw)))$value
    }
    .oe_pilot_extra_fit(base, train, holdout, truth, model, baseline, calibration)
  }
  .oe_pilot_fit_records(stock, seeds, cache, data$value, data$error,
    models = .oe_pilot_extra_models, fit = fit)
}

.oe_pilot_reference_draws <- function(stock, i, cache, old_cache, archived) {
  draws <- vector("list", length(.oe_pilot_sd_years) * .oe_pilot_sd_draws)
  pop <- NULL
  for (year in .oe_pilot_sd_years) {
    files <- sprintf("draw-%02d-cal-%02d.rds", year, seq_len(.oe_pilot_sd_draws))
    paths <- file.path(archived, files)
    if (!all(file.exists(paths))) {
      if (is.null(pop)) {
        admit()
        pop <- oe_pilot_population(stock$population_seed, stock$field_seed, grid)
      }
      snapshot <- oe_snapshot(pop, year)
      removals <- required(readRDS(file.path(old_cache, sprintf("removals-%02d.rds", year))))
    }
    for (replicate in seq_along(files)) {
      # Match the independent draws used by the original mean diagnostic.
      seed <- 910000000L + 10000L * i + 50L * (year - 1L) + replicate
      reused <- file.exists(paths[replicate])
      attempt <- if (reused) readRDS(paths[replicate]) else {
        .oe_pilot_draw(snapshot, removals, c(Gillnet = seed, RV = seed + 2000L),
          cache, sprintf("cal-%02d", replicate))
      }
      rows <- if (is.null(attempt$error)) {
        do.call(rbind, lapply(attempt$value, function(d) d[c("component", "year", "age", "obs", "truth")]))
      }
      j <- (match(year, .oe_pilot_sd_years) - 1L) * .oe_pilot_sd_draws + replicate
      draws[[j]] <- list(rows = if (!is.null(rows)) transform(rows,
          population_id = stock$population_id, regime = stock$regime, replicate = replicate),
        status = data.frame(population_id = stock$population_id, year = year,
          replicate = replicate, gillnet_seed = seed, rv_seed = seed + 2000L,
          reused = reused, usable = !is.null(rows),
          warnings = paste(attempt$warnings, collapse = " | "),
          error = if (is.null(attempt$error)) "" else attempt$error))
    }
  }
  list(rows = .oe_pilot_stack(draws, "rows"), attempts = .oe_pilot_stack(draws, "status"))
}

oe_pilot_extend <- function() {
  stopifnot(identical(grid, 40L))
  limit <- as.integer(Sys.getenv("TINYAM_OE_POPULATION_LIMIT"))
  stopifnot(limit %in% c(3L, 24L))
  files <- file.path(study, c("observation_error_pilot_populations.R",
    "observation_error_fit.R", "observation_error_design_precision.R",
    "observation_error_pilot_rv_offsets.R"))
  for (file in c(files, diag_file)) source(file, encoding = "UTF-8")
  original <- file.path(out, "population-manifest.rds")
  previous <- readRDS(original)
  manifest <- previous$seed_design$populations
  immutable <- setdiff(names(previous$sources), pilot_file)
  stopifnot(identical(tools::md5sum(immutable), previous$sources[immutable]),
    identical(basename(simsurvey), paste0("SimSurvey-", previous$simsurvey_revision)),
    identical(unname(tools::md5sum(sort(list.files(file.path(simsurvey, "R"), "[.]R$", full.names = TRUE)))), unname(previous$simsurvey_sources)),
    identical(vapply(names(previous$versions), function(p) {
      as.character(utils::packageVersion(p))
    }, character(1)), previous$versions),
    identical(as.character(getRversion()), previous$R))
  calibration <- oe_pilot_cadigan_calibration(file.path(study, "results",
    "observation_error_seed21", "study-grid40", "summary.rds"))
  extension <- file.path(out, "extension")
  dir.create(extension, showWarnings = FALSE)
  archive <- file.path(out, "mean-diagnostic")
  verify_manifest(list(original = tools::md5sum(original),
    sources = tools::md5sum(c(files, diag_file, pilot_file, extension_file)),
    calibration = calibration, years = .oe_pilot_sd_years, draws = .oe_pilot_sd_draws),
    file.path(extension, "manifest.rds"))

  for (i in seq_len(limit)) {
    stock <- manifest[i, , drop = FALSE]
    cache <- file.path(extension, stock$population_id)
    old_cache <- file.path(out, stock$population_id)
    dir.create(cache, showWarnings = FALSE)
    cat("Additional fits for", stock$population_id, "\n")
    flush.console()
    fit_path <- file.path(cache, "fits.rds")
    if (!file.exists(fit_path)) {
      save_checkpoint(.oe_pilot_extra_records(stock, cache, old_cache, calibration), fit_path)
    }
  }
  for (i in seq_len(limit)) {
    stock <- manifest[i, , drop = FALSE]
    cache <- file.path(extension, stock$population_id)
    old_cache <- file.path(out, stock$population_id)
    cat("Independent SD reference for", stock$population_id, "\n")
    flush.console()
    ref_path <- file.path(cache, "sd-reference.rds")
    if (!file.exists(ref_path)) {
      save_checkpoint(.oe_pilot_reference_draws(stock, i, cache, old_cache,
        file.path(archive, stock$population_id)), ref_path)
    }
    invisible(gc())
  }
}

oe_pilot_extension_report <- function(out) {
  extension <- file.path(out, "extension")
  manifest <- readRDS(file.path(out, "population-manifest.rds"))$seed_design
  design <- readRDS(file.path(extension, "manifest.rds"))
  load_results <- function(root, file) {
    paths <- file.path(root, manifest$populations$population_id, file)
    lapply(paths[file.exists(paths)], readRDS)
  }
  original <- load_results(out, "result.rds")
  added <- load_results(extension, "fits.rds")
  reference <- load_results(extension, "sd-reference.rds")
  combined <- c(original, added)
  source(file.path(study, "observation_error_pilot_evaluation.R"), encoding = "UTF-8")
  summary <- list(manifest = manifest$populations, truth = transform(manifest$ssb, SSB = ssb),
    models = c(.oe_pilot_models, .oe_pilot_extra_models),
    fits = .oe_pilot_stack(combined, "fits"), annual = .oe_pilot_stack(combined, "annual"),
    holdout = .oe_pilot_stack(combined, "holdout"), fisheries = .oe_pilot_stack(original, "fisheries"),
    complete = all(lengths(list(original, added, reference)) == nrow(manifest$populations)),
    cadigan_calibration = design$calibration, sd_reference_draws = design$draws,
    sd_reference_years = design$years, sd_reference_attempts = .oe_pilot_stack(reference, "attempts"))
  means <- file.path(out, "mean-diagnostic.rds")
  if (file.exists(means)) summary$mean_diagnostic <- readRDS(means)
  summary$evaluation <- oe_pilot_evaluate(summary$fits, summary$annual, summary$holdout,
    sd_validation = .oe_pilot_stack(reference, "rows"))
  path <- file.path(extension, "summary.rds")
  save_checkpoint(summary, path)
  rmarkdown::render(file.path(study, "observation_error_pilot_report.Rmd"),
    params = list(result = normalizePath(path, winslash = "/")), output_dir = dirname(normalizePath(path)),
    quiet = TRUE)
}
