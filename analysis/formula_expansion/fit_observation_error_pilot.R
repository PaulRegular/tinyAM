# Sourced by the population-pilot phases of simulate_observation_error.R.
.oe_pilot_models <- c("constant", "gillnet_groups", "gillnet_age", "rv_pooled", "rv_age")

.oe_pilot_spec <- function(input, model, training_rv) {
  settings <- .oe_fit_settings("", "constant")
  settings$index_settings$q_form <- ~ 0 + logistic(age)
  settings$F_settings$mu_form <- ~ factor(year) + double_logistic(age)
  fallback <- 0L
  if (model == "gillnet_groups") {
    input$catch$sd_age_group <- cut(input$catch$age, c(0, 4, 10, Inf),
      labels = c("1–4", "5–10", "11–14+"))
    settings$catch_settings$sd_form <- ~ sd_age_group
  } else if (model == "gillnet_age") {
    settings$catch_settings$sd_form <- ~ s(age, k = 5)
  } else if (model == "rv_age") {
    settings$index_settings$sd_form <- ~ s(age, k = 5)
  } else if (model == "rv_pooled") {
    rv <- training_rv
    sd <- numeric(nrow(rv))
    for (year in sort(unique(rv$year))) {
      test <- rv$year == year
      pooled <- oe_pooled_rv_sd(rv[!test, ], rv$age[test])
      sd[test] <- pooled
      fallback <- fallback + sum(attr(pooled, "fallback") & is.finite(rv$obs[test]) & rv$obs[test] > 0)
    }
    input$index$rv_log_sd <- sd[match(.oe_fit_key(input$index), .oe_fit_key(rv))]
    input$index$rv_sd_scale <- 1
    settings$index_settings$sd_supplied <- ~ rv_log_sd
    settings$index_settings$sd_form <- ~ 0 + rv_sd_scale
  }
  list(input = input, settings = settings, fallback = fallback)
}

.oe_pilot_draw <- function(snapshot, removals, seeds, cache, label) {
  checkpoint(file.path(cache, sprintf("draw-%02d-%s.rds", snapshot$years, label)), function() {
    rv <- oe_rv_draw(snapshot, seeds[["RV"]], efforts = "standard")$rows
    gillnet <- oe_gillnet_draw(removals, seeds[["Gillnet"]], efforts = "standard")
    rv <- subset(as.data.frame(rv), bin == "modeled")
    gillnet <- subset(as.data.frame(gillnet), bin == "modeled")
    stopifnot(all(rv$n_primary == 160L), setequal(rv$age, 1:14), setequal(gillnet$age, 1:14),
      nrow(rv) == 14L, nrow(gillnet) == 14L, all(rv$year == snapshot$years),
      all(gillnet$year == snapshot$years))
    list(catch = gillnet, index = rv)
  })
}

.oe_pilot_seed_pairs <- function(stock) {
  role <- function(draw, purpose) {
    field <- paste0(if (draw == "repeat") "repeat_" else "", purpose, "_seed")
    c(Gillnet = stock[[paste0("gillnet_", field)]], RV = stock[[paste0("rv_", field)]])
  }
  draws <- if (stock$repeat_monitoring) c("primary", "repeat") else "primary"
  stats::setNames(lapply(draws, function(draw) {
    list(fit = role(draw, "fit"), holdout = role(draw, "holdout"))
  }), draws)
}

.oe_pilot_stack <- function(values, name) {
  do.call(rbind, lapply(values, function(x) x[[name]]))
}

.oe_pilot_removals <- function(snapshot, seed) {
  # Capture probabilities depend only on age within each cell. Binomial
  # capture and sampling without replacement have the same age margins
  # after pooling lengths. The F_full target is in numbers of fish.
  ages <- snapshot$ages
  biology <- SimSurvey::stock_biology(snapshot)$at_age
  weights <- biology$mean_weight_kg[match(ages, biology$age)]
  key <- diag(length(ages))
  dimnames(key) <- list(length = ages, age = ages)
  fishing <- snapshot
  fishing$sim_length <- function(age, length_age_key) key[, as.character(age), drop = FALSE]
  fishing$sim_weight <- function(length) weights[match(length, ages)]
  result <- oe_gillnet_truth(fishing, seed)
  p <- result$population
  stopifnot(identical(as.integer(p$age), as.integer(ages)),
    isTRUE(all.equal(p$weight, weights, tolerance = 1e-12)),
    isTRUE(all.equal(p$SSB, as.numeric(snapshot$S) / 1000, tolerance = 1e-12)),
    all(p$N == as.numeric(snapshot$N) / 1000), all(p$Z == as.numeric(snapshot$Z)))
  # The study uses counts and operation-level age samples. The grouped
  # simulator's biomass would use mean weights rather than individual sizes.
  result$fishery_annual <- as.data.frame(result$fishery_annual)[c("year", "target", "achieved",
    "shortfall", "n_operations", "status")]
  result$fishery_reference <- NULL
  result
}

.oe_pilot_fit <- function(base, train, holdout, truth, model,
    specification = .oe_pilot_spec) {
  input <- .oe_fit_inputs(base$data, train, train, NULL, "constant")$data
  spec <- specification(input, model, train$index)
  fit <- do.call(tinyAM::fit_tam, c(list(data = spec$input, silent = TRUE), spec$settings))
  checks <- tinyAM::check_tam(fit)
  ans <- list(usable = isTRUE(checks$is_converged), numerical = checks$numerical,
    rv_fallback = spec$fallback)
  if (!ans$usable) return(ans)
  ans$annual <- .oe_fit_population(fit, truth)$annual
  if (!is.null(holdout)) {
    ans$holdout <- tryCatch(.oe_fit_predictions(fit, train, holdout, truth),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
        NULL
      })
  }
  ans
}

.oe_pilot_fit_records <- function(stock, seeds, cache, data, data_error,
    models = .oe_pilot_models,
    fit = function(base, train, holdout, truth, model, draw) {
      .oe_pilot_fit(base, train, holdout, truth, model)
    }) {
  keys <- function(x, draw, model) {
    if (is.null(x)) return(NULL)
    transform(x, population_id = stock$population_id, regime = stock$regime,
      draw = draw, model = model)
  }
  records <- lapply(names(seeds), function(draw) {
    pair <- if (is.null(data)) NULL else data$histories[[draw]]
    lapply(models, function(model) {
      result <- checkpoint(file.path(cache, sprintf("fit-%s-%s.rds", draw, model)), function() {
        if (!is.null(data_error)) stop(data_error, call. = FALSE)
        if (is.null(pair$fit)) stop("Training history is incomplete. See the draw checkpoints.")
        fit(data$base, pair$fit, pair$holdout, data$truth, model, draw)
      })
      x <- result$value
      list(fits = keys(data.frame(usable = isTRUE(x$usable),
          warnings = paste(result$warnings, collapse = " | "),
          error = if (is.null(result$error)) "" else result$error,
          rv_fallback = if (is.null(x$rv_fallback)) NA_integer_ else x$rv_fallback,
          elapsed = result$elapsed), draw, model),
        annual = keys(x$annual, draw, model),
        holdout = keys(x$holdout, draw, model))
    })
  })
  records <- unlist(records, recursive = FALSE)
  stats::setNames(lapply(c("fits", "annual", "holdout"), function(nm) {
    .oe_pilot_stack(records, nm)
  }), c("fits", "annual", "holdout"))
}

.oe_pilot_inputs <- function(stock, seeds, cache) {
  pop <- oe_pilot_population(stock$population_seed, stock$field_seed, grid)
  snapshots <- vector("list", length(pop$years))
  paths <- lapply(seeds, function(purposes) lapply(purposes, function(x) vector("list", length(pop$years))))
  for (i in seq_along(pop$years)) {
    admit()
    year <- pop$years[i]
    snapshot <- oe_snapshot(pop, year)
    removals <- required(checkpoint(file.path(cache, sprintf("removals-%02d.rds", year)),
      function() .oe_pilot_removals(snapshot, stock$removals_seed + year - 1L)))
    snapshots[[i]] <- removals[c("population", "annual", "fishery_annual")]
    for (draw in names(seeds)) {
      for (purpose in names(seeds[[draw]])) {
        label <- paste(draw, purpose, sep = "-")
        paths[[draw]][[purpose]][[i]] <- .oe_pilot_draw(snapshot, removals,
          seeds[[draw]][[purpose]] + year - 1L, cache, label)
      }
    }
    rm(snapshot, removals)
    invisible(gc())
  }
  truth <- oe_assessment_inputs(snapshots)
  history <- function(attempts) {
    if (!all(vapply(attempts, function(x) is.null(x$error), logical(1)))) return(NULL)
    values <- lapply(attempts, "[[", "value")
    lapply(stats::setNames(c("catch", "index"), c("catch", "index")), function(nm) {
      .oe_pilot_stack(values, nm)
    })
  }
  list(truth = truth, base = .oe_fit_data(truth),
    fisheries = .oe_pilot_stack(snapshots, "fishery_annual"),
    histories = lapply(paths, function(x) lapply(x, history)))
}

oe_pilot_run <- function() {
  stopifnot(identical(grid, 40L))
  limit <- as.integer(Sys.getenv("TINYAM_OE_POPULATION_LIMIT"))
  stopifnot(limit %in% c(3L, 24L))
  files <- file.path(study, c("observation_error_pilot_populations.R",
    "observation_error_fit.R", "observation_error_design_precision.R"))
  for (file in c(files, diag_file)) source(file, encoding = "UTF-8")
  design <- oe_pilot_population_manifest()
  manifest <- design$populations
  provenance <- list(grid = grid, seed_design = design,
    sources = tools::md5sum(c(files, pilot_file, sim_file, diag_file,
      "DESCRIPTION", "NAMESPACE", sort(list.files("R", "[.]R$", full.names = TRUE)))),
    simsurvey_revision = .observation_error_revision,
    simsurvey_sources = tools::md5sum(sort(list.files(file.path(simsurvey, "R"), "[.]R$", full.names = TRUE))),
    versions = vapply(c("tinyAM", "SimSurvey", "RTMB", "data.table", "mgcv"),
      function(p) as.character(utils::packageVersion(p)), character(1)),
    R = as.character(getRversion()))
  verify_manifest(provenance, file.path(out, "population-manifest.rds"))
  for (i in seq_len(limit)) {
    stock <- manifest[i, , drop = FALSE]
    cache <- file.path(out, stock$population_id)
    dir.create(cache, showWarnings = FALSE)
    result_path <- file.path(cache, "result.rds")
    if (file.exists(result_path)) next
    cat("Population", stock$population_id, "of", nrow(manifest), stock$regime, "\n")
    flush.console()
    seeds <- .oe_pilot_seed_pairs(stock)
    inputs <- checkpoint(file.path(cache, "inputs.rds"),
      function() .oe_pilot_inputs(stock, seeds, cache))
    if (!is.null(inputs$value)) {
      expected <- subset(design$ssb, population_id == stock$population_id)
      truth <- inputs$value$truth$annual
      stopifnot(setequal(truth$year, expected$year),
        isTRUE(all.equal(truth$SSB, expected$ssb[match(truth$year, expected$year)], tolerance = 1e-5)))
    }
    result <- .oe_pilot_fit_records(stock, seeds, cache, inputs$value, inputs$error)
    if (!is.null(inputs$value)) {
      result$fisheries <- transform(as.data.frame(inputs$value$fisheries),
        population_id = stock$population_id, regime = stock$regime)
    }
    result$population_id <- stock$population_id
    save_checkpoint(result, result_path)
    cat("Completed", stock$population_id, sum(result$fits$usable), "usable fits of",
      nrow(result$fits), "\n")
    rm(inputs, result)
    invisible(gc())
  }
}

oe_pilot_report <- function(out) {
  manifest <- readRDS(file.path(out, "population-manifest.rds"))$seed_design
  results <- lapply(file.path(out, manifest$populations$population_id, "result.rds"), function(p) {
    if (file.exists(p)) readRDS(p)
  })
  results <- Filter(Negate(is.null), results)
  stopifnot(length(results) > 0L)
  source(file.path(study, "observation_error_pilot_evaluation.R"), encoding = "UTF-8")
  summary <- list(manifest = manifest$populations, truth = transform(manifest$ssb, SSB = ssb),
    fits = .oe_pilot_stack(results, "fits"), annual = .oe_pilot_stack(results, "annual"),
    holdout = .oe_pilot_stack(results, "holdout"), fisheries = .oe_pilot_stack(results, "fisheries"),
    complete = length(results) == nrow(manifest$populations))
  summary$evaluation <- oe_pilot_evaluate(summary$fits, summary$annual, summary$holdout)
  path <- file.path(out, "summary.rds")
  save_checkpoint(summary, path)
  rmarkdown::render(file.path(study, "observation_error_pilot_report.Rmd"),
    params = list(result = normalizePath(path, winslash = "/")), output_dir = dirname(normalizePath(path)),
    quiet = TRUE)
}
