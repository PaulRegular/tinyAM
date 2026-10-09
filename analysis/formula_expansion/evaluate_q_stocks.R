pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
source(file.path(root, "q_validation.R"))
source(file.path("analysis", "comp_assessments", "R", "run_assessment.R"))
database <- read_committed_database()
out <- file.path(root, "results", "q_stocks")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

candidates <- list(
  afsc_cod_goa_2026 = list(
    logistic = ~ 1 + logistic(age),
    rw = ~ q_block + rw(year),
    ar1 = ~ q_block + ar1(year)),
  afsc_pollock_ebs_2024 = list(
    iid = ~ 0 + q_key + iid(year, by = survey),
    rw = ~ 0 + q_key + rw(year, by = survey),
    ar1 = ~ 0 + q_key + ar1(year, by = survey)),
  afsc_pollock_goa_2024 = list(
    iid = ~ 0 + q_key + environmental_effect + iid(adfg_year, by = adfg_indicator),
    rw = ~ 0 + q_key + environmental_effect + rw(adfg_time, by = adfg_indicator),
    ar1 = ~ 0 + q_key + environmental_effect + ar1(adfg_time, by = adfg_indicator)),
  dfo_cod_4t4vn_2019 = list(
    rw = ~ 0 + q_key + rw(year, by = survey),
    ar1 = ~ 0 + q_key + ar1(year, by = survey)),
  dfo_herring_4tvn_spring_2024 = list(
    rw = ~ 0 + mono(age, by = survey) + cpue_period + rw(year, by = survey),
    ar1 = ~ 0 + mono(age, by = survey) + cpue_period + ar1(year, by = survey),
    logistic = ~ 0 + cpue_period + logistic(age, by = survey)))

diagnostics <- summaries <- effect_parameters <- list()
baseline_settings <- list()
for (id in names(candidates)) {
  cli::cli_inform("Fitting settled baseline: {id}.")
  source <- read_assessment(id, database)
  stock <- new.env(parent = environment())
  sys.source(file.path(.assessment_root, "scripts", "translation", "stocks", paste0(id, ".R")), stock)
  translated <- stock$translate_stock(source)
  baseline_settings[[id]] <- translated$settings
  baseline <- run_assessment(id, database = database)
  row <- baseline$diagnostics
  row$model <- "baseline"
  row$q_formula <- paste(deparse(translated$settings$index_settings$q_form), collapse = " ")
  row$warnings <- ""
  diagnostics[[paste(id, "baseline")]] <- row
  fits <- if (inherits(baseline$fit, "tam_fit")) list(baseline = baseline$fit) else list()
  if (!length(fits)) next
  dir.create(file.path(out, id), showWarnings = FALSE)
  saveRDS(baseline$fit, file.path(out, id, "baseline.rds"))
  if (!is.null(baseline$summary)) {
    summary <- baseline$summary
    summary$model <- "baseline"
    summaries[[paste(id, "baseline")]] <- summary
  }
  obs <- translated$obs
  if (id == "afsc_pollock_goa_2024") {
    adfg <- obs$index$survey == "ADF&G crab/groundfish trawl"
    obs$index$adfg_indicator <- as.numeric(adfg)
    obs$index$adfg_time <- ifelse(adfg, obs$index$year, min(obs$index$year[adfg]))
  }
  for (model in names(candidates[[id]])) {
    settings <- translated$settings
    settings$index_settings$q_form <- candidates[[id]][[model]]
    started <- Sys.time()
    result <- q_capture(function() {
      do.call(fit_tam, c(list(data = obs, years = translated$years, ages = translated$ages,
        silent = TRUE, start_par = tinyAM:::.tam_parameter_summary(baseline$fit, "Estimate")), settings))
    })
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    fitted <- if (inherits(result, "tam_fit")) result else NULL
    row <- .assessment_diagnostics(id, database,
      if (is.null(fitted)) "fit_failed" else if (fitted$is_converged) "converged" else "not_converged",
      fitted, elapsed, if (is.null(result$error)) "" else result$error)
    row$model <- model
    row$q_formula <- paste(deparse(candidates[[id]][[model]]), collapse = " ")
    row$warnings <- paste(result$warnings, collapse = " | ")
    diagnostics[[paste(id, model)]] <- row
    cli::cli_inform("{id}: {model}, {row$status}; gradient {signif(row$max_abs_gradient, 3)}.")
    if (is.null(fitted)) next
    fitted <- .assessment_catch_reporting(fitted, translated$catch_reporting)
    fits[[model]] <- fitted
    saveRDS(fitted, file.path(out, id, paste0(model, ".rds")))
    selected <- grepl("^(sd_q_|phi_q_|q_a50_|q_slope_)", fitted$fixed_par$par)
    effect_parameters[[paste(id, model)]] <- cbind(assessment_id = id, model,
      fitted$fixed_par[selected, , drop = FALSE])
    ref <- q_capture(function() database_to_tam_ref(id,
      if (is.null(translated$comparison_outputs)) database$outputs else translated$comparison_outputs,
      obs = obs, years = translated$years, ages = translated$ages,
      terminal_year = source$assessment$terminal_year[[1]],
      age_plus_group = translated$age_plus_group, comparison_scales = translated$comparison_scales,
      template = fitted, assumptions = source$assumptions,
      comparison_aggregates = translated$comparison_aggregates,
      comparison_age_groups = translated$comparison_age_groups,
      comparison_definitions = translated$comparison_definitions))
    if (inherits(ref, "tam_ref")) {
      summary <- .assessment_comparison_summary(.assessment_percent_differences(ref), id)
      summary$model <- model
      summaries[[paste(id, model)]] <- summary
    }
  }
  saveRDS(list(diagnostics = do.call(rbind, diagnostics), summaries = do.call(rbind, summaries),
    parameters = effect_parameters, database_revision = database$commit,
    baseline_settings = baseline_settings, session = sessionInfo()), file.path(out, "results.rds"))
  labels <- vapply(fits, function(fit) if (isTRUE(fit$is_converged)) "" else " (not converged)", character(1))
  names(fits) <- paste0(names(fits), labels)
  models <- if (inherits(baseline$ref, "tam_ref")) c(list(Accepted = baseline$ref), fits) else fits
  vis_tam(model_list = models, output_file = file.path(out, id, "dashboard.html"),
    open_file = FALSE, background = c(translated$background, "",
      "#### Catchability experiment", "",
      "These candidates change only the catchability formula. Baseline is the settled tinyAM translation. Candidates labelled not converged are diagnostic examples, not supported fitted assessments. Agreement with Accepted trajectories is not a model-selection rule."),
    render_args = list(quiet = TRUE))
}
utils::write.csv(do.call(rbind, diagnostics), file.path(out, "diagnostics.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, summaries), file.path(out, "comparison_summary.csv"), row.names = FALSE)
