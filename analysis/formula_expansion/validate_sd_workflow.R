pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion", "results")
cases <- c("iid_catch", "rw_index")
rows <- list()
for (case in cases) {
  fit <- readRDS(file.path(root, paste0("sd_full_", case, "_1.rds")))
  before <- serialize(fit$dat, NULL)
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    projection <- update(fit, proj_settings = list(n_proj = 3, n_mean = 1, F_mult = 1),
      silent = TRUE, start_par = tinyAM:::.tam_parameter_summary(fit, "Estimate"))
    stopifnot(sum(projection$dat$is_proj) == 3L, identical(before, serialize(fit$dat, NULL)))
    retro <- fit_retro(projection, folds = 1, start_from_fit = TRUE, progress = FALSE)
    terminal <- max(fit$dat$years[!fit$dat$is_proj])
    stopifnot(all(as.integer(names(retro$fits)) <= terminal))
    for (fold in retro$fits) for (term in tinyAM:::.formula_terms(fold$dat)) {
      actual <- names(tinyAM:::.tam_parameter_summary(fold, "Estimate")[[term$parameter]])
      expected <- unlist(lapply(term$groups, `[[`, "states"), use.names = FALSE)
      stopifnot(identical(actual, expected), !anyDuplicated(actual))
    }
    hindcast <- fit_hindcast(projection, folds = 1, start_from_fit = TRUE, progress = FALSE)
    stopifnot(all(vapply(hindcast$fits, function(fold) sum(fold$dat$is_proj) == 1L, logical(1))))
    simulated <- sim_tam(projection, n = 2, par_uncertainty = "none", redraw_random = TRUE,
      seed = 640001 + match(case, cases), progress = FALSE)
    stopifnot(length(simulated) > 0L, identical(before, serialize(fit$dat, NULL)))
    effects <- tidy_tam(model_list = list(original = fit, projected = projection))$formula_effects$levels
    stopifnot(length(effects) > 0L, all(vapply(effects, function(d)
      all(is.finite(d$est)) && any(startsWith(d$component, "sd_")), logical(1))))
    for (effect in effects) {
      plot <- plotly::plotly_build(tinyAM:::.plot_formula_effect(effect))
      stopifnot(plot$x$layout$yaxis$title == "Effect on log SD")
    }
    vis_tam(model_list = list(original = fit, projected = projection),
      output_file = normalizePath(file.path(root, paste0("sd_", case, "_dashboard.html")), mustWork = FALSE),
      open_file = FALSE, render_args = list(quiet = TRUE))
    list(projection = projection$is_converged, retro = length(retro$fits), hindcast = length(hindcast$fits))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  rows[[case]] <- data.frame(case, projection_converged = isTRUE(result$projection),
    retro_folds = if (is.null(result$retro)) 0L else result$retro,
    hindcast_folds = if (is.null(result$hindcast)) 0L else result$hindcast,
    error = if (is.null(result$error)) "" else result$error,
    warnings = paste(unique(warnings), collapse = " | "))
  cli::cli_inform("SD workflow: {case} complete.")
}
saveRDS(do.call(rbind, rows), file.path(root, "sd_workflow.rds"))
