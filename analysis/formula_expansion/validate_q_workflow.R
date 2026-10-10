pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion", "results")
source(file.path("analysis", "formula_expansion", "q_validation.R"))
cases <- c("rw_log", "ar1_logit", "intercept_log", "logistic_logit")
rows <- list()
for (case in cases) {
  fit <- readRDS(file.path(root, paste0("q_full_", case, "_replicated_1.rds")))
  before <- serialize(fit$dat, NULL)
  result <- q_capture(function() {
    projection <- update(fit, proj_settings = list(n_proj = 3, n_mean = 1, F_mult = 1),
      silent = TRUE, start_par = tinyAM:::.tam_parameter_summary(fit, "Estimate"))
    stopifnot(sum(projection$dat$is_proj) == 3L, identical(before, serialize(fit$dat, NULL)))
    retro <- fit_retro(projection, folds = 1, start_from_fit = TRUE, progress = FALSE)
    terminal <- max(fit$dat$years[!fit$dat$is_proj])
    stopifnot(all(as.integer(names(retro$fits)) <= terminal))
    for (fold in retro$fits) for (term in fold$dat$q_terms) {
      if (term$type == "logistic") next
      actual <- names(tinyAM:::.tam_parameter_summary(fold, "Estimate")[[term$parameter]])
      expected <- unlist(lapply(term$groups, `[[`, "states"), use.names = FALSE)
      stopifnot(identical(actual, expected), !anyDuplicated(actual))
    }
    hindcast <- fit_hindcast(projection, folds = 1, start_from_fit = TRUE, progress = FALSE)
    stopifnot(all(vapply(hindcast$fits, function(fold) sum(fold$dat$is_proj) == 1L, logical(1))))
    set.seed(510001 + match(case, cases))
    simulated <- sim_tam(projection, n = 2, par_uncertainty = "none", redraw_random = TRUE,
      seed = 510001 + match(case, cases), progress = FALSE)
    stopifnot(length(simulated) > 0L)
    list(projection = check_tam(projection)$is_converged, retro = length(retro$fits),
      hindcast = length(hindcast$fits))
  })
  rows[[case]] <- data.frame(case, projection_converged = isTRUE(result$projection),
    retro_folds = if (is.null(result$retro)) 0L else result$retro,
    hindcast_folds = if (is.null(result$hindcast)) 0L else result$hindcast,
    error = if (is.null(result$error)) "" else result$error,
    warnings = paste(result$warnings, collapse = " | "))
  cli::cli_inform("Workflow: {case} complete.")
}
saveRDS(do.call(rbind, rows), file.path(root, "q_workflow.rds"))
