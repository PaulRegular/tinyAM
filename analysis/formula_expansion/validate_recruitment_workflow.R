pkgload::load_all(quiet = TRUE)
future::plan(future::sequential)
root <- file.path("analysis", "formula_expansion", "results")
cases <- c("bh_iid", "ricker_ar1", "cov_rw")
if (!all(file.exists(file.path(root, paste0("recruitment_", cases, ".rds"))))) {
  cli::cli_abort("Complete the recruitment recovery study before running workflow validation.")
}
rows <- list()
for (case in cases) {
  fit <- readRDS(file.path(root, paste0("recruitment_", case, ".rds")))
  before <- serialize(fit$dat, NULL)
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    projection <- update(fit, proj_settings = list(n_proj = 3, n_mean = 1, F_mult = 1),
      silent = TRUE, start_par = tinyAM:::.tam_parameter_summary(fit, "Estimate"))
    stopifnot(sum(projection$dat$is_proj) == 3L, identical(before, serialize(fit$dat, NULL)))
    retro <- fit_retro(projection, folds = 2, start_from_fit = TRUE, progress = FALSE)
    for (fold in retro$fits) {
      p <- tinyAM:::.tam_parameter_summary(fold, "Estimate")
      stopifnot(identical(names(p$log_r), as.character(fold$dat$years[fold$dat$rec$eligible])))
      if (length(fold$dat$rec$boundary) > 1L) stopifnot(identical(names(p$log_r_init),
        as.character(fold$dat$years[fold$dat$rec$boundary[-1L]])))
    }
    hindcast <- fit_hindcast(projection, folds = 2, start_from_fit = TRUE, progress = FALSE)
    stopifnot(all(vapply(hindcast$fits, function(fold) sum(fold$dat$is_proj) == 1L, logical(1))))
    simulated <- sim_tam(projection, n = 2, par_uncertainty = "none", redraw_random = TRUE,
      seed = 810001L + match(case, cases), progress = FALSE)
    stopifnot(length(simulated) > 0L, identical(before, serialize(fit$dat, NULL)))
    tab <- tidy_recruitment(projection)
    stopifnot(all(is.finite(tab$residual$est)), all(is.finite(tab$prediction$est)))
    if (!is.null(fit$dat$rec$curve)) {
      stopifnot(all(tab$pairs$parent_year == tab$pairs$year - fit$dat$rec$curve$lag))
      if (isTRUE(projection$sdrep$pdHess)) stopifnot(all(is.finite(tab$curve$se))) else
        stopifnot(all(is.na(tab$curve$se)))
    }
    dashboard <- normalizePath(file.path(root, paste0("recruitment_", case, "_dashboard.html")), mustWork = FALSE)
    vis_tam(model_list = list(historical = fit, projected = projection),
      output_file = dashboard,
      open_file = FALSE, render_args = list(quiet = TRUE))
    connection <- file(dashboard, "rb")
    html <- paste(readLines(connection, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    close(connection)
    stopifnot(grepl('id="recruitment-residuals"', html, fixed = TRUE),
      grepl('id="stockrecruitment"', html, fixed = TRUE) == !is.null(fit$dat$rec$curve))
    list(projection = projection$is_converged,
      retro = vapply(retro$fits, function(fold) fold$is_converged, logical(1)),
      hindcast = vapply(hindcast$fits, function(fold) fold$is_converged, logical(1)))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  rows[[case]] <- data.frame(case, projection_converged = isTRUE(result$projection),
    retro_requested = 3L, retro_converged = sum(result$retro),
    hindcast_requested = 3L, hindcast_converged = sum(result$hindcast),
    error = if (is.null(result$error)) "" else result$error,
    warnings = paste(unique(warnings), collapse = " | "))
  cli::cli_inform("Recruitment workflow: {case} complete.")
  saveRDS(do.call(rbind, rows), file.path(root, "recruitment_workflow.rds"))
}
