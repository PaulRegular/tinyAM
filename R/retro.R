.retro_projection <- function(settings) {
  defaults <- list(n_proj = 1, n_mean = 1, catch = "observed", F_mult = NULL)
  if (is.null(settings)) return(defaults)
  out <- .merge_settings(defaults, settings, "proj_settings")
  if (out$n_proj < 1) cli::cli_abort("Hindcasts require n_proj of at least one.")
  key <- if (is.null(out$catch)) "F_mult" else "catch"
  if (!identical(out[[key]], "observed") || key != "catch") {
    out[[key]] <- unname(.projection_vector(unname(out[[key]]), seq_len(out$n_proj), key))
  }
  out
}

# Resolve removals from raw biological ages before plus-group pooling.
.hindcast_controls <- function(fit, years, settings) {
  kind <- if (is.null(settings$catch)) "F_mult" else "catch"
  observed <- identical(settings$catch, "observed")
  out <- data.frame(year = years, horizon = seq_along(years), control = kind,
    target = if (observed) rep(NA_real_, length(years)) else settings[[kind]][seq_along(years)],
    target_source = if (observed) "observed catch-at-age x supplied weight" else if (kind == "catch") "supplied catch" else "F multiplier",
    target_unit = if (kind == "catch") .population_unit(fit$dat, "total_yield") else "multiplier", error = NA_character_)
  source <- fit$refit_args$data
  ages <- sort(unique(c(fit$dat$ages, fit$dat$plus_ages)))
  if (observed) for (i in seq_along(years)) {
    d <- source$catch[source$catch$year == years[i], ]
    w <- source$weight[source$weight$year == years[i], ]
    counts <- d$obs[match(ages, d$age)]; weights <- w$obs[match(ages, w$age)]
    if (any(!is.finite(counts) | counts < 0) || any(!is.finite(weights) | weights < 0)) {
      out$error[i] <- "Observed catch biomass unavailable: incomplete catch-at-age or weights; supply catch targets or F_mult."
    } else out$target[i] <- sum(counts * weights)
  }
  out
}

.retro_fold <- function(fit, terminal, final_year, settings, start_from_fit, grad_tol) {
  reference <- terminal == final_year
  n_proj <- if (reference || is.null(settings)) 0 else min(settings$n_proj, final_year - terminal)
  diagnostic <- data.frame(fold = terminal, is_reference = reference,
    requested_n_proj = if (is.null(settings) || reference) 0 else settings$n_proj,
    available_n_proj = n_proj, n_mean = if (is.null(settings)) NA_real_ else settings$n_mean,
    is_converged = FALSE, n_scored_index = 0L, error = NA_character_, warnings = "")
  controls <- data.frame(); warnings <- character()
  result <- withCallingHandlers(tryCatch({
    projection <- NULL
    if (n_proj > 0) {
      controls <- .hindcast_controls(fit, terminal + seq_len(n_proj), settings)
      if (any(!is.na(controls$error))) cli::cli_abort(paste(unique(stats::na.omit(controls$error)), collapse = "; "))
      projection <- list(n_proj = n_proj, n_mean = settings$n_mean)
      projection[[controls$control[1]]] <- controls$target
    }
    r <- if (reference && !any(fit$dat$is_proj)) fit else update(fit,
      years = min(fit$dat$years[!fit$dat$is_proj]):terminal, ages = fit$dat$ages,
      proj_settings = projection, start_par = if (start_from_fit) fit else NULL, silent = TRUE)
    r$grad_tol <- grad_tol
    r$diagnostics <- check_tam(r, grad_tol = grad_tol)
    r$is_converged <- isTRUE(r$diagnostics$is_converged)
    diagnostic$is_converged <- r$is_converged
    if (!r$is_converged) diagnostic$error <- "Numerical convergence checks did not pass."
    if (n_proj && r$is_converged) {
      d <- r$obs_pred$index; d$fold <- terminal; d$horizon <- d$year - terminal
      p <- .forecast_pairs(d, fit$obs_pred$index)
      diagnostic$n_scored_index <- sum(is.finite(p$pred) & p$pred > 0 & is.finite(p$obs) & p$obs > 0)
    }
    r
  }, error = function(e) { diagnostic$error <<- conditionMessage(e); NULL }),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  diagnostic$warnings <- paste(unique(warnings), collapse = "\n")
  list(fit = result, diagnostic = diagnostic, controls = controls)
}

.retro_labels <- function(models, expressions) {
  labels <- names(models)
  if (is.null(labels)) labels <- rep("", length(models))
  empty <- is.na(labels) | !nzchar(labels)
  labels[empty] <- vapply(expressions[empty], function(x) paste(deparse(x), collapse = ""), character(1))
  if (anyDuplicated(labels)) cli::cli_abort("Supply distinct names for the fitted models, e.g. fit_retro(fit, alternative = fit2).")
  labels
}

.run_retros <- function(models, folds, start_from_fit, hindcast, proj_settings, grad_tol, progress, globals) {
  for (model in models) .require_tam_fit(model, arg = "fit")
  for (flag in list(start_from_fit, hindcast, progress)) {
    if (!is.logical(flag) || length(flag) != 1 || is.na(flag)) cli::cli_abort("Retrospective flags must be single TRUE/FALSE values.")
  }
  if (!hindcast && !is.null(proj_settings)) cli::cli_abort("proj_settings apply only when hindcast = TRUE.")
  if (!is.null(grad_tol) && (!is.numeric(grad_tol) || length(grad_tol) != 1 || !is.finite(grad_tol) || grad_tol <= 0)) cli::cli_abort("grad_tol must be positive and finite.")
  settings <- if (hindcast) .retro_projection(proj_settings) else NULL
  jobs <- .bind_rows(lapply(seq_along(models), function(i) {
    fit <- models[[i]]; final <- .terminal_year(fit)
    first <- min(fit$dat$years[!fit$dat$is_proj])
    if (!is.numeric(folds) || length(folds) != 1 || !is.finite(folds) || folds != trunc(folds) || folds < 0 || folds > final - first - 1) {
      cli::cli_abort("folds must be a non-negative integer leaving at least two historical years in every fit.")
    }
    data.frame(model_index = i, fold = seq.int(final - folds, final), final = final)
  }))
  results <- progressr::with_progress({
    tick <- progressr::progressor(steps = nrow(jobs))
    furrr::future_map(seq_len(nrow(jobs)), function(i) {
      fit <- models[[jobs$model_index[i]]]
      tolerance <- if (!is.null(grad_tol)) grad_tol else if (!is.null(fit$grad_tol)) fit$grad_tol else .01
      result <- .retro_fold(fit, jobs$fold[i], jobs$final[i], settings, start_from_fit, tolerance)
      tick(); result
    }, .options = furrr::furrr_options(seed = 1, packages = "tinyAM", globals = globals))
  }, enable = progress)
  out <- list(fits = list(), pop = list(), obs_pred = list(), fold_diagnostics = data.frame(),
    projection_controls = data.frame(), mohns_rho = data.frame(), scores = data.frame(),
    forecast_pairs = list(index = data.frame(), aggregate_index = data.frame(), catch = data.frame(), mean_age = data.frame()),
    hindcast = hindcast, proj_settings = settings)
  stacks <- list(pop = list(), obs_pred = list()); signatures <- list()
  for (i in seq_along(models)) {
    label <- names(models)[i]; fit <- models[[i]]
    attempt <- results[jobs$model_index == i]; terminals <- jobs$fold[jobs$model_index == i]
    diagnostics <- .bind_rows(lapply(attempt, `[[`, "diagnostic")); diagnostics$model <- label
    out$fold_diagnostics <- .bind_rows(list(out$fold_diagnostics, diagnostics))
    controls <- .bind_rows(lapply(seq_along(attempt), function(j) {
      d <- attempt[[j]]$controls
      if (nrow(d)) { d$model <- label; d$fold <- terminals[j] }; d
    }))
    out$projection_controls <- .bind_rows(list(out$projection_controls, controls))
    all_fits <- lapply(attempt, `[[`, "fit"); names(all_fits) <- terminals
    out$fits[[label]] <- all_fits
    retained <- all_fits[diagnostics$is_converged]
    tidy <- if (length(retained)) tidy_tam(model_list = retained, label = "fold") else list(pop = list(), obs_pred = list())
    for (section in names(stacks)) for (metric in names(tidy[[section]])) {
      d <- tidy[[section]][[metric]]
      if (!nrow(d)) next
      d$fold <- as.numeric(as.character(d$fold)); d$model <- label
      d$horizon <- ifelse(d$is_proj, d$year - d$fold, 0)
      stacks[[section]][[metric]] <- c(stacks[[section]][[metric]], list(d))
    }
    for (metric in names(tidy$pop)) {
      d <- tidy$pop[[metric]]
      if (!nrow(d)) next
      d$fold <- as.numeric(as.character(d$fold))
      groups <- if ("age" %in% names(d)) split(d, d$age) else list(d)
      rows <- lapply(groups, function(x) data.frame(model = label, metric = metric,
        age = if ("age" %in% names(x)) x$age[1] else NA_real_,
        rho = if (as.character(.terminal_year(fit)) %in% names(retained)) compute_mohns_rho(x) else NA_real_))
      out$mohns_rho <- .bind_rows(c(list(out$mohns_rho), rows))
    }
    pairs <- list(index = data.frame(), aggregate_index = data.frame(), catch = data.frame(), mean_age = data.frame())
    if (hindcast && length(retained)) for (type in c("index", "catch")) {
      p <- tidy$obs_pred[[type]]
      if (is.null(p) || !nrow(p)) next
      p$fold <- as.numeric(as.character(p$fold)); p$horizon <- ifelse(p$is_proj, p$year - p$fold, 0)
      matched <- .forecast_pairs(p, fit$obs_pred[[type]])
      scored <- if (type == "index") .index_score_tables(matched, fit$obs_pred$index) else .catch_compositions(matched, fit$obs_pred$catch, fit$dat$ages)
      pairs[[type]] <- scored$pairs
      if (type == "index") pairs$aggregate_index <- scored$aggregate else pairs$mean_age <- scored$mean_age
      s <- scored$scores
      if (nrow(s)) { s$type <- type; s$model <- label; out$scores <- .bind_rows(list(out$scores, s)) }
    }
    for (type in names(pairs)) {
      d <- pairs[[type]]
      if (nrow(d)) d$model <- label
      out$forecast_pairs[[type]] <- .bind_rows(list(out$forecast_pairs[[type]], d))
    }
    signatures[[label]] <- lapply(pairs, function(d) {
      if (!nrow(d)) return(data.frame())
      valid <- is.finite(d$obs) & is.finite(d$pred)
      fields <- intersect(c("fold", "horizon", "year", "age", "survey", "fleet", "obs", "baseline", "baseline_year", "unit"), names(d))
      signature <- d[valid, fields, drop = FALSE]
      if (nrow(signature)) signature <- signature[do.call(order, signature), , drop = FALSE]
      rownames(signature) <- NULL; signature
    })
  }
  for (section in names(stacks)) out[[section]] <- lapply(stacks[[section]], .bind_rows)
  out$comparable <- length(models) == 1 || all(vapply(signatures[-1], identical, logical(1), signatures[[1]]))
  if (any(!out$fold_diagnostics$is_converged)) cli::cli_warn("Some retrospective fits or projections failed. See fold_diagnostics; failed fits are excluded from scores.")
  if (hindcast && !out$comparable) cli::cli_warn("Models have different forecast observations, units or retained folds. Scores are not directly comparable; no ranking is provided.")
  class(out) <- c("tam_retro", "list"); out
}

#' Retrospective fits and observation hindcasts
#'
#' Refit models after peeling terminal years. Hindcasts additionally forecast
#' held-out observations, conditional on observed removals by default.
#' @param fit,... One or more `tam_fit` objects. Name additional models for distinct
#'   labels. All other controls must be named.
#' @param folds Number of peels, leaving at least two historical years per fit.
#' @param start_from_fit Use each supplied fit to initialize its peels?
#' @param hindcast Project peeled fits? The only default difference between functions.
#' @param proj_settings Hindcast settings: `n_proj` (default 1), `n_mean` (default 1),
#'   and either `catch` (default `"observed"`) or `F_mult`. Numeric controls are one
#'   value recycled, or `n_proj` values ordered by forecast horizon. Observed catch
#'   is complete raw catch numbers times supplied weights. Forecasts are shortened
#'   to available validation years. Ordinary retrospectives reject this argument.
#' @param grad_tol Numerical check tolerance; defaults to each fit's stored value.
#' @param progress Show progress using progressr?
#' @param globals Global names passed to furrr workers. Configure your own future plan.
#' @details Only observed years determine folds. Historical references have no
#'   forecasts. Held-out indices and catch compositions do not enter fitting;
#'   observed-catch conditioning uses only aggregate removals. Forecast biology
#'   uses training years and `n_mean`.
#'
#' Survey log-RMSE pools observation errors. MASE compares absolute log errors
#' with persistence on identical cases; see [compute_hindcast_mase()]. Overall
#' MASE averages ages within survey, surveys within horizon, then horizons.
#' Aggregate indices require complete, consistent age coverage. Catch proportions
#' and mean modeled age (plus age coded by its lower bound) are scored separately
#' on their natural scale; genuine zeros are retained. Overlapping forecasts are
#' distinct origin/horizon cases, not independent observations. These are
#' conditional point-forecast diagnostics, not probability scores.
#'
#' Compare equivalent observations and age definitions. No ranking is produced.
#' Missing baselines or zero persistence error make MASE unavailable. Failed
#' numerical fits remain in `fits`, but not in stacked results or scores.
#' Supplied fits and fitted objectives are unchanged.
#' @return A `tam_retro` list with `fits` nested by model and terminal year;
#'   stacked `pop` and `obs_pred` with model, fold and horizon (zero historically);
#'   `mohns_rho`; detailed `scores`; matched `forecast_pairs`; `fold_diagnostics`;
#'   and resolved `projection_controls`. Failed references make rho unavailable.
#'   All-failed runs still return diagnostics.
#' @examples
#' \dontrun{
#' retro <- fit_retro(fit, alternative = other_fit, folds = 5)
#' hind <- fit_hindcast(fit, folds = 10,
#'   proj_settings = list(n_proj = 3, n_mean = 3, catch = "observed"))
#' print(hind)
#' hind$scores
#' vis_tam(hind$fits[[1]])
#' }
#' @seealso [fit_tam()], [compute_mohns_rho()]
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr with_progress progressor
#' @importFrom stats update
#' @export
fit_retro <- function(fit, ..., folds = 2, start_from_fit = FALSE, hindcast = FALSE,
                      proj_settings = NULL, grad_tol = NULL, progress = TRUE, globals = NULL) {
  models <- list(fit, ...)
  names(models) <- .retro_labels(models, as.list(substitute(list(fit, ...)))[-1])
  .run_retros(models, folds, start_from_fit, hindcast, proj_settings, grad_tol, progress, globals)
}

#' @rdname fit_retro
#' @export
fit_hindcast <- function(fit, ..., folds = 2, start_from_fit = FALSE, hindcast = TRUE,
                         proj_settings = NULL, grad_tol = NULL, progress = TRUE, globals = NULL) {
  models <- list(fit, ...)
  names(models) <- .retro_labels(models, as.list(substitute(list(fit, ...)))[-1])
  .run_retros(models, folds, start_from_fit, hindcast, proj_settings, grad_tol, progress, globals)
}

#' @rdname fit_retro
#' @param x A `tam_retro` object.
#' @export
print.tam_retro <- function(x, ...) {
  cli::cli_h1(if (x$hindcast) "Retrospectives and hindcasts" else "Retrospectives")
  models <- names(x$fits); table <- data.frame(model = models)
  for (metric in c("ssb", "recruitment", "F_bar")) table[[paste0("rho_", metric)]] <- vapply(models, function(model) {
    d <- x$mohns_rho
    values <- d$rho[d$model == model & d$metric == metric & is.na(d$age)]
    if (length(values)) values[1] else NA_real_
  }, numeric(1))
  if (x$hindcast) for (stat in c("RMSE", "MASE", "n")) table[[if (stat == "RMSE") "log_RMSE" else stat]] <- vapply(models, function(model) {
    d <- x$scores
    if (!nrow(d)) return(NA_real_)
    d <- d[d$model == model & d$type == "index" & !is.na(d$level) & d$level == "overall" & is.na(d$horizon), ]
    if (nrow(d)) d[[stat]][1] else NA_real_
  }, numeric(1))
  table$scored_folds <- vapply(models, function(m) sum(x$fold_diagnostics$model == m & x$fold_diagnostics$n_scored_index > 0), integer(1))
  table$failed_folds <- vapply(models, function(m) sum(x$fold_diagnostics$model == m & !x$fold_diagnostics$is_converged), integer(1))
  print(table, row.names = FALSE, digits = 4)
  if (x$hindcast) cli::cli_inform(if (!nrow(x$scores)) "No eligible forecast scores are available." else if (x$comparable) "Scores use matching validation cases across models." else "Validation cases differ; scores are not directly comparable.")
  invisible(x)
}
