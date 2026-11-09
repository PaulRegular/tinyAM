# Optional assessment diagnostics for the spatial × mortality experiment.
# Source observation_error_sim.R and observation_error_fit.R first.
.oe_z_fit_observations <- function(input, sd_reference, seed) {
  series <- .oe_fit_series()
  stopifnot(is.data.frame(sd_reference), all(c("component", "age", "sd") %in% names(sd_reference)),
    !anyDuplicated(sd_reference[c("component", "age")]))
  actual <- lapply(names(series), function(type) {
    d <- as.data.frame(input$training[[type]])
    d <- d[order(d$year, d$age), ]
    s <- sd_reference[sd_reference$component == series[[type]], ]
    known <- is.finite(s$sd) & s$sd > 0
    if (!any(known)) stop("No positive calibration log SD for ", series[[type]], ".")
    d$obs_sd <- s$sd[match(d$age, s$age)]
    d$sd_fallback <- !is.finite(d$obs_sd) | d$obs_sd <= 0
    d$obs_sd[d$sd_fallback] <- stats::median(s$sd[known])
    # Apply the same positive-observation mask to both protocols.
    d$unsupported_target <- is.finite(d$obs) & d$obs > 0 &
      (!is.finite(d$truth) | d$truth <= 0)
    d$obs[d$unsupported_target] <- NA_real_
    d
  })
  names(actual) <- names(series)
  independent <- .observation_error_with_seed(seed, lapply(actual, function(d) {
    positive <- is.finite(d$obs) & d$obs > 0
    # The known target is a log-scale mean (arithmetic-scale median).
    error <- stats::rnorm(nrow(d), sd = d$obs_sd)
    d$obs[positive] <- exp(log(d$truth[positive]) + error[positive])
    stopifnot(all(is.finite(d$obs[positive]) & d$obs[positive] > 0))
    d
  }))
  list(SimSurvey = actual, Independent = independent)
}

.oe_z_fit_residuals <- function(d, component, stage, pred, sd) {
  err <- .oe_log_error(d$obs, pred)
  data.frame(year = d$year, age = d$age, component = component, stage = stage,
    obs = d$obs, target = d$truth, pred = pred, sd = sd, log_error = err,
    residual = err / sd, sd_fallback = d$sd_fallback,
    unsupported_target = d$unsupported_target)
}

.oe_z_fit_spec <- function(truth, training, model, heldout_year) {
  if (model == "constant_M") {
    # A nonphysical realized M in the generator should only fail the truth-M fit.
    # These replacements supply the template, not the target observations.
    truth$population <- as.data.frame(truth$population)
    truth$population$M <- 0.3
    truth$population$Z <- truth$population$F + 0.3
  }
  data <- .oe_fit_data(truth)$data
  data <- .oe_fit_inputs(data, training, training, NULL, "constant")$data
  for (type in names(training)) {
    keys <- transform(as.data.frame(data[[type]]), age = pmin(age, 14L))
    data[[type]]$obs_sd <- training[[type]]$obs_sd[
      match(.oe_fit_key(keys), .oe_fit_key(training[[type]]))]
  }
  if (!is.na(heldout_year)) data$index$obs[data$index$year == heldout_year] <- NA_real_
  settings <- .oe_fit_settings("", "constant")
  settings$F_settings$mu_form <- ~ factor(year) + double_logistic(age)
  settings$catch_settings <- list(sd_form = ~ 0, sd_supplied = ~ obs_sd)
  settings$index_settings <- list(q_form = ~ 0 + logistic(age),
    sd_form = ~ 0, sd_supplied = ~ obs_sd)
  if (model == "constant_M") settings$M_settings$mu_supplied <- ~ I(0.3)
  c(list(data = data, silent = TRUE), settings)
}

.oe_z_fit_once <- function(truth, training, model, heldout_year) {
  fit <- do.call(tinyAM::fit_tam, .oe_z_fit_spec(truth, training, model, heldout_year))
  checks <- tinyAM::check_tam(fit)
  if (!isTRUE(checks$is_converged)) {
    failed <- subset(checks$numerical, status != "pass")
    return(list(usable = FALSE,
      numerical = paste(paste(failed$check, failed$detail, sep = ": "), collapse = " | ")))
  }
  series <- .oe_fit_series()
  types <- if (is.na(heldout_year)) names(series) else "index"
  residuals <- lapply(types, function(type) {
    d <- training[[type]]
    fitted <- fit$obs_pred[[type]]
    pred_rows <- match(.oe_fit_key(d), .oe_fit_key(fitted))
    stopifnot(!anyNA(pred_rows), !anyDuplicated(.oe_fit_key(fitted)))
    if (is.na(heldout_year)) {
      ans <- .oe_z_fit_residuals(d, series[[type]], "fitted", fitted$pred[pred_rows],
        fitted$sd[pred_rows])
      ans$residual <- fitted$std_res[pred_rows]
      return(ans)
    }
    keep <- d$year == heldout_year
    pred_rows <- pred_rows[keep]
    stopifnot(any(keep), all(is.na(fitted$obs[pred_rows])))
    # The observation SD is known; parameter and state uncertainty are excluded.
    .oe_z_fit_residuals(d[keep, ], series[[type]], "heldout", fitted$pred[pred_rows],
      fitted$sd[pred_rows])
  })
  list(usable = TRUE, residuals = do.call(rbind, residuals))
}

oe_z_fit_history <- function(input, cache, sd_reference, seed, heldout_years = c(5L, 16L)) {
  observations <- .oe_z_fit_observations(input, sd_reference, seed)
  keys <- c("scenario", "spatial", "mortality", "history_id")
  meta <- c(as.list(input$scenario), as.list(input$history))
  stopifnot(all(keys %in% names(meta)), !anyDuplicated(heldout_years),
    all(heldout_years %in% input$training$index$year))
  annotate <- function(rows, observation, model) {
    if (is.null(rows)) return(NULL)
    for (key in keys) rows[[key]] <- meta[[key]]
    rows$observation <- observation
    rows$model <- model
    rows
  }
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  attempts <- list()
  residuals <- list()
  for (observation in names(observations)) {
    training <- observations[[observation]]
    for (type in names(training)) {
      d <- training[[type]]
      residuals[[length(residuals) + 1L]] <- annotate(
        .oe_z_fit_residuals(d, .oe_fit_series()[[type]], "oracle", d$truth, d$obs_sd),
        observation, "oracle")
    }
    for (model in c("truth_M", "constant_M")) {
      for (year in c(NA_integer_, heldout_years)) {
        fold <- if (is.na(year)) "full" else sprintf("leave-%02d", year)
        path <- file.path(cache, paste0("fit-", tolower(observation), "-", model, "-", fold, ".rds"))
        attempt <- checkpoint(path, function() .oe_z_fit_once(input$truth, training, model, year))
        fit <- attempt$value
        usable <- isTRUE(fit$usable)
        status <- if (usable) "success" else if (is.null(attempt$error)) "numerical failure" else "error"
        error <- if (is.null(attempt$error)) "" else attempt$error
        attempts[[length(attempts) + 1L]] <- annotate(data.frame(heldout_year = year,
          usable = usable, status = status, warnings = paste(attempt$warnings, collapse = " | "),
          error = error, numerical = paste(fit$numerical, collapse = " | "),
          elapsed = attempt$elapsed), observation, model)
        if (usable) residuals[[length(residuals) + 1L]] <- annotate(fit$residuals, observation, model)
      }
    }
  }
  list(attempts = do.call(rbind, attempts), residuals = do.call(rbind, residuals))
}

oe_z_residual_summary <- function(residuals, heldout_years = c(5L, 16L)) {
  if (!nrow(residuals)) return(data.frame())
  group_keys <- c("scenario", "spatial", "mortality", "history_id", "component")
  groups <- split(residuals, residuals[group_keys], drop = TRUE)
  summaries <- lapply(groups, function(group) {
    all_years <- group[group$stage != "heldout", ]
    scopes <- list(all_years = all_years)
    if (group$component[1L] == "RV") {
      scopes$heldout_years <- group[group$year %in% heldout_years, ]
    }
    do.call(rbind, lapply(names(scopes), function(scope) {
      d <- scopes[[scope]]
      series <- split(d, d[c("observation", "model", "stage")], drop = TRUE)
      if (!length(series)) return(NULL)
      # Match actual/control and available fitted models on finite year–age cells.
      keys <- Reduce(intersect, lapply(series, function(x) .oe_fit_key(x)[is.finite(x$residual)]))
      complete <- length(series) == if (scope == "all_years") 6L else 10L
      if (scope == "heldout_years") complete <- complete &&
        all(vapply(series, function(x) all(heldout_years %in% x$year), logical(1)))
      do.call(rbind, lapply(series, function(x) {
        xy <- x[match(keys, .oe_fit_key(x)), c("year", "age", "residual")]
        pairs <- merge(xy, xy, by = "year", suffixes = c("_a", "_b"))
        pairs <- pairs[pairs$age_b > pairs$age_a, ]
        by_gap <- split(pairs, pairs$age_b - pairs$age_a)
        if (!length(by_gap)) return(NULL)
        out <- do.call(rbind, lapply(names(by_gap), function(gap) {
          p <- by_gap[[gap]]
          # Descriptive Pearson association; age pairs within a year overlap.
          estimate <- if (nrow(p) > 2L && stats::sd(p$residual_a) > 0 && stats::sd(p$residual_b) > 0) {
            stats::cor(p$residual_a, p$residual_b)
          } else NA_real_
          data.frame(age_gap = as.integer(gap), estimate = estimate, n_pairs = nrow(p),
            n_years = length(unique(p$year)))
        }))
        cbind(x[rep(1L, nrow(out)), c(group_keys, "observation", "model", "stage")], out,
          scope = scope, complete = complete, n_series = length(series))
      }))
    }))
  })
  out <- do.call(rbind, summaries)
  if (is.null(out)) data.frame() else out
}

oe_z_assessments <- function(z, out, heldout_years = c(5L, 16L)) {
  inputs <- as.data.frame(z$assessment_inputs)
  scenarios <- unique(inputs$scenario)
  histories <- sort(unique(inputs$history_id))
  scales <- stats::setNames(lapply(scenarios, function(id) {
    tryCatch(oe_z_residual_scale(z$rows[z$rows$scenario == id, ]), error = identity)
  }), scenarios)
  unavailable <- function(case, problem) {
    d <- expand.grid(heldout_year = c(NA_integer_, heldout_years),
      observation = c("SimSurvey", "Independent"), model = c("truth_M", "constant_M"),
      stringsAsFactors = FALSE)
    for (key in c("scenario", "spatial", "mortality", "history_id")) d[[key]] <- case[[key]]
    d$usable <- FALSE
    d$status <- "input unavailable"
    d$warnings <- d$numerical <- ""
    d$error <- problem
    d$elapsed <- 0
    list(attempts = d, residuals = NULL)
  }
  runs <- lapply(seq_len(nrow(inputs)), function(i) {
    case <- inputs[i, ]
    if (!isTRUE(case$complete) || !file.exists(case$file)) {
      return(unavailable(case, "The simulated assessment history is incomplete."))
    }
    scale <- scales[[case$scenario]]
    if (inherits(scale, "error")) return(unavailable(case, conditionMessage(scale)))
    tryCatch(oe_z_fit_history(readRDS(case$file),
      file.path(out, "assessment-fits", case$scenario, case$history_id), scale,
      seed = 55100L + match(case$history_id, histories), heldout_years = heldout_years),
      error = function(e) unavailable(case, conditionMessage(e)))
  })
  attempts <- do.call(rbind, lapply(runs, `[[`, "attempts"))
  residuals <- do.call(rbind, lapply(runs, `[[`, "residuals"))
  if (is.null(attempts)) attempts <- data.frame()
  if (is.null(residuals)) residuals <- data.frame()
  list(attempts = attempts, residuals = residuals,
    correlation = oe_z_residual_summary(residuals, heldout_years), heldout_years = heldout_years)
}
