# Compact summaries only: no fitting and no changes to stored results.
delta_report_colours <- c(standard_ln = "#f05c28",
  standard_gg = "#b4461e", poisson_ln = "#028f8d", poisson_gg = "#735391",
  tweedie = "#a17400")
delta_report_keys <- c("experiment", "driver", "effort", "age", "mortality", "id", "family")

delta_report_table <- function(x, defaults) {
  if (is.null(x)) x <- data.frame()
  if (!is.data.frame(x)) stop("A report summary entry is not a data frame.")
  for (name in setdiff(names(defaults), names(x))) x[[name]] <- rep(defaults[[name]], nrow(x))
  x
}

delta_report_read <- function(path) {
  if (is.null(path) || length(path) != 1L || !file.exists(path)) {
    stop("Supply the current study's summary.rds as params$result.")
  }
  result <- readRDS(path)
  if (is.null(result$fits) && !is.null(result$stages)) {
    stop("This is a frozen earlier-study summary. Keep its original report; ",
      "render this report from the new flat summary without modifying old results.")
  }
  keys <- list(experiment = "unknown", driver = "none", effort = "standard",
    age = NA_integer_, mortality = "known", id = NA_integer_, family = "unknown")
  result$fits <- delta_report_table(result$fits, c(keys, list(converged = NA,
    error = "", warnings = "", failed_checks = "", elapsed = NA_real_)))
  result$recovery <- delta_report_table(result$recovery, c(keys, list(metric = "", window = "full",
    absolute = NA_real_, bias = NA_real_, coverage = NA_real_, estimate = NA_real_, truth = NA_real_)))
  result$scores <- delta_report_table(result$scores, c(keys, list(component = "unknown",
    age_group = "All ages", n = 0L, mean_zero = NA_real_, observed_zero = NA_real_,
    brier = NA_real_, n_positive = 0L, log_rmse = NA_real_, positive_coverage = NA_real_)))
  result$gates <- delta_report_table(result$gates, list())
  result$audits <- delta_report_table(result$audits, list())
  result$reference <- delta_report_table(result$reference, list(experiment = "unknown",
    driver = "none", id = NA_integer_, year = NA_integer_, age = NA_integer_,
    N = NA_real_, Z = NA_real_, F_realized = NA_real_, M_realized = NA_real_, ssb = NA_real_))
  result$profile <- delta_report_table(result$profile, c(keys,
    list(component = "unknown", n = 0L, mean_zero = NA_real_, observed_zero = NA_real_)))
  result$sampling <- delta_report_table(result$sampling, list())
  result$parameters <- delta_report_table(result$parameters, list())
  result$observations <- delta_report_table(result$observations, list(experiment = "unknown",
    driver = "none", effort = "standard", id = NA_integer_, component = "unknown",
    year = NA_integer_, age = NA_integer_, obs = NA_real_))
  for (table in c("fits", "recovery", "scores", "profile")) {
    result[[table]]$family <- as.character(result[[table]]$family)
    result[[table]]$experiment <- as.character(result[[table]]$experiment)
    result[[table]]$age <- as.integer(result[[table]]$age)
  }
  if (is.list(result$fits$warnings)) result$fits$warnings <- vapply(result$fits$warnings,
    function(x) paste(unique(as.character(x)), collapse = "; "), character(1L))
  result$fits$warnings[is.na(result$fits$warnings)] <- ""
  result$fits$error[is.na(result$fits$error)] <- ""
  result$fits$status <- ifelse(is.na(result$fits$converged), "Pending",
    ifelse(result$fits$converged, "Passed", "Failed"))
  if (anyDuplicated(delta_report_key(result$fits))) stop("Duplicate fit identifiers; refusing to count a history twice.")
  result
}

delta_report_key <- function(x, keys = delta_report_keys) {
  if (!nrow(x)) return(character())
  do.call(paste, c(lapply(x[keys], function(y) ifelse(is.na(y), "<NA>", as.character(y))), sep = "\034"))
}

delta_report_passed <- function(x, fits) {
  passed <- delta_report_key(fits[fits$status == "Passed", , drop = FALSE])
  x[delta_report_key(x) %in% passed, , drop = FALSE]
}

delta_report_primary <- function(x) {
  x$experiment == "spatial" & x$effort == "standard" & !is.na(x$age) &
    x$age == 20L & x$mortality %in% c("known", "estimated")
}

delta_report_analysis_scope <- function(x) {
  # The stopped high-F AR1-M trial is retained in fit/convergence records,
  # but its partial family coverage does not define a performance comparison.
  !(x$experiment == "spatial" & x$mortality == "estimated" & x$driver == "high_f")
}

delta_report_interval <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  centre <- if (n) mean(x) else NA_real_
  half <- if (n >= 2L) stats::qt(0.975, n - 1L) * stats::sd(x) / sqrt(n) else NA_real_
  c(n = n, mean = centre, lower = centre - half, upper = centre + half)
}

delta_report_summarize <- function(data, value, extra = character()) {
  if ("window" %in% names(data)) extra <- unique(c(extra, "window"))
  keys <- c(setdiff(delta_report_keys, "id"), extra)
  if (!nrow(data)) {
    out <- data[FALSE, keys, drop = FALSE]
    out$n_histories <- integer()
    out$mean <- out$lower <- out$upper <- numeric()
    return(out)
  }
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(n_histories = delta_report_interval(.data[[value]])[["n"]],
      mean = delta_report_interval(.data[[value]])[["mean"]],
      lower = delta_report_interval(.data[[value]])[["lower"]],
      upper = delta_report_interval(.data[[value]])[["upper"]], .groups = "drop")
}

delta_report_pairs <- function(data, value, baseline = "lognormal", extra = character(),
    count = NULL, support = NULL) {
  if ("window" %in% names(data)) extra <- unique(c(extra, "window"))
  keys <- c(setdiff(delta_report_keys, "family"), extra)
  selected <- unique(c(keys, "family", value, count, support))
  alternatives <- data[data$family != baseline, selected, drop = FALSE]
  reference <- data[data$family == baseline, setdiff(selected, "family"), drop = FALSE]
  if (anyDuplicated(delta_report_key(alternatives, c(keys, "family"))) ||
      anyDuplicated(delta_report_key(reference, keys))) stop("Duplicate scoring identifiers; refusing non-independent pairs.")
  joined <- merge(alternatives, reference, by = keys, suffixes = c("_arm", "_base"))
  joined$delta <- joined[[paste0(value, "_arm")]] - joined[[paste0(value, "_base")]]
  joined$matched_support <- rep(TRUE, nrow(joined))
  for (column in c(count, support)) {
    arm <- joined[[paste0(column, "_arm")]]
    base <- joined[[paste0(column, "_base")]]
    joined$matched_support <- joined$matched_support & !is.na(arm) & !is.na(base) & arm == base
    if (identical(column, count)) joined$matched_support <- joined$matched_support & arm > 0
  }
  omitted <- sum(!joined$matched_support)
  joined <- joined[joined$matched_support & is.finite(joined$delta), , drop = FALSE]
  summary <- delta_report_summarize(joined, "delta", extra)
  names(summary)[names(summary) == "mean"] <- "change"
  summary$baseline <- rep(baseline, nrow(summary))
  list(summary = summary, values = joined, unmatched_support = omitted)
}

delta_report_mortality_pairs <- function(recovery) {
  d <- recovery[recovery$experiment == "spatial" & recovery$driver == "low_f" &
    recovery$effort == "standard" & recovery$age == 20L &
    recovery$mortality %in% c("known", "estimated") &
    recovery$metric %in% c("Age 14+ abundance", "SSB", "Fishing mortality",
      "Natural mortality", "Total mortality"), , drop = FALSE]
  keys <- c(setdiff(delta_report_keys, "mortality"), "metric", "window")
  ar1 <- d[d$mortality == "estimated", c(keys, "absolute"), drop = FALSE]
  fixed <- d[d$mortality == "known", c(keys, "absolute"), drop = FALSE]
  if (anyDuplicated(delta_report_key(ar1, keys)) || anyDuplicated(delta_report_key(fixed, keys))) {
    stop("Duplicate mortality-comparison identifiers; refusing non-independent pairs.")
  }
  joined <- merge(ar1, fixed, by = keys, suffixes = c("_ar1", "_fixed"))
  joined$delta <- joined$absolute_ar1 - joined$absolute_fixed
  joined <- joined[is.finite(joined$delta), , drop = FALSE]
  joined$mortality <- rep("estimated", nrow(joined))
  summary <- delta_report_summarize(joined, "delta", extra = c("metric", "window"))
  names(summary)[names(summary) == "mean"] <- "change"
  summary$baseline <- rep("Fixed mean M", nrow(summary))
  list(summary = summary, values = joined)
}

delta_report_runtime <- function(fits) {
  if (!nrow(fits)) return("No fit attempts are recorded yet.")
  completed <- fits$status != "Pending"
  seconds <- sum(fits$elapsed[completed & is.finite(fits$elapsed)])
  pending <- sum(fits$status == "Pending")
  sprintf("%d of %d completed attempts passed; %d additional %s marked pending. Recorded fitting time totals %.1f worker-hours.",
    sum(fits$status == "Passed"), sum(completed), pending,
    if (pending == 1L) "entry is" else "entries are", seconds / 3600)
}

delta_report_zero_figure <- function(observations) {
  by_history <- observations |>
    dplyr::filter(experiment == "spatial", !is.na(obs)) |>
    dplyr::mutate(program = dplyr::case_when(
      component == "index" ~ "RV survey (both shares)",
      driver == "low_f" ~ "Catch (10% fishing share)",
      driver == "high_f" ~ "Catch (90% fishing share)"
    )) |>
    # Both fishing shares use the same survey draw.
    dplyr::distinct(id, effort, program, year, age, obs) |>
    dplyr::group_by(id, effort, program) |>
    dplyr::summarise(fraction = mean(obs == 0), .groups = "drop")

  means <- by_history |>
    dplyr::group_by(effort, program) |>
    dplyr::summarise(interval = list(delta_report_interval(fraction)), .groups = "drop")
  means <- cbind(means[c("effort", "program")], as.data.frame(do.call(rbind, means$interval)))
  means$effort <- factor(means$effort, levels = c("low", "standard", "high"),
    labels = c("Low", "Standard", "High"))
  colours <- c("RV survey (both shares)" = "#003e56", "Catch (10% fishing share)" = "#b4461e",
    "Catch (90% fishing share)" = "#007976")
  means$program <- factor(means$program, levels = names(colours))

  ggplot2::ggplot(means, ggplot2::aes(effort, mean, group = program, colour = program, shape = program)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.07, linewidth = 0.65) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = "Monitoring effort", y = "Age-by-year observations recorded as zero",
      colour = NULL, shape = NULL) +
    ggplot2::theme(legend.position = "top", panel.grid.major.x = ggplot2::element_blank())
}

delta_report_outcome_figure <- function(positive_scores, older_recovery) {
  positives <- positive_scores |>
    dplyr::filter(age_group == "All ages") |>
    dplyr::transmute(dplyr::across(dplyr::all_of(delta_report_keys)),
      outcome = ifelse(component == "index", "Positive survey index\nlog-RMSE", "Positive catch\nlog-RMSE"),
      error = log_rmse, observations = n_positive)
  older <- older_recovery |>
    dplyr::filter(window == "full", metric == "Age 14+ abundance") |>
    dplyr::transmute(dplyr::across(dplyr::all_of(delta_report_keys)),
      outcome = "Age 14+ abundance\nabsolute log error", error = absolute)
  positive_changes <- delta_report_pairs(positives, "error", extra = "outcome", count = "observations")
  older_changes <- delta_report_pairs(older, "error", extra = "outcome")
  # The descriptive baseline uses every eligible lognormal history once.
  # Paired changes below still use the histories shared with each alternative.
  baseline <- dplyr::bind_rows(positives, older) |>
    dplyr::filter(family == "lognormal", is.finite(error)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(setdiff(delta_report_keys, "id"))), outcome) |>
    dplyr::summarise(n_histories = dplyr::n(), score = mean(error), change = 0, .groups = "drop")
  settings <- c("90% fishing\nFixed mean M", "10% fishing\nFixed mean M", "10% fishing\nVariable M")
  families <- c(lognormal = "Ordinary lognormal (reference)",
    standard_ln = "Logit-link delta lognormal", standard_gg = "Logit-link delta gen. gamma",
    poisson_ln = "Poisson-link delta lognormal", poisson_gg = "Poisson-link delta gen. gamma",
    tweedie = "Tweedie")
  data <- dplyr::bind_rows(baseline, positive_changes$summary, older_changes$summary) |>
    dplyr::filter(experiment == "spatial", age == 20L, effort == "standard") |>
    dplyr::mutate(setting = dplyr::case_when(
      driver == "high_f" & mortality == "known" ~ settings[1],
      driver == "low_f" & mortality == "known" ~ settings[2],
      driver == "low_f" & mortality == "estimated" ~ settings[3]
    )) |>
    dplyr::filter(!is.na(setting), family %in% names(families))
  data$setting <- factor(data$setting, levels = settings)
  data$outcome <- factor(data$outcome, levels = c("Positive survey index\nlog-RMSE",
    "Positive catch\nlog-RMSE", "Age 14+ abundance\nabsolute log error"))

  # Use one label per family and setting. A range covers unequal scoring support.
  rows <- data |>
    dplyr::group_by(setting, family) |>
    dplyr::summarise(count = paste(unique(range(n_histories)), collapse = "-"), .groups = "drop")
  rows$key <- paste(rows$setting, rows$family)
  rows$label <- paste0(families[rows$family], "  (", rows$count, ")")
  rows <- rows[order(rows$setting, match(rows$family, rev(names(families)))), ]
  data$row <- factor(paste(data$setting, data$family), levels = rows$key)
  reference <- data[data$family == "lognormal", ]
  alternatives <- data[data$family != "lognormal", ]

  ggplot2::ggplot(data, ggplot2::aes(change, row, colour = family)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45", linetype = "dashed") +
    ggplot2::geom_errorbar(data = alternatives, ggplot2::aes(xmin = lower, xmax = upper),
      orientation = "y", width = 0.17, linewidth = 0.65) +
    ggplot2::geom_point(data = alternatives, size = 2.5) +
    ggplot2::geom_point(data = reference, shape = 18, size = 3.4) +
    ggplot2::geom_text(data = reference, ggplot2::aes(label = sprintf("%.2f", score)),
      hjust = 1.35, size = 3.3, fontface = "bold") +
    ggplot2::facet_grid(rows = ggplot2::vars(setting), cols = ggplot2::vars(outcome), scales = "free") +
    ggplot2::scale_y_discrete(labels = stats::setNames(rows$label, rows$key)) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    ggplot2::scale_colour_manual(values = c(lognormal = "#003e56", delta_report_colours), guide = "none") +
    ggplot2::labs(x = "Change in error from ordinary lognormal (lower is better)", y = NULL,
      subtitle = "Diamonds mark the reference. The number beside each diamond is the lognormal model's mean error.") +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0), strip.text = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 9), panel.grid.major.y = ggplot2::element_blank())
}
