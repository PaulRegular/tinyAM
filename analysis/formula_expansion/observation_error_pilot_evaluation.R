# Summarize the stock, rather than each year, as the independent replicate.
.oe_pilot_fit_keys <- c("population_id", "regime", "draw", "model")
.oe_pilot_regimes <- c("stable", "decline", "rise_fall")
.oe_pilot_recovery_metrics <- data.frame(
  metric = c("error_gain_pp", "terminal_gain_pp", "coverage_change_pp",
    "bias_change_pp", "terminal_bias_change_pp", "width_change_pp",
    "terminal_coverage_change_pp"),
  field = c("mean_abs_pct_error", "terminal_abs_pct_error", "coverage_pct",
    "mean_pct_bias", "terminal_pct_error", "mean_relative_width_pct",
    "terminal_coverage_pct"),
  direction = c(-1, -1, 1, 1, 1, 1, 1)
)
.oe_pilot_holdout_metrics <- data.frame(
  metric = c("nll_gain", "log_rmse_gain", "residual_coverage_change_pp"),
  field = c("nll", "log_rmse", "residual_coverage_pct"),
  direction = c(-1, -1, 1)
)

.oe_pilot_attach_fit <- function(d, fits, row_keys) {
  if (is.null(d) || nrow(d) == 0L) return(data.frame())
  stopifnot(!anyNA(d[c(.oe_pilot_fit_keys, row_keys)]),
    !anyDuplicated(d[c(.oe_pilot_fit_keys, row_keys)]))
  keys <- fits[.oe_pilot_fit_keys]
  keys$.fit_id <- seq_len(nrow(fits))
  linked <- merge(d, keys, by = .oe_pilot_fit_keys, sort = FALSE)
  stopifnot(nrow(linked) == nrow(d))
  linked
}

.oe_pilot_stock <- function(fits, annual, years) {
  quantities <- c("ssb", "recruitment", "F_bar")
  annual <- .oe_pilot_attach_fit(annual, fits, c("quantity", "year"))
  blocks <- if (nrow(annual)) {
    stopifnot(all(annual$quantity %in% quantities), all(annual$year %in% years))
    truth_blocks <- split(annual$truth, interaction(
      annual[c("population_id", "quantity", "year")], drop = TRUE))
    stopifnot(all(vapply(truth_blocks, function(x) {
      length(unique(x[is.finite(x)])) <= 1L
    }, logical(1))))
    split(annual, paste(annual$.fit_id, annual$quantity))
  }
  fits$.fit_id <- seq_len(nrow(fits))
  stock <- merge(fits[c(.oe_pilot_fit_keys, ".fit_id", "usable")],
    data.frame(quantity = quantities))
  names(stock)[names(stock) == "usable"] <- "fit_usable"
  fields <- .oe_pilot_recovery_metrics$field
  for (field in fields) stock[[field]] <- NA_real_
  stock$n_years <- 0L
  stock$eligible <- FALSE
  for (i in seq_len(nrow(stock))) {
    d <- blocks[[paste(stock$.fit_id[i], stock$quantity[i])]]
    if (is.null(d)) next
    stock$n_years[i] <- nrow(d)
    if (!stock$fit_usable[i] || stock$n_years[i] != length(years)) next
    d <- d[match(years, d$year), ]
    if (!all(is.finite(d$est) & is.finite(d$truth) & d$truth > 0 &
        is.finite(d$lwr) & is.finite(d$upr) & d$lwr <= d$upr)) next
    error <- 100 * (d$est / d$truth - 1)
    covered <- d$lwr <= d$truth & d$truth <= d$upr
    width <- 100 * (d$upr - d$lwr) / d$truth
    stock[i, fields] <- list(mean(abs(error)), abs(tail(error, 1L)),
      100 * mean(covered), mean(error), tail(error, 1L), mean(width),
      100 * tail(covered, 1L))
    stock$eligible[i] <- TRUE
  }
  stock$.fit_id <- NULL
  stock
}

.oe_pilot_pairs <- function(d, id, metrics, extra = character(), baseline_model = "constant",
                            candidate_models = setdiff(d$model, baseline_model)) {
  fields <- c("eligible", metrics$field, extra)
  baseline <- d[d$model == baseline_model, c(id, fields), drop = FALSE]
  candidate <- d[d$model %in% candidate_models, c(id, "model", fields), drop = FALSE]
  pairs <- merge(candidate, baseline, by = id, all.x = TRUE, sort = FALSE,
    suffixes = c("_candidate", "_baseline"))
  pairs$usable_pair <- (pairs$eligible_candidate %in% TRUE) &
    (pairs$eligible_baseline %in% TRUE)
  for (i in seq_len(nrow(metrics))) {
    field <- metrics$field[i]
    delta <- pairs[[paste0(field, "_candidate")]] - pairs[[paste0(field, "_baseline")]]
    pairs[[metrics$metric[i]]] <- ifelse(pairs$usable_pair, metrics$direction[i] * delta, NA_real_)
  }
  pairs
}

.oe_pilot_comparisons <- function(pairs, by, metrics) {
  groups <- split(pairs, interaction(pairs[by], drop = TRUE))
  if (!length(groups)) return(data.frame())
  do.call(rbind, lapply(groups, function(group) {
    label <- group[1L, by, drop = FALSE]
    paired <- group[group$usable_pair, , drop = FALSE]
    do.call(rbind, lapply(c("overall", .oe_pilot_regimes), function(regime) {
      regimes <- if (regime == "overall") .oe_pilot_regimes else regime
      cells <- split(paired, factor(paired$regime, levels = regimes))
      n <- vapply(cells, nrow, integer(1))
      do.call(rbind, lapply(seq_len(nrow(metrics)), function(i) {
        metric <- metrics$metric[i]
        field <- metrics$field[i]
        values <- lapply(cells, "[[", metric)
        cell_mean <- function(name) {
          if (any(n == 0L)) return(NA_real_)
          mean(vapply(cells, function(d) mean(d[[name]]), numeric(1)))
        }
        estimate <- cell_mean(metric)
        lower <- upper <- mcse <- NA_real_
        if (all(n >= 2L)) {
          contributions <- vapply(values, stats::var, numeric(1)) / n / length(regimes)^2
          variance <- sum(contributions)
          mcse <- sqrt(variance)
          # Satterthwaite degrees of freedom for independent regime means.
          half_width <- if (variance == 0) 0 else {
            df <- variance^2 / sum(contributions^2 / (n - 1L))
            stats::qt(.975, df) * mcse
          }
          lower <- estimate - half_width
          upper <- estimate + half_width
        }
        cbind(label, regime = regime, metric = metric, n_pairs = sum(n),
          estimate = estimate, lower = lower, upper = upper, mcse = mcse,
          baseline = cell_mean(paste0(field, "_baseline")),
          candidate = cell_mean(paste0(field, "_candidate")))
      }))
    }))
  }))
}

.oe_pilot_direct_comparisons <- function(stock) {
  contrasts <- data.frame(
    model = c("gillnet_smooth3", "gillnet_smooth3", "rv_smooth3", "rv_design_annual"),
    reference_model = c("gillnet_groups", "gillnet_age", "rv_age", "rv_pooled")
  )
  d <- stock[stock$draw == "primary" & stock$quantity == "ssb", ]
  metrics <- .oe_pilot_recovery_metrics[.oe_pilot_recovery_metrics$metric %in%
    c("error_gain_pp", "coverage_change_pp"), ]
  id <- c("population_id", "regime", "draw", "quantity")
  available <- contrasts$model %in% d$model & contrasts$reference_model %in% d$model
  pairs <- do.call(rbind, lapply(which(available), function(i) {
    reference <- contrasts$reference_model[i]
    p <- .oe_pilot_pairs(d, id, metrics, baseline_model = reference,
      candidate_models = contrasts$model[i])
    transform(p, reference_model = reference)
  }))
  if (!NROW(pairs)) return(data.frame())
  .oe_pilot_comparisons(pairs, c("reference_model", "model", "quantity"), metrics)
}

.oe_pilot_sample_size <- function(pairs, comparisons) {
  ssb <- pairs[pairs$quantity == "ssb", ]
  primary <- comparisons[comparisons$regime == "overall" & comparisons$quantity == "ssb", ]
  do.call(rbind, lapply(split(ssb, ssb$model), function(d) {
    by_class <- split(d, factor(d$regime, levels = .oe_pilot_regimes))
    usable <- lapply(by_class, function(x) x[x$usable_pair, , drop = FALSE])
    counts <- vapply(usable, nrow, integer(1))
    planned <- vapply(by_class, nrow, integer(1))
    paired_success <- sum(counts) / sum(planned)
    result <- data.frame(model = d$model[1L], n_pairs = sum(counts),
      smallest_class_n = min(counts), paired_success = paired_success,
      status = "Investigate fit reliability", error_usable_total = NA_integer_,
      coverage_usable_total = NA_integer_, balanced_usable_total = NA_integer_,
      balanced_scheduled_total = NA_integer_)
    if (sum(counts) < 20L || any(counts < 6L)) return(result)

    # Equal weighting and balanced allocation imply Var(mean) = sum(s_r^2)/(9*n).
    # The planning floor provides at least 20 balanced usable pairs.
    per_regime <- function(metric, half_width) {
      variances <- vapply(usable, function(x) stats::var(x[[metric]]), numeric(1))
      max(ceiling(20 / length(counts)),
        ceiling((stats::qnorm(.975) / half_width)^2 * sum(variances) / length(counts)^2))
    }
    error_n <- per_regime("error_gain_pp", 1)
    coverage_n <- per_regime("coverage_change_pp", 5)
    balanced_n <- max(error_n, coverage_n)
    result$error_usable_total <- length(counts) * error_n
    result$coverage_usable_total <- length(counts) * coverage_n
    result$balanced_usable_total <- length(counts) * balanced_n
    # Allow for paired fit failures at the lowest observed regime success rate.
    result$balanced_scheduled_total <- length(counts) * ceiling(balanced_n / min(counts / planned))

    change <- primary[primary$model == result$model, ]
    error <- change[change$metric == "error_gain_pp", ]
    coverage <- change[change$metric == "coverage_change_pp", ]
    result$status <- if (coverage$upper < -5) "Investigate coverage loss" else
      if (error$upper < 2) "Retain constant SD" else
      if (error$estimate >= 2 && error$lower > 0 && coverage$lower >= -5) {
        "Advance to confirmation"
      } else "More populations needed"
    result
  }))
}

.oe_pilot_holdout <- function(fits, holdout) {
  d <- .oe_pilot_attach_fit(holdout, fits, c("component", "year", "age"))
  components <- c("Gillnet", "RV")
  if (nrow(d)) stopifnot(all(d$component %in% components))
  blocks <- if (nrow(d)) split(d, paste(d$.fit_id, d$component))
  fits$.fit_id <- seq_len(nrow(fits))
  scores <- merge(fits[c(.oe_pilot_fit_keys, ".fit_id", "usable")],
    data.frame(component = components))
  scores$eligible <- FALSE
  scores$n_cells <- 0L
  scores$nll <- scores$log_rmse <- scores$residual_coverage_pct <- NA_real_
  scores$support <- rep(list(NULL), nrow(scores))
  for (i in seq_len(nrow(scores))) {
    d <- blocks[[paste(scores$.fit_id[i], scores$component[i])]]
    if (is.null(d)) next
    d <- d[is.finite(d$holdout_obs) & d$holdout_obs > 0, ]
    d <- d[order(d$year, d$age), ]
    scores$n_cells[i] <- nrow(d)
    support <- d[c("year", "age", "holdout_obs")]
    rownames(support) <- NULL
    scores$support[i] <- list(support)
    if (!scores$usable[i] || !nrow(d) ||
        !all(is.finite(d$holdout_log_error) & is.finite(d$sd) & d$sd > 0)) next
    nll <- -stats::dnorm(d$holdout_log_error, sd = d$sd, log = TRUE)
    if (!all(is.finite(nll))) next
    scores$nll[i] <- mean(nll)
    scores$log_rmse[i] <- sqrt(mean(d$holdout_log_error^2))
    scores$residual_coverage_pct[i] <- 100 * mean(abs(d$holdout_log_error) <= stats::qnorm(.975) * d$sd)
    scores$eligible[i] <- TRUE
  }
  pairs <- .oe_pilot_pairs(scores, c("population_id", "regime", "draw", "component"),
    .oe_pilot_holdout_metrics, extra = c("n_cells", "support"))
  pairs$same_support <- vapply(seq_len(nrow(pairs)), function(i) {
    identical(pairs$support_candidate[[i]], pairs$support_baseline[[i]])
  }, logical(1))
  pairs$usable_pair <- pairs$usable_pair & pairs$same_support
  pairs$n_cells <- ifelse(pairs$same_support, pairs$n_cells_baseline, NA_integer_)
  pairs[!pairs$usable_pair, .oe_pilot_holdout_metrics$metric] <- NA_real_
  pairs[c("support_candidate", "support_baseline")] <- NULL
  pairs
}

.oe_pilot_sd_cells <- function(validation, fits) {
  keys <- c("population_id", "component", "year", "age")
  stopifnot(all(c(keys, "replicate", "obs", "truth") %in% names(validation)),
    !anyNA(validation[c(keys, "replicate")]),
    !anyDuplicated(validation[c(keys, "replicate")]),
    all(validation$component %in% c("Gillnet", "RV")))
  cells <- do.call(rbind, lapply(split(validation, interaction(validation[keys], drop = TRUE)),
    function(d) {
      truth <- unique(d$truth)
      stopifnot(length(truth) == 1L)
      values <- d$obs[is.finite(d$obs) & d$obs >= 0]
      positive <- values[values > 0]
      n_positive <- length(positive)
      supported <- n_positive >= 5L && is.finite(truth) && truth > 0
      empirical_sd <- if (supported) stats::sd(log(positive) - log(truth)) else NA_real_
      arithmetic_sd <- stats::sd(values)
      cbind(d[1L, keys, drop = FALSE], truth = truth, n_draws = nrow(d),
        n_finite = length(values), n_positive = n_positive,
        n_zero = sum(values == 0), zero_fraction = mean(values == 0),
        empirical_sd = empirical_sd, arithmetic_sd = arithmetic_sd,
        arithmetic_cv = if (length(positive)) arithmetic_sd / mean(values) else NA_real_,
        relative_mcse = if (is.finite(empirical_sd) && empirical_sd > 0) {
          1 / sqrt(2 * (n_positive - 1L))
        } else NA_real_)
    }))
  metadata <- unique(fits[c("population_id", "regime")])
  cells <- merge(cells, metadata, by = "population_id", all.x = TRUE)
  stopifnot(!anyNA(cells$regime))
  cells
}

.oe_pilot_sd_profiles <- function(d) {
  keys <- c("population_id", "regime", "component", "model", "age")
  do.call(rbind, lapply(split(d, d$age), function(a) {
    variance <- a$empirical_sd^2
    empirical_sd <- sqrt(mean(variance))
    fitted_sd <- sqrt(mean(a$fitted_sd^2))
    baseline_sd <- sqrt(mean(a$baseline_sd^2))
    # Delta-method precision for the RMS of independent year-level SD estimates.
    relative_mcse <- if (empirical_sd > 0) {
      sqrt(sum(2 * variance^2 / (a$n_positive - 1L))) / (2 * sum(variance))
    } else NA_real_
    cbind(a[1L, keys, drop = FALSE], n_years = nrow(a),
      n_positive_min = min(a$n_positive), empirical_sd = empirical_sd,
      fitted_sd = fitted_sd, baseline_sd = baseline_sd, relative_mcse = relative_mcse,
      log_ratio = if (empirical_sd > 0) log(fitted_sd / empirical_sd) else NA_real_,
      log_ratio_baseline = if (empirical_sd > 0) log(baseline_sd / empirical_sd) else NA_real_)
  }))
}

.oe_pilot_sd_validation <- function(fits, holdout, validation) {
  if (!NROW(validation)) return(NULL)
  cells <- .oe_pilot_sd_cells(validation, fits)
  reliability <- do.call(rbind, lapply(split(cells, cells$component), function(d) {
    supported <- is.finite(d$empirical_sd)
    data.frame(component = d$component[1L], n_cells = nrow(d), n_supported = sum(supported),
      median_draws = stats::median(d$n_draws), median_positive = stats::median(d$n_positive),
      zero_fraction = sum(d$n_zero) / sum(d$n_finite),
      median_relative_mcse = stats::median(d$relative_mcse, na.rm = TRUE))
  }))
  answer <- list(cells = cells, reliability = reliability, profiles = data.frame(),
    pairs = data.frame(), comparisons = data.frame())
  if (!NROW(holdout)) return(answer)
  keys <- c("population_id", "component", "year", "age")
  stopifnot(all(c(keys, "draw", "model", "obs", "sd", "truth") %in% names(holdout)))

  # Compare candidate and baseline on the same positive training cells and
  # independent Monte Carlo draws. The empirical SD is conditional on positives.
  primary <- holdout[holdout$draw == "primary" & is.finite(holdout$obs) & holdout$obs > 0 &
    is.finite(holdout$sd) & holdout$sd > 0, ]
  primary <- merge(primary, fits[fits$usable, .oe_pilot_fit_keys], by = .oe_pilot_fit_keys)
  baseline <- primary[primary$model == "constant", c(keys, "obs", "sd", "truth")]
  names(baseline)[match(c("obs", "sd", "truth"), names(baseline))] <-
    c("baseline_obs", "baseline_sd", "baseline_truth")
  candidate <- primary[(startsWith(primary$model, "gillnet_") & primary$component == "Gillnet") |
    (startsWith(primary$model, "rv_") & primary$component == "RV"), c(keys, "model", "obs", "sd", "truth")]
  names(candidate)[match(c("sd", "truth"), names(candidate))] <- c("fitted_sd", "fitted_truth")
  matched <- merge(candidate, baseline, by = keys)
  matched <- merge(matched, cells[is.finite(cells$empirical_sd), ], by = keys)
  if (!nrow(matched)) return(answer)
  stopifnot(isTRUE(all.equal(matched$obs, matched$baseline_obs, tolerance = 0)),
    isTRUE(all.equal(matched$truth, matched$baseline_truth, tolerance = 1e-12)),
    isTRUE(all.equal(matched$truth, matched$fitted_truth, tolerance = 1e-12)))

  metrics <- data.frame(metric = c("sd_nll_gain", "abs_log_sd_gain"),
    field = c("sd_nll", "abs_log_sd"))
  id <- c("population_id", "regime", "model", "component")
  scores <- lapply(split(matched, interaction(matched[id], drop = TRUE)), function(d) {
    profile <- .oe_pilot_sd_profiles(d)
    loss <- function(sigma) {
      # Expectation under the empirical centered log-error distribution.
      # Unbiased sample variance estimates its second centered moment.
      log(sigma) + 0.5 * log(2 * pi) + 0.5 * (d$empirical_sd / sigma)^2
    }
    pairs <- cbind(d[1L, id, drop = FALSE], n_cells = nrow(d),
      n_ages = sum(is.finite(profile$log_ratio)),
      sd_nll_baseline = mean(loss(d$baseline_sd)), sd_nll_candidate = mean(loss(d$fitted_sd)),
      abs_log_sd_baseline = mean(abs(profile$log_ratio_baseline), na.rm = TRUE),
      abs_log_sd_candidate = mean(abs(profile$log_ratio), na.rm = TRUE))
    list(profile = profile, pairs = pairs)
  })
  answer$profiles <- do.call(rbind, lapply(scores, "[[", "profile"))
  pairs <- do.call(rbind, lapply(scores, "[[", "pairs"))
  pairs$usable_pair <- is.finite(pairs$sd_nll_baseline) & is.finite(pairs$sd_nll_candidate) &
    is.finite(pairs$abs_log_sd_baseline) & is.finite(pairs$abs_log_sd_candidate)
  pairs$sd_nll_gain <- pairs$sd_nll_baseline - pairs$sd_nll_candidate
  pairs$abs_log_sd_gain <- pairs$abs_log_sd_baseline - pairs$abs_log_sd_candidate
  answer$pairs <- pairs
  answer$comparisons <- .oe_pilot_comparisons(pairs, c("model", "component"), metrics)
  answer
}

oe_pilot_evaluate <- function(fits, annual, holdout = NULL, years = 1:20,
                              sd_validation = NULL) {
  stopifnot(nrow(fits) > 0L, length(years) > 0L, !anyNA(years),
    !is.unsorted(years, strictly = TRUE), !anyNA(fits[c(.oe_pilot_fit_keys, "usable")]),
    is.logical(fits$usable), !anyDuplicated(fits[.oe_pilot_fit_keys]),
    !anyDuplicated(unique(fits[c("population_id", "regime")])$population_id),
    all(fits$regime %in% .oe_pilot_regimes), all(fits$draw %in% c("primary", "repeat")))
  stock <- .oe_pilot_stock(fits, annual, years)
  pairs <- .oe_pilot_pairs(stock, c("population_id", "regime", "draw", "quantity"),
    .oe_pilot_recovery_metrics)
  repeat_pairs <- pairs[pairs$draw == "repeat", , drop = FALSE]
  pairs <- pairs[pairs$draw == "primary", , drop = FALSE]
  comparisons <- .oe_pilot_comparisons(pairs, c("model", "quantity"),
    .oe_pilot_recovery_metrics)

  ssb <- stock[stock$quantity == "ssb", ]
  failure_counts <- function(d) {
    stats::aggregate(cbind(attempted = rep(1L, nrow(d)),
        numerical_failures = !d$fit_usable, usable = d$eligible, unusable = !d$eligible),
      by = d[c("draw", "regime", "model")], FUN = sum)
  }
  failures <- rbind(failure_counts(ssb), failure_counts(transform(ssb, regime = "overall")))
  holdout_pairs <- .oe_pilot_holdout(fits, holdout)
  list(stock = stock, pairs = pairs, repeat_pairs = repeat_pairs,
    comparisons = comparisons, direct_comparisons = .oe_pilot_direct_comparisons(stock),
    sample_size = .oe_pilot_sample_size(pairs, comparisons),
    failures = failures, holdout = holdout_pairs,
    sd_validation = .oe_pilot_sd_validation(fits, holdout, sd_validation),
    holdout_comparisons = .oe_pilot_comparisons(holdout_pairs,
      c("draw", "model", "component"), .oe_pilot_holdout_metrics))
}
