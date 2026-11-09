# The independent unit is a mortality/spatial history. Sampling draws within
# a fixed state estimate monitoring covariance and its finite-draw correction.
.oe_z_columns <- c("scenario", "spatial", "mortality", "phi_year", "component")

.oe_z_state <- function(rows) {
  ages <- sort(unique(rows$age))
  draws <- sort(unique(rows$replicate))
  key <- rows[1L, c(.oe_z_columns, "history_id", "year")]
  rownames(key) <- NULL
  stopifnot(identical(as.integer(ages), 1:14), !anyDuplicated(rows[c("age", "replicate")]),
    all(rows$obs >= 0 | is.na(rows$obs)))
  target <- rows$truth[match(ages, rows$age)]
  stopifnot(all(is.finite(target)), all(target[match(rows$age, ages)] == rows$truth))

  obs <- matrix(NA_real_, nrow = length(draws), ncol = length(ages))
  obs[cbind(match(rows$replicate, draws), match(rows$age, ages))] <- rows$obs
  arithmetic <- obs[rowSums(is.finite(obs)) == length(ages), , drop = FALSE]
  n_complete <- nrow(arithmetic)

  log_obs <- obs
  log_obs[!is.finite(obs) | obs <= 0 | matrix(target <= 0, nrow(obs), ncol(obs), byrow = TRUE)] <- NA_real_
  log_obs <- log(log_obs)
  joint <- crossprod(is.finite(log_obs) * 1)
  log_cor <- suppressWarnings(stats::cor(log_obs, use = "pairwise.complete.obs"))
  log_cor[joint < 10L] <- NA_real_

  list(key = key, ages = ages, target = target, log_cor = log_cor, joint = joint,
    n_complete = n_complete, n_draws = length(draws),
    center = if (n_complete >= 2L) colMeans(arithmetic),
    within = if (n_complete >= 2L) stats::cov(arithmetic),
    n_obs = sum(is.finite(obs)), n_zero = sum(obs == 0, na.rm = TRUE),
    n_missing = sum(!is.finite(obs)))
}

.oe_z_average <- function(x) Reduce("+", x) / length(x)

.oe_z_covariance <- function(states, index) {
  states <- states[index]
  states <- Filter(function(x) x$n_complete >= 2L, states)
  if (length(states) < 2L) return(NULL)
  means <- do.call(rbind, lapply(states, "[[", "center"))
  targets <- do.call(rbind, lapply(states, "[[", "target"))
  within <- .oe_z_average(lapply(states, "[[", "within"))
  mean_noise <- .oe_z_average(lapply(states, function(x) x$within / x$n_complete))
  between <- stats::cov(means) - mean_noise
  list(within = within, between = between, target = stats::cov(targets),
    histories = length(states))
}

.oe_z_decomposition <- function(states, index) {
  years <- vapply(states[index], function(x) x$key$year, numeric(1))
  by_year <- lapply(split(index, years), function(i) .oe_z_covariance(states, i))
  by_year <- Filter(Negate(is.null), by_year)
  if (!length(by_year)) return(NULL)

  within <- .oe_z_average(lapply(by_year, "[[", "within"))
  between <- .oe_z_average(lapply(by_year, "[[", "between"))
  target <- .oe_z_average(lapply(by_year, "[[", "target"))
  total <- within + between
  ages <- states[[index[1L]]]$ages
  p <- which(upper.tri(total), arr.ind = TRUE)
  sd_products <- function(cov) outer(sqrt(pmax(diag(cov), 0)), sqrt(pmax(diag(cov), 0)))[p]
  divide <- function(num, den) ifelse(den > 0, num / den, NA_real_)
  den <- sd_products(total)
  data.frame(states[[index[1L]]]$key[.oe_z_columns], age_a = ages[p[, 1L]],
    age_b = ages[p[, 2L]], age_gap = ages[p[, 2L]] - ages[p[, 1L]],
    within_cov = within[p], between_cov = between[p], total_cov = total[p],
    target_cov = target[p], within_corr = divide(within[p], sd_products(within)),
    within_contribution = divide(within[p], den),
    between_contribution = divide(between[p], den), total_corr = divide(total[p], den),
    target_corr = divide(target[p], sd_products(target)),
    n_years = length(by_year), min_histories = min(vapply(by_year, "[[", integer(1), "histories")))
}

.oe_z_log_pairs <- function(state) {
  p <- which(upper.tri(state$log_cor) & is.finite(state$log_cor), arr.ind = TRUE)
  if (!nrow(p)) return(NULL)
  age_a <- state$ages[p[, 1L]]
  age_b <- state$ages[p[, 2L]]
  data.frame(state$key, age_a = age_a, age_b = age_b, age_gap = age_b - age_a,
    estimate = state$log_cor[p], n_draw_pairs = state$joint[p])
}

.oe_z_interval <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  half <- if (n >= 2L) stats::qt(0.975, n - 1) * stats::sd(x) / sqrt(n) else NA_real_
  data.frame(estimate = if (n) mean(x) else NA_real_, lower = mean(x) - half,
    upper = mean(x) + half, n_histories = n)
}

.oe_z_summarize_log <- function(states) {
  age_pairs <- do.call(rbind, lapply(states, .oe_z_log_pairs))
  if (is.null(age_pairs)) return(list(lags = data.frame(), summary = data.frame(),
    contrasts = data.frame()))
  snapshot_keys <- c(.oe_z_columns, "history_id", "year", "age_gap")
  snapshot <- do.call(rbind, lapply(split(age_pairs, age_pairs[snapshot_keys], drop = TRUE), function(x) {
    data.frame(x[1L, snapshot_keys], estimate = mean(x$estimate), n_pairs = nrow(x),
      n_draw_pairs = sum(x$n_draw_pairs))
  }))
  keys <- c(.oe_z_columns, "history_id", "age_gap")
  lags <- do.call(rbind, lapply(split(snapshot, snapshot[keys], drop = TRUE), function(x) {
    data.frame(x[1L, keys], estimate = mean(x$estimate), n_years = nrow(x),
      n_pairs = sum(x$n_pairs), n_draw_pairs = sum(x$n_draw_pairs))
  }))
  keys <- c(.oe_z_columns, "age_gap")
  summary <- do.call(rbind, lapply(split(lags, lags[keys], drop = TRUE), function(x) {
    data.frame(x[1L, keys], .oe_z_interval(x$estimate), n_pairs = sum(x$n_pairs),
      n_draw_pairs = sum(x$n_draw_pairs))
  }))

  pair_keys <- c("spatial", "component", "history_id", "year", "age_a", "age_b")
  reference <- age_pairs[age_pairs$mortality == "iid" & age_pairs$phi_year == 0,
    c(pair_keys, "estimate", "scenario")]
  names(reference)[match(c("estimate", "scenario"), names(reference))] <-
    c("reference_estimate", "reference_scenario")
  stopifnot(!anyDuplicated(reference[pair_keys]))
  pairs <- merge(age_pairs, reference, by = pair_keys, sort = FALSE)
  pairs <- pairs[pairs$scenario != pairs$reference_scenario, , drop = FALSE]
  pairs$difference <- pairs$estimate - pairs$reference_estimate
  contrast_keys <- c(keys, "reference_scenario")
  snapshot_keys <- c(contrast_keys, "history_id", "year")
  paired_year <- do.call(rbind, lapply(split(pairs, pairs[snapshot_keys], drop = TRUE), function(x) {
    data.frame(x[1L, snapshot_keys], difference = mean(x$difference), n_pairs = nrow(x))
  }))
  contrasts <- if (NROW(paired_year)) {
    do.call(rbind, lapply(split(paired_year, paired_year[contrast_keys], drop = TRUE), function(x) {
      effects <- vapply(split(x$difference, x$history_id), mean, numeric(1))
      data.frame(x[1L, contrast_keys], .oe_z_interval(effects), n_pairs = sum(x$n_pairs))
    }))
  } else data.frame()
  list(lags = lags, summary = summary, contrasts = contrasts)
}

oe_z_summarize <- function(rows) {
  stopifnot(is.data.frame(rows), nrow(rows) > 0L)
  rows <- as.data.frame(rows)
  keys <- c(.oe_z_columns, "history_id", "year")
  states <- lapply(split(rows, rows[keys], drop = TRUE), .oe_z_state)
  state_keys <- do.call(rbind, lapply(states, "[[", "key"))
  groups <- split(seq_along(states), state_keys[.oe_z_columns], drop = TRUE)
  decomposition <- do.call(rbind, lapply(groups, function(i) .oe_z_decomposition(states, i)))
  support <- do.call(rbind, lapply(groups, function(i) {
    columns <- c("n_obs", "n_zero", "n_missing", "n_complete", "n_draws")
    values <- vapply(columns, function(name) {
      sum(vapply(states[i], "[[", numeric(1), name))
    }, numeric(1))
    data.frame(states[[i[1L]]]$key[.oe_z_columns], as.list(values), n_states = length(i))
  }))
  log <- .oe_z_summarize_log(states)
  list(decomposition = decomposition, log_lags = log$lags, log_summary = log$summary,
    contrasts = log$contrasts, support = support)
}

oe_z_residual_scale <- function(rows) {
  rows <- as.data.frame(rows)
  rows <- rows[rows$replicate > 2L, ]
  keys <- c("component", "age", "history_id", "year")
  grouped <- split(rows$obs, rows[keys], drop = TRUE)
  variances <- lapply(grouped, function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) >= 5L) stats::var(log(x)) else NA_real_
  })
  by_state <- lapply(split(rows, rows[keys], drop = TRUE), function(x) x[1L, keys])
  state <- do.call(rbind, by_state)
  state$variance <- unlist(variances, use.names = FALSE)
  do.call(rbind, lapply(split(state, state[c("component", "age")], drop = TRUE), function(x) {
    positive <- x$variance[is.finite(x$variance) & x$variance > 0]
    data.frame(x[1L, c("component", "age")], sd = if (length(positive)) sqrt(mean(positive)) else NA_real_,
      n_states = length(positive))
  }))
}
