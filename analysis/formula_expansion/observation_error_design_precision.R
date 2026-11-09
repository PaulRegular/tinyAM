# Year-held-out RV design precision. Source after observation_error_diagnostics.R.
oe_pooled_rv_sd <- function(train, ages) {
  stopifnot(all(c("year", "age", "obs", "se") %in% names(train)))
  reported <- subset(as.data.frame(train), is.finite(obs) & obs > 0 & is.finite(se) & se > 0)
  reported$variance <- oe_diag_reported_log_sd(reported$obs, reported$se)^2
  reported <- subset(reported, is.finite(variance) & variance > 0)
  if (!nrow(reported)) stop("No positive design-based RV variances in training years.")

  # Give years with valid SEs equal weight at each age.
  yearly <- aggregate(variance ~ year + age, reported, mean)
  pooled <- aggregate(variance ~ age, yearly, mean)
  variance <- pooled$variance[match(ages, pooled$age)]
  fallback <- is.na(variance)
  variance[fallback] <- median(pooled$variance)
  structure(sqrt(variance), fallback = fallback)
}

oe_pooled_score_year <- function(train, test, held_year) {
  train <- subset(train, year != held_year & is.finite(centered))
  test <- subset(test, year == held_year & is.finite(centered))
  if (!nrow(train) || !nrow(test)) stop("An RV year split has no positive observations.")

  raw_sd <- oe_pooled_rv_sd(train, c(train$age, test$age))
  constant <- sqrt(mean(train$centered^2))
  multiplier <- sqrt(mean((train$centered / head(raw_sd, nrow(train)))^2))
  pooled <- multiplier * tail(raw_sd, nrow(test))
  stopifnot(is.finite(constant), constant > 0, all(is.finite(pooled) & pooled > 0))
  cutoff <- qnorm(0.975)
  data.frame(replicate = test$replicate,
    delta = dnorm(test$centered, sd = constant, log = TRUE) -
      dnorm(test$centered, sd = pooled, log = TRUE),
    coverage_constant = abs(test$centered) <= cutoff * constant,
    coverage_pooled = abs(test$centered) <= cutoff * pooled,
    fallback = tail(attr(raw_sd, "fallback"), nrow(test)))
}

oe_pooled_year_scores <- function(rows, calibration, validation) {
  groups <- c("om", "population_id", "effort")
  keys <- c(groups, "year", "age", "replicate")
  required <- c(keys, "component", "bin", "obs", "truth", "se")
  stopifnot(all(required %in% names(rows)), length(calibration) >= 2L,
    length(validation) >= 2L, !length(intersect(calibration, validation)))
  rv <- subset(as.data.frame(rows), component == "RV" & bin == "modeled" &
    replicate %in% c(calibration, validation))
  stopifnot(nrow(rv) > 0L, !anyNA(rv[keys]), !anyDuplicated(rv[keys]))
  positive <- is.finite(rv$obs) & rv$obs > 0 & is.finite(rv$truth) & rv$truth > 0
  rv$log_error <- NA_real_
  rv$log_error[positive] <- log(rv$obs[positive]) - log(rv$truth[positive])
  rv$cell <- paste(rv$year, rv$age)
  min_n <- if (length(calibration) < 10L) 2L else 10L
  scores <- do.call(rbind, lapply(split(rv, rv[groups], drop = TRUE), function(d) {
    # Same-year independent calibration estimates the log center only.
    # The SD shape and its multiplier are estimated using other years.
    z <- oe_diag_center(subset(d, replicate %in% calibration),
      subset(d, replicate %in% validation), min_n)
    out <- do.call(rbind, lapply(sort(unique(d$year)), function(year) {
      oe_pooled_score_year(z$train, z$test, year)
    }))
    out$effort <- d$effort[1L]
    out
  }))
  result <- do.call(rbind, lapply(split(scores, scores$effort, drop = TRUE), function(d) {
    draw <- aggregate(cbind(delta, coverage_constant, coverage_pooled) ~ replicate, d, mean)
    n <- nrow(draw)
    delta <- mean(draw$delta)
    half_width <- if (n > 1L) qt(.975, n - 1L) * sd(draw$delta) / sqrt(n) else NA_real_
    data.frame(effort = d$effort[1L], delta = delta, lower = delta - half_width,
      upper = delta + half_width, n_draws = n, n_scored = nrow(d),
      n_fallback = sum(d$fallback), better_draws = sum(draw$delta < 0),
      coverage_constant = mean(draw$coverage_constant), coverage_pooled = mean(draw$coverage_pooled))
  }))
  rownames(result) <- NULL
  result
}
