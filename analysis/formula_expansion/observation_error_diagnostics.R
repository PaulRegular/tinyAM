# Summarize repeated observations against fixed truth, outside the assessment.
oe_diag_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

oe_diag_sd <- function(x) stats::sd(x[is.finite(x)])

oe_diag_interval <- function(x, bootstrap) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || bootstrap < 1L) return(c(NA_real_, NA_real_))
  samples <- matrix(sample(x, length(x) * bootstrap, replace = TRUE), nrow = length(x))
  stats::quantile(colMeans(samples), c(0.025, 0.975), names = FALSE)
}

oe_diag_cells <- function(x, summarize) {
  do.call(rbind, lapply(split(x, x$cell), function(d) {
    cbind(d[1L, c("year", "age"), drop = FALSE], summarize(d))
  }))
}

oe_diag_age_year <- function(d) {
  observed <- is.finite(d$obs)
  true_positive <- is.finite(d$truth) & d$truth > 0
  n_observed <- sum(observed)
  n_zero <- sum(observed & d$obs == 0)
  n_missing <- sum(!observed)
  n_log <- sum(is.finite(d$log_error))
  zero_rate <- if (n_observed > 0L) n_zero / n_observed else NA_real_
  missing_rate <- n_missing / nrow(d)
  mean_obs <- oe_diag_mean(d$obs)
  sd_obs <- oe_diag_sd(d$obs)
  sd_log_error <- oe_diag_sd(d$log_error)
  out <- data.frame(
    n_draws = nrow(d), n_observed = n_observed, n_zero = n_zero,
    n_missing = n_missing, n_positive = sum(observed & d$obs > 0),
    zero_rate = zero_rate, missing_rate = missing_rate,
    zero_rate_mcse = sqrt(zero_rate * (1 - zero_rate) / n_observed),
    missing_rate_mcse = sqrt(missing_rate * (1 - missing_rate) / nrow(d)),
    n_log = n_log, mean_obs = mean_obs, sd_obs = sd_obs,
    cv_obs = if (is.finite(mean_obs) && mean_obs > 0) sd_obs / mean_obs else NA_real_,
    mcse_obs = sd_obs / sqrt(n_observed),
    mean_truth = oe_diag_mean(d$truth),
    bias = oe_diag_mean(d$obs - d$truth),
    relative_bias = oe_diag_mean((d$obs[true_positive] - d$truth[true_positive]) /
      d$truth[true_positive]),
    mean_log_error = oe_diag_mean(d$log_error), sd_log_error = sd_log_error,
    mcse_log_error = sd_log_error / sqrt(n_log),
    mean_selectivity = oe_diag_mean(d$selectivity),
    mean_N_frame = oe_diag_mean(d$N_frame),
    mean_expected_encounters = oe_diag_mean(d$expected_encounters)
  )
  for (column in intersect(c("n_caught_at_age", "n_aged_at_age", "n_primary",
    "n_equivalent", "n_caught", "n_measured", "n_aged"), names(d))) {
    out[[paste0("mean_", column)]] <- oe_diag_mean(d[[column]])
    out[[paste0(column, "_draws")]] <- sum(is.finite(d[[column]]))
  }
  out
}

oe_diag_reported_log_sd <- function(obs, se) {
  # Stable evaluation of sqrt(log1p((se / obs)^2)).
  log_cv <- log(se) - log(obs)
  sqrt(pmax(2 * log_cv, 0) + log1p(exp(-abs(2 * log_cv))))
}

oe_diag_rv_precision <- function(d) {
  arithmetic <- is.finite(d$obs) & is.finite(d$truth) & d$truth > 0 &
    is.finite(d$se) & d$se >= 0
  logged <- is.finite(d$log_error) & is.finite(d$se) & d$se > 0
  supplied_sd <- oe_diag_reported_log_sd(d$obs[logged], d$se[logged])
  z <- stats::qnorm(0.975)
  data.frame(
    n_arithmetic = sum(arithmetic),
    n_log = sum(logged),
    empirical_sd = oe_diag_sd(d$obs[arithmetic]),
    rms_se = sqrt(oe_diag_mean(d$se[arithmetic]^2)),
    empirical_log_sd = oe_diag_sd(d$log_error[logged]),
    rms_log_sd = sqrt(oe_diag_mean(supplied_sd^2)),
    arithmetic_coverage = oe_diag_mean(abs(d$obs[arithmetic] - d$truth[arithmetic]) <=
      z * d$se[arithmetic]),
    log_coverage = oe_diag_mean(abs(d$log_error[logged]) <= z * supplied_sd)
  )
}

oe_diag_rv_precision_scores <- function(train, test, bootstrap) {
  eligible <- function(d) {
    d <- d[is.finite(d$centered) & is.finite(d$obs) & d$obs > 0 &
      is.finite(d$se) & d$se > 0, ]
    d$reported_log_sd <- oe_diag_reported_log_sd(d$obs, d$se)
    d[is.finite(d$reported_log_sd) & d$reported_log_sd > 0, ]
  }
  train <- eligible(train)
  test <- eligible(test)
  constant_sd <- sqrt(mean(train$centered^2))
  multiplier <- sqrt(mean((train$centered / train$reported_log_sd)^2))
  sds <- list(
    constant = rep(constant_sd, nrow(test)),
    supplied_fixed = test$reported_log_sd,
    supplied_scaled = multiplier * test$reported_log_sd
  )
  losses <- lapply(sds, function(sigma) {
    if (any(!is.finite(sigma) | sigma <= 0)) return(rep(NA_real_, length(sigma)))
    loss <- -stats::dnorm(test$centered, sd = sigma, log = TRUE)
    if (any(!is.finite(loss))) loss[] <- NA_real_
    loss
  })
  by_draw <- function(x) vapply(split(x, test$replicate), oe_diag_mean, numeric(1))
  q <- stats::qnorm(0.975)
  do.call(rbind, lapply(names(sds), function(model) {
    sigma <- sds[[model]]
    loss <- losses[[model]]
    score <- by_draw(loss)
    delta <- by_draw(loss - losses$constant)
    covered <- by_draw(ifelse(is.finite(loss), abs(test$centered) <= q * sigma, NA))
    score_ci <- oe_diag_interval(score, bootstrap)
    delta_ci <- oe_diag_interval(delta, bootstrap)
    covered_ci <- oe_diag_interval(covered, bootstrap)
    status <- if (!nrow(test)) "no_eligible_heldout" else if (
      any(!is.finite(sigma) | sigma <= 0)) "invalid_sd" else if (
      any(!is.finite(loss))) "nonfinite_scores" else "ok"
    data.frame(
      model = model, status = status,
      constant_sd = if (model == "constant") constant_sd else NA_real_,
      multiplier = switch(model, supplied_fixed = 1, supplied_scaled = multiplier, NA_real_),
      n_calibration = nrow(train), n_validation = nrow(test),
      n_validation_draws = length(unique(test$replicate)), n_scored = sum(is.finite(loss)),
      nll = oe_diag_mean(score), nll_lower = score_ci[1L], nll_upper = score_ci[2L],
      delta = oe_diag_mean(delta), delta_lower = delta_ci[1L], delta_upper = delta_ci[2L],
      coverage = oe_diag_mean(covered),
      coverage_lower = covered_ci[1L], coverage_upper = covered_ci[2L]
    )
  }))
}

oe_diag_sd_result <- function(sigma, status = "ok", warnings = "") {
  if (!length(sigma) && status == "ok") status <- "no_eligible_heldout"
  if (status == "ok" && any(!is.finite(sigma) | sigma <= 0)) status <- "invalid_sd"
  if (status != "ok") sigma[] <- NA_real_
  list(sigma = sigma, status = status, warnings = warnings)
}

oe_diag_sd_grouped <- function(error, groups, predicted) {
  sigma <- vapply(split(error, groups), function(x) sqrt(mean(x^2)), numeric(1))
  estimate <- unname(sigma[match(predicted, names(sigma))])
  oe_diag_sd_result(estimate, if (anyNA(estimate)) "unsupported_age" else "ok")
}

oe_diag_age_pair <- function(age) {
  # Pool 1–2 through 11–12, retain 13 and 14+, and pair older native ages.
  ifelse(age %in% 13:14, paste0("single_", age), paste0("pair_", (age - 1L) %/% 2L))
}

oe_diag_spline_design <- function(age) {
  # The null space supplies the fixed linear component. The penalized basis
  # is transformed to independent standard-scale random coefficients.
  smooth <- mgcv::smoothCon(mgcv::s(age, k = 5, bs = "cr"),
    data.frame(age = age), absorb.cons = TRUE)[[1L]]
  eigen_s <- eigen(smooth$S[[1L]], symmetric = TRUE)
  penalized <- eigen_s$values > max(eigen_s$values) * 1e-10
  list(x = cbind(1, smooth$X %*% eigen_s$vectors[, !penalized, drop = FALSE]),
    z = sweep(smooth$X %*% eigen_s$vectors[, penalized, drop = FALSE],
      2L, sqrt(eigen_s$values[penalized]), "/"))
}

oe_diag_sd_latent <- function(train, test, smooth = FALSE) {
  ages <- sort(unique(train$age))
  blank <- rep(NA_real_, nrow(test))
  minimum_ages <- if (smooth) 5L else 2L
  if (length(ages) < minimum_ages) {
    return(oe_diag_sd_result(blank, "insufficient_age_support"))
  }
  if (any(!test$age %in% ages)) return(oe_diag_sd_result(blank, "unsupported_age"))
  # This is the exact conditional Normal likelihood, up to its constant,
  # compressed to age-level sufficient statistics for the AD calculation.
  by_age <- split(train$centered, factor(train$age, levels = ages))
  n <- lengths(by_age)
  ss <- vapply(by_age, function(y) sum(y^2), numeric(1))
  scale <- sqrt(sum(ss) / sum(n))
  if (!is.finite(scale) || scale <= 0) return(oe_diag_sd_result(blank, "zero_training_spread"))
  warning_text <- character()
  attempt <- tryCatch(withCallingHandlers({
    design <- if (smooth) oe_diag_spline_design(ages) else
      list(x = matrix(1, length(ages)), z = diag(length(ages)))
    objective <- function(p) {
      log_sd <- drop(design$x %*% p$beta + design$z %*% p$u)
      sum(n * log_sd + 0.5 * ss * exp(-2 * log_sd)) -
        sum(RTMB::dnorm(p$u, 0, exp(p$log_tau), log = TRUE))
    }
    obj <- RTMB::MakeADFun(objective,
      list(beta = c(log(scale), rep(0, ncol(design$x) - 1L)),
        u = numeric(ncol(design$z)), log_tau = log(0.2)), random = "u", silent = TRUE)
    lower <- ifelse(names(obj$par) == "log_tau", -12, -Inf)
    upper <- ifelse(names(obj$par) == "log_tau", 8, Inf)
    fit <- stats::nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper)
    grad <- obj$gr(fit$par)
    grad[fit$par == lower & grad > 0 | fit$par == upper & grad < 0] <- 0
    if (fit$convergence != 0L || !is.finite(fit$objective) ||
      any(!is.finite(grad)) || max(abs(grad)) > 0.01) {
      oe_diag_sd_result(blank, "optimization_failure", fit$message)
    } else {
      obj$fn(fit$par)
      p <- obj$env$parList()
      sigma <- exp(drop(design$x %*% p$beta + design$z %*% p$u))
      out <- oe_diag_sd_result(sigma[match(test$age, ages)])
      out$log_tau <- unname(p$log_tau)
      out$boundary <- isTRUE(fit$par["log_tau"] %in% c(-12, 8))
      out
    }
  }, warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) oe_diag_sd_result(blank, "fit_error", conditionMessage(e)))
  attempt$warnings <- paste(unique(c(attempt$warnings, warning_text)), collapse = " | ")
  attempt
}

oe_diag_cadigan <- function(train, test) {
  # Model the variance of the index mean in fish per tow, then transform
  # back to the supplied index units (thousand fish).
  valid <- function(d, m) is.finite(m) & m > 0 & is.finite(d$n_equivalent) &
    d$n_equivalent > 0 & is.finite(d$n_primary) & d$n_primary > 0
  variance <- function(d, mean, theta) {
    mean_tow <- 1000 * mean / d$n_equivalent
    (d$n_equivalent / 1000)^2 * (mean_tow + theta * mean_tow^2) / d$n_primary
  }
  is_standard <- train$effort %in% c("standard", "1")
  ref <- train[is_standard & is.finite(train$obs), , drop = FALSE]
  key <- function(d) paste(d$year, d$age, sep = "/")
  means <- vapply(split(ref$obs, key(ref)), mean, numeric(1))
  methods <- c(oracle = "truth", calibration = "cadigan_mean")
  ref$cadigan_mean <- unname(means[key(ref)])
  test$cadigan_mean <- unname(means[key(test)])
  train_mean <- unname(means[key(train)])
  train_ok <- valid(train, train$truth) & valid(train, train_mean)
  test_ok <- valid(test, test$truth) & valid(test, test$cadigan_mean)
  common <- valid(ref, ref$truth) & valid(ref, ref$cadigan_mean) & is.finite(ref$arith_centered)
  ref <- ref[common, , drop = FALSE]
  fits <- lapply(methods, function(field) {
    fail <- function(status) list(
      original = oe_diag_sd_result(rep(NA_real_, nrow(test)), status),
      log = oe_diag_sd_result(rep(NA_real_, nrow(test)), status), theta = NA_real_)
    if (nrow(ref) < 2L) return(fail("insufficient_standard_calibration"))
    objective <- function(theta) {
      v <- variance(ref, ref[[field]], theta)
      sum(log(v) + ref$arith_centered^2 / v) / 2
    }
    fit <- tryCatch(stats::nlminb(0.01, objective, lower = 0), error = function(e) NULL)
    if (is.null(fit) || fit$convergence != 0L || !is.finite(fit$objective)) {
      return(fail("optimization_failure"))
    }
    m <- test[[field]]
    ok <- valid(test, m)
    v <- rep(NA_real_, nrow(test))
    v[ok] <- variance(test[ok, , drop = FALSE], m[ok], fit$par)
    list(original = list(sigma = sqrt(v), status = "ok", warnings = ""),
      log = list(sigma = oe_diag_reported_log_sd(m, sqrt(v)), status = "ok", warnings = ""),
      theta = unname(fit$par))
  })
  fits$n_standard_calibration <- nrow(ref)
  fits$train_ok <- train_ok
  fits$test_ok <- test_ok
  fits$calibration_mean <- test$cadigan_mean
  fits
}

oe_diag_score_variance <- function(fits, test, error, n_train, bootstrap) {
  by_draw <- function(x) vapply(split(x, test$replicate), oe_diag_mean, numeric(1))
  losses <- lapply(fits, function(fit) -stats::dnorm(error, sd = fit$sigma, log = TRUE))
  draws <- lapply(names(fits), function(model) {
    nll <- by_draw(losses[[model]])
    data.frame(model = rep(model, length(nll)), replicate = as.integer(names(nll)),
      nll = unname(nll), delta = unname(nll - by_draw(losses$constant)))
  })
  out <- lapply(seq_along(draws), function(i) {
    model <- names(fits)[i]
    fit <- fits[[i]]
    d <- draws[[i]]
    ci <- oe_diag_interval(d$delta, bootstrap)
    data.frame(model = model, status = fit$status, warnings = fit$warnings,
      n_calibration = n_train, n_validation = nrow(test),
      n_validation_draws = length(unique(test$replicate)),
      n_scored = sum(is.finite(losses[[i]])), nll = oe_diag_mean(d$nll),
      delta = oe_diag_mean(d$delta), delta_lower = ci[1L], delta_upper = ci[2L],
      coverage = oe_diag_mean(ifelse(is.finite(losses[[i]]),
        abs(error) <= stats::qnorm(0.975) * fit$sigma, NA)),
      log_tau = if (is.null(fit$log_tau)) NA_real_ else fit$log_tau,
      boundary = if (is.null(fit$boundary)) NA else fit$boundary)
  })
  list(scores = do.call(rbind, out), draws = do.call(rbind, draws))
}

oe_diag_variance <- function(train, test, bootstrap, extra = list()) {
  keep <- is.finite(test$centered)
  train <- train[is.finite(train$centered), , drop = FALSE]
  test <- test[keep, , drop = FALSE]
  constant <- oe_diag_sd_result(rep(sqrt(mean(train$centered^2)), nrow(test)))
  fits <- list(constant = constant,
    paired_age = oe_diag_sd_grouped(train$centered, oe_diag_age_pair(train$age),
      oe_diag_age_pair(test$age)),
    iid_age = oe_diag_sd_latent(train, test),
    smooth_age = oe_diag_sd_latent(train, test, smooth = TRUE),
    free_age = oe_diag_sd_grouped(train$centered, train$age, test$age))
  fits <- c(fits, lapply(extra, function(fit) {
    fit$sigma <- fit$sigma[keep]
    fit
  }))
  oe_diag_score_variance(fits, test, test$centered, nrow(train), bootstrap)
}

oe_diag_cadigan_subset <- function(cadigan, scale, i) {
  methods <- c("oracle", "calibration")
  fits <- stats::setNames(lapply(methods, function(method) {
    fit <- cadigan[[method]][[scale]]
    oe_diag_sd_result(fit$sigma[i], fit$status, fit$warnings)
  }), paste0("cadigan_", methods))
  fits
}

oe_diag_cadigan_table <- function(cadigan, train, test, training, validation, bootstrap) {
  theta <- c(constant = NA_real_, cadigan_oracle = cadigan$oracle$theta,
    cadigan_calibration = cadigan$calibration$theta)
  n_standard <- c(constant = NA_integer_,
    cadigan_oracle = cadigan$n_standard_calibration,
    cadigan_calibration = cadigan$n_standard_calibration)
  do.call(rbind, lapply(c("original", "log"), function(scale) {
    centered <- if (scale == "original") "arith_centered" else "centered"
    j <- training & cadigan$train_ok & is.finite(train[[centered]])
    i <- validation & cadigan$test_ok & is.finite(test[[centered]])
    y <- test[[centered]][i]
    fits <- c(list(constant = oe_diag_sd_result(rep(sqrt(mean(train[[centered]][j]^2)), length(y)))),
      oe_diag_cadigan_subset(cadigan, scale, i))
    tab <- oe_diag_score_variance(fits, test[i, , drop = FALSE], y, sum(j), bootstrap)$scores
    tab$theta <- unname(theta[tab$model])
    tab$n_standard_calibration <- unname(n_standard[tab$model])
    tab$scale <- scale
    tab
  }))
}

oe_diag_oracle <- function(d) {
  valid <- is.finite(d$oracle_obs) & is.finite(d$truth) &
    is.finite(d$oracle_se) & d$oracle_se >= 0
  error <- d$oracle_obs[valid] - d$truth[valid]
  data.frame(n_oracle = sum(valid), oracle_mean = oe_diag_mean(d$oracle_obs[valid]),
    oracle_bias = oe_diag_mean(error), oracle_empirical_sd = oe_diag_sd(error),
    oracle_design_sd = sqrt(oe_diag_mean(d$oracle_se[valid]^2)),
    oracle_coverage = oe_diag_mean(abs(error) <= stats::qnorm(0.975) * d$oracle_se[valid]))
}

oe_diag_variance_selection <- function(scores, draws, bootstrap) {
  # Equal weight per monitoring draw, then across efforts and realizations
  # within OM, and finally across the three OMs. Retain only common support.
  complexity <- c(constant = 0, cadigan_calibration = 1, paired_age = 2,
    iid_age = 3, smooth_age = 4, free_age = 5, cadigan_oracle = Inf)
  practical <- c("constant", "paired_age", "iid_age", "free_age", "cadigan_calibration")
  if (!nrow(draws)) {
    draws <- scores[FALSE, c("om", "population_id", "component", "bin", "effort", "model", "delta")]
    draws$replicate <- integer()
  }
  output <- lapply(split(scores, interaction(scores$component, scores$bin, drop = TRUE)),
    function(d) {
      s <- draws[draws$component == d$component[1L] & draws$bin == d$bin[1L], , drop = FALSE]
      models <- names(complexity)[names(complexity) %in% d$model]
      base <- s[s$model == "constant", c("om", "population_id", "effort", "replicate", "delta")]
      base <- base[is.finite(base$delta), , drop = FALSE]
      key <- function(x) do.call(paste, c(x[c("om", "population_id", "effort", "replicate")], sep = "/"))
      deltas <- vector("list", length(models))
      names(deltas) <- models
      rows <- lapply(models, function(model) {
        x <- s[s$model == model, , drop = FALSE]
        x <- x[match(key(base), key(x)), , drop = FALSE]
        complete <- nrow(base) > 0L && all(is.finite(x$delta)) &&
          all(d$status[d$model == model & d$n_validation > 0L] == "ok")
        by_om <- if (complete) stats::aggregate(delta ~ om + replicate, x, mean) else NULL
        paired <- if (complete) vapply(split(by_om$delta, by_om$replicate), mean, numeric(1)) else numeric()
        om_means <- if (complete) vapply(split(by_om$delta, by_om$om), mean, numeric(1)) else numeric()
        ci <- oe_diag_interval(paired, bootstrap)
        deltas[[model]] <<- paired
        improved <- sum(om_means < 0)
        data.frame(component = d$component[1L], bin = d$bin[1L], model = model,
          status = if (complete) "ok" else "incomplete_support",
          delta = oe_diag_mean(paired), delta_lower = ci[1L], delta_upper = ci[2L],
          mcse = if (length(paired) > 1L) stats::sd(paired) / sqrt(length(paired)) else NA_real_,
          n_om_improved = improved, n_om = length(om_means),
          practical = model %in% practical,
          eligible = isTRUE(ci[2L] < 0) && improved >= 2L &&
            length(om_means) == length(unique(d$om)),
          within_best_mcse = FALSE, selected = FALSE, advance = FALSE)
      })
      tab <- do.call(rbind, rows)
      candidates <- which(tab$practical & tab$eligible)
      selected <- match("constant", tab$model)
      if (length(candidates)) {
        best <- candidates[which.min(tab$delta[candidates])]
        within <- vapply(candidates, function(i) {
          difference <- deltas[[tab$model[i]]] - deltas[[tab$model[best]]]
          isTRUE(mean(difference) <= stats::sd(difference) / sqrt(length(difference)))
        }, logical(1))
        tab$within_best_mcse[candidates] <- within
        simple <- candidates[within]
        selected <- simple[which.min(complexity[tab$model[simple]])]
      }
      tab$selected[selected] <- TRUE
      tab$advance[selected] <- tab$bin[selected] == "modeled" && tab$model[selected] != "constant"
      tab
    })
  do.call(rbind, output)
}

oe_diag_pairs <- function(test, scale = c("log", "original")) {
  scale <- match.arg(scale)
  pairs <- lapply(split(test, test$year), function(d) {
    ages <- sort(unique(d$age))
    if (length(ages) < 2L) return(NULL)
    reps <- sort(unique(d$replicate))
    errors <- matrix(NA_real_, nrow = length(reps), ncol = length(ages))
    error <- if (scale == "log") d$log_error else d$obs - d$truth
    errors[cbind(match(d$replicate, reps), match(d$age, ages))] <- error
    pairs <- utils::combn(seq_along(ages), 2L)
    values <- vapply(seq_len(ncol(pairs)), function(i) {
      x <- errors[, pairs[1L, i]]
      y <- errors[, pairs[2L, i]]
      keep <- is.finite(x) & is.finite(y)
      n <- sum(keep)
      if (n < 2L) return(c(n, rep(NA_real_, 6L)))
      x <- x[keep]
      y <- y[keep]
      cross <- (x - mean(x)) * (y - mean(y))
      cov <- sum(cross) / (n - 1)
      sd_x <- stats::sd(x)
      sd_y <- stats::sd(y)
      corr <- if (sd_x > 0 && sd_y > 0) {
        max(-1, min(1, cov / (sd_x * sd_y)))
      } else NA_real_
      # Pointwise, approximate intervals across independent observation draws.
      cov_delta <- if (n > 2L) stats::qt(0.975, n - 1L) * stats::sd(cross) *
        sqrt(n) / (n - 1L) else NA_real_
      corr_ci <- if (n > 3L) {
        tanh(atanh(corr) + c(-1, 1) * stats::qnorm(0.975) / sqrt(n - 3))
      } else c(NA_real_, NA_real_)
      c(n, cov, cov - cov_delta, cov + cov_delta, corr, corr_ci)
    }, numeric(7))
    data.frame(
      scale = scale, year = d$year[1L], age_a = ages[pairs[1L, ]],
      age_b = ages[pairs[2L, ]], gap = ages[pairs[2L, ]] - ages[pairs[1L, ]],
      n_draws = length(reps), n = values[1L, ], covariance = values[2L, ],
      covariance_lower = values[3L, ], covariance_upper = values[4L, ],
      correlation = values[5L, ], correlation_lower = values[6L, ],
      correlation_upper = values[7L, ]
    )
  })
  out <- do.call(rbind, pairs)
  if (is.null(out)) data.frame() else out
}

oe_diag_lags <- function(pairs) {
  if (!nrow(pairs)) return(data.frame())
  do.call(rbind, lapply(split(pairs, interaction(pairs$scale, pairs$gap, drop = TRUE)), function(d) {
    keep <- is.finite(d$correlation)
    data.frame(scale = d$scale[1L], gap = d$gap[1L],
      n_age_year_pairs = sum(keep), n_draw_pairs = sum(d$n[keep]),
      mean_correlation = oe_diag_mean(d$correlation[keep]),
      mean_covariance = oe_diag_mean(d$covariance[keep]))
  }))
}

oe_diag_level_share <- function(draws, min_n = 10L) {
  draws <- subset(draws, om == "coupled" & bin == "modeled" & effort == "standard" &
    component %in% c("Gillnet", "RV"))
  snapshots <- split(draws, draws[c("component", "population_id", "year")], drop = TRUE)
  results <- lapply(snapshots, function(d) {
    ages <- sort(unique(d$age))
    reps <- sort(unique(d$replicate))
    key <- cbind(match(d$replicate, reps), match(d$age, ages))
    grid <- function(x) {
      y <- matrix(NA_real_, length(reps), length(ages))
      y[key] <- x
      y[!is.finite(y) | y < 0] <- NA_real_
      y
    }
    positive_log <- function(x) {
      x[x == 0 & !is.na(x)] <- NA_real_
      log(x)
    }
    truth <- grid(d$truth)
    obs <- stats::setNames(list(grid(d$obs)), d$component[1L])
    if (d$component[1L] == "RV" && any(is.finite(d$oracle_obs))) {
      obs[["RV capture only"]] <- grid(d$oracle_obs)
    }
    errors <- lapply(obs, function(x) {
      total <- rowSums(x)
      level <- positive_log(x) - positive_log(truth)
      share <- level - (positive_log(total) - positive_log(rowSums(truth)))
      list(level = level, share = share, total = total)
    })
    total_cv <- do.call(rbind, lapply(names(errors), function(name) {
      total <- errors[[name]]$total
      total <- total[is.finite(total)]
      data.frame(component = name, population_id = d$population_id[1L], year = d$year[1L],
        n = length(total), cv = if (length(total) > 1L && mean(total) > 0) {
          stats::sd(total) / mean(total)
        } else NA_real_)
    }))
    pairs <- lapply(which(diff(ages) == 1L), function(i) {
      pair <- i + 0:1
      values <- lapply(errors, function(x) cbind(x$level[, pair], x$share[, pair]))
      keep <- Reduce(`&`, lapply(values, function(x) rowSums(!is.finite(x)) == 0L))
      n <- sum(keep)
      values <- lapply(values, function(x) x[keep, , drop = FALSE])
      usable <- n >= min_n && all(vapply(values, function(x) {
        all(apply(x, 2L, stats::sd) > 0)
      }, logical(1)))
      # Compare levels and shares, including the RV capture benchmark, on the
      # same positive draw pairs. A missing age anywhere precludes a draw total.
      do.call(rbind, lapply(names(values), function(name) {
        x <- values[[name]]
        data.frame(component = name, population_id = d$population_id[1L], year = d$year[1L],
          age_a = ages[i], age_b = ages[i + 1L], n = n,
          correlation_level = if (usable) stats::cor(x[, 1L], x[, 2L]) else NA_real_,
          correlation_share = if (usable) stats::cor(x[, 3L], x[, 4L]) else NA_real_)
      }))
    })
    list(level_share = do.call(rbind, pairs), total_cv = total_cv)
  })
  tables <- lapply(c("level_share", "total_cv"), function(name) {
    tab <- do.call(rbind, lapply(results, `[[`, name))
    if (is.null(tab)) return(data.frame())
    rownames(tab) <- NULL
    tab
  })
  stats::setNames(tables, c("level_share", "total_cv"))
}

oe_diag_ar1_data <- function(d) {
  d <- d[is.finite(d$centered) & is.finite(d$cell_sd) & d$cell_sd > 0, ]
  d <- d[order(d$replicate, d$year, d$age), ]
  d$z <- d$centered / d$cell_sd
  block <- paste(d$replicate, d$year)
  following <- which(block[-1L] == head(block, -1L)) + 1L
  list(d = d, following = following, gap = d$age[following] - d$age[following - 1L])
}

oe_diag_ar1_loss <- function(x, rho) {
  d <- x$d
  loss <- -stats::dnorm(d$z, log = TRUE) + log(d$cell_sd)
  current <- x$following
  phi <- rho^x$gap
  loss[current] <- -stats::dnorm(d$z[current], phi * d$z[current - 1L],
    sqrt(1 - phi^2), log = TRUE) + log(d$cell_sd[current])
  loss
}

oe_diag_covariance <- function(train, test, bootstrap) {
  a <- oe_diag_ar1_data(train)
  b <- oe_diag_ar1_data(test)
  rho <- NA_real_
  sign_identified <- NA
  if (length(a$following)) {
    loss <- function(rho) sum(oe_diag_ar1_loss(a, rho))
    # Even gaps give identical likelihoods for opposite signs.
    sign_identified <- any(a$gap %% 2L == 1L)
    grid <- seq(if (sign_identified) -0.99 else 0, 0.99, length.out = 41L)
    values <- vapply(grid, loss, numeric(1))
    best <- which.min(values)
    opt <- stats::optimize(loss, range(grid[pmax(1L, pmin(length(grid), best + c(-1L, 1L)))]))
    rho <- if (opt$objective < values[best]) opt$minimum else grid[best]
  }
  stopifnot(isTRUE(all.equal(oe_diag_ar1_loss(b, 0),
    -stats::dnorm(b$d$z, log = TRUE) + log(b$d$cell_sd), tolerance = 1e-12)))
  independent <- oe_diag_ar1_loss(b, 0)
  correlated <- if (is.finite(rho)) oe_diag_ar1_loss(b, rho) else rep(NA_real_, length(independent))
  draws <- do.call(rbind, lapply(split(seq_len(nrow(b$d)), b$d$replicate), function(i) {
    data.frame(replicate = b$d$replicate[i[1L]], n_obs = length(i),
      nll_independent = mean(independent[i]), nll_ar1 = mean(correlated[i]),
      delta = mean(correlated[i] - independent[i]))
  }))
  if (is.null(draws)) draws <- data.frame()
  ci <- oe_diag_interval(draws$delta, bootstrap)
  list(
    covariance_scores = data.frame(
      status = if (is.finite(rho)) "ok" else "insufficient_calibration",
      rho = rho, rho_sign_identified = sign_identified,
      n_calibration_pairs = length(a$following),
      n_validation_draws = nrow(draws), n_obs = nrow(b$d),
      n_blocks = length(unique(paste(b$d$replicate, b$d$year))),
      n_validation_pairs = length(b$following),
      nll_independent = oe_diag_mean(draws$nll_independent),
      nll_ar1 = oe_diag_mean(draws$nll_ar1), delta = oe_diag_mean(draws$delta),
      delta_lower = ci[1L], delta_upper = ci[2L]
    ),
    covariance_score_draws = draws
  )
}

oe_diag_center <- function(train, test, min_n) {
  moments <- function(y) {
    n <- sum(is.finite(y))
    if (n < min_n) return(c(n = n, mean = NA_real_, sd = NA_real_))
    c(n = n, mean = oe_diag_mean(y), sd = oe_diag_sd(y))
  }
  center <- function(field, result) {
    cells <- split(train[[field]], train$cell)
    m <- matrix(vapply(cells, moments, numeric(3)), ncol = 3L, byrow = TRUE,
      dimnames = list(names(cells), c("n", "mean", "sd")))
    j <- match(train$cell, rownames(m))
    i <- match(test$cell, rownames(m))
    # With a fitted cell mean, compensate for the variance lost from the
    # centered calibration residuals. Held-out residuals remain unchanged.
    train[[result]] <<- (train[[field]] - m[j, "mean"]) *
      sqrt(ifelse(m[j, "n"] >= 2L, m[j, "n"] / (m[j, "n"] - 1), NA_real_))
    test[[result]] <<- test[[field]] - m[i, "mean"]
    list(train = unname(m[j, "sd"]), test = unname(m[i, "sd"]))
  }
  log <- center("log_error", "centered")
  center("obs", "arith_centered")
  train$cell_sd <- log$train
  test$cell_sd <- log$test
  list(train = train, test = test)
}

oe_diagnostics <- function(draws, calibration = 1:30, validation = 31:60, bootstrap = 200L) {
  draws <- as.data.frame(draws)
  defaults <- c(population_id = "seed17", om = "reference", effort = "standard")
  for (name in setdiff(names(defaults), names(draws))) draws[[name]] <- rep(defaults[[name]], nrow(draws))
  id <- c("om", "population_id", "component", "bin")
  required <- c(id, "effort", "replicate", "year", "age", "obs", "truth", "se")
  stopifnot(all(required %in% names(draws)), length(calibration) >= 2L, length(validation) > 0L,
    !length(intersect(calibration, validation)),
    !anyNA(draws[c(id, "effort", "replicate", "year", "age")]),
    all(draws$age == as.integer(draws$age)),
    !anyDuplicated(draws[c(id, "effort", "replicate", "year", "age")]))
  if (!nrow(draws)) return(list())

  auxiliary <- intersect(c("selectivity", "expected_encounters", "N_frame", "n_caught_at_age",
    "n_aged_at_age", "n_primary", "n_equivalent", "n_caught", "n_measured", "n_aged",
    "oracle_obs", "oracle_se"), names(draws))
  groups <- split(draws[c(required, auxiliary)],
    interaction(draws[id], drop = TRUE, lex.order = TRUE))
  seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (is.null(seed)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", seed, envir = .GlobalEnv))
  set.seed(41051L)
  results <- lapply(groups, function(d) {
    keep <- is.finite(d$obs) & d$obs > 0 & is.finite(d$truth) & d$truth > 0
    d$log_error <- NA_real_
    d$log_error[keep] <- log(d$obs[keep]) - log(d$truth[keep])
    d$cell <- paste(d$year, d$age, d$effort)
    fixed <- vapply(split(d$truth, paste(d$year, d$age)), function(x) {
      all(is.finite(x)) && diff(range(x)) <= 1e-10 * max(1, abs(x))
    }, logical(1))
    if (!all(fixed)) stop("Observation truth changed across draws or efforts for ",
      d$population_id[1L], " / ", d$om[1L], ".")
    min_n <- if (length(calibration) < 10L) 2L else 10L
    centered <- oe_diag_center(d[d$replicate %in% calibration, , drop = FALSE],
      d[d$replicate %in% validation, , drop = FALSE], min_n)
    train <- centered$train
    test <- centered$test
    rv <- d$component[1L] == "RV"
    cadigan <- if (rv && all(c("n_equivalent", "n_primary") %in% names(d))) oe_diag_cadigan(train, test)
    lapply(unique(d$effort), function(effort) {
      training <- train$effort == effort
      validation <- test$effort == effort
      a <- train[training, , drop = FALSE]
      b <- test[validation, , drop = FALSE]
      if (nrow(a)) a$cell <- paste(a$year, a$age)
      if (nrow(b)) b$cell <- paste(b$year, b$age)
      pairs <- rbind(oe_diag_pairs(b, "original"), oe_diag_pairs(b, "log"))
      has_oracle <- rv && all(c("oracle_obs", "oracle_se") %in% names(b))
      offsets <- NULL
      if (!is.null(cadigan)) {
        v <- validation & cadigan$test_ok
        offsets <- unique(data.frame(year = test$year[v], age = test$age[v],
          sd_supplied = cadigan$calibration$log$sigma[v],
          theta = rep(cadigan$calibration$theta, sum(v)),
          calibration_mean = cadigan$calibration_mean[v], n_primary = test$n_primary[v],
          n_equivalent = test$n_equivalent[v]))
        if (anyDuplicated(offsets[c("year", "age")])) {
          stop("Cadigan offsets require one planned set count and area per prediction cell.")
        }
        variance <- oe_diag_variance(train[training & cadigan$train_ok, , drop = FALSE],
          test[v, , drop = FALSE], bootstrap, oe_diag_cadigan_subset(cadigan, "log", v))
      } else {
        variance <- oe_diag_variance(a, b, bootstrap)
      }
      tables <- c(list(
        age_year = oe_diag_cells(b, oe_diag_age_year),
        rv_precision = if (rv) oe_diag_cells(b, oe_diag_rv_precision),
        rv_precision_scores = if (rv) oe_diag_rv_precision_scores(a, b, bootstrap),
        variance = variance$scores,
        variance_score_draws = variance$draws,
        cadigan = if (!is.null(cadigan)) oe_diag_cadigan_table(cadigan, train, test,
          training, validation, bootstrap),
        cadigan_offsets = offsets,
        oracle = if (has_oracle) oe_diag_cells(b, oe_diag_oracle),
        covariance_pairs = pairs,
        covariance_lags = oe_diag_lags(pairs)
      ), oe_diag_covariance(a, b, bootstrap))
      lapply(tables, function(tab) {
        if (is.null(tab) || !nrow(tab)) return(NULL)
        cbind(d[rep.int(1L, nrow(tab)), id, drop = FALSE], effort = effort, tab)
      })
    })
  })
  results <- do.call(c, results)
  out <- lapply(names(results[[1L]]), function(name) {
    tab <- do.call(rbind, lapply(results, function(x) x[[name]]))
    if (is.null(tab)) return(data.frame())
    rownames(tab) <- NULL
    tab
  })
  names(out) <- names(results[[1L]])
  out$variance_selection <- oe_diag_variance_selection(out$variance, out$variance_score_draws, bootstrap)
  c(out, oe_diag_level_share(draws, min_n = if (length(calibration) < 10L) 2L else 10L))
}

oe_diagnostic_checks <- function() {
  ages <- c(1L, 3L, 4L, 7L)
  x <- oe_diag_ar1_data(data.frame(
    age = ages, year = 1L, replicate = 1L, centered = c(0.1, -0.6, 0.8, -0.2),
    cell_sd = c(0.2, 0.7, 1.4, 0.5)
  ))
  stopifnot(isTRUE(all.equal(oe_diag_ar1_loss(x, 0),
    -stats::dnorm(x$d$centered, sd = x$d$cell_sd, log = TRUE), tolerance = 1e-12)))
  for (rho in c(-0.6, 0.6)) {
    r <- rho^abs(outer(ages, ages, "-"))
    chol_r <- chol(r)
    z <- backsolve(chol_r, x$d$z, transpose = TRUE)
    expected <- sum(log(diag(chol_r)) + log(x$d$cell_sd) + (log(2 * pi) + z^2) / 2)
    stopifnot(isTRUE(all.equal(sum(oe_diag_ar1_loss(x, rho)), expected, tolerance = 1e-12)))
  }
  sparse <- expand.grid(year = 1L, age = ages, replicate = 1:4)
  sparse$component <- "RV"
  sparse$scenario <- "baseline"
  sparse$bin <- "modeled"
  sparse$truth <- 10
  sparse$selectivity <- 1
  sparse$expected_encounters <- 5
  sparse$obs <- NA_real_
  sparse$se <- NA_real_
  out <- oe_diagnostics(sparse, calibration = 1:2, validation = 3:4, bootstrap = 0L)
  stopifnot(all(out$age_year$population_id == "seed17"), all(out$age_year$n_missing == 2L),
    all(is.na(out$age_year$zero_rate)), all(out$variance$n_scored == 0L),
    out$covariance_scores$n_obs == 0L, is.na(out$covariance_scores$rho))

  # A zero contributes to arithmetic moments and support, while an unavailable
  # observation contributes only to the missingness denominator.
  moments <- oe_diag_age_year(data.frame(obs = c(0, 2, NA), truth = 2,
    log_error = c(NA, 0, NA), selectivity = 1, expected_encounters = 2,
    n_caught_at_age = c(0, 4, NA), n_aged_at_age = c(0, 3, NA)))
  stopifnot(moments$n_observed == 2L, moments$mean_obs == 1,
    abs(moments$sd_obs - sqrt(2)) < 1e-12, abs(moments$mcse_obs - 1) < 1e-12,
    moments$zero_rate == 0.5, moments$missing_rate == 1 / 3,
    moments$mean_n_caught_at_age == 2, moments$n_aged_at_age_draws == 2L)

  paired <- expand.grid(age = c(1L, 3L), year = 1L, replicate = 1:4)
  paired$truth <- 2
  paired$obs <- c(0, 1, 1, 2, 3, 3, 2, NA)
  paired$log_error <- NA_real_
  positive <- is.finite(paired$obs) & paired$obs > 0
  paired$log_error[positive] <- log(paired$obs[positive] / paired$truth[positive])
  pairs <- rbind(oe_diag_pairs(paired, "original"), oe_diag_pairs(paired, "log"))
  stopifnot(identical(pairs$scale, c("original", "log")), all(pairs$n_draws == 4L),
    identical(as.integer(pairs$n), c(3L, 2L)), abs(pairs$covariance[1L] - 1.5) < 1e-12)

  shared <- expand.grid(age = 1:3, replicate = 1:5)
  shared$year <- 1L
  shared$population_id <- "seed17"
  shared$om <- "coupled"
  shared$effort <- "standard"
  shared$bin <- "modeled"
  shared$component <- "RV"
  shared$truth <- 5
  shared$obs <- c(1, 3, 4, 3, 4, 8, 8, 5, 13, 10, 7, NA, 11, 7, 18)
  shared$oracle_obs <- shared$obs
  shared$oracle_obs[shared$replicate == 5L & shared$age == 2L] <- 0
  adjusted <- oe_diag_level_share(shared, min_n = 3L)
  stopifnot(nrow(adjusted$level_share) == 4L, all(adjusted$level_share$n == 3L),
    all(adjusted$total_cv$n == 4L),
    all(adjusted$level_share$correlation_level > 0.9),
    adjusted$level_share$correlation_share[1L] < 0)
  shared$oracle_obs <- NULL
  adjusted <- oe_diag_level_share(shared, min_n = 3L)
  stopifnot(nrow(adjusted$level_share) == 2L, all(adjusted$level_share$n == 4L))

  paired$centered <- paired$log_error
  paired$cell_sd <- 1
  covariance <- oe_diag_covariance(paired, paired, bootstrap = 0L)$covariance_scores
  stopifnot(!covariance$rho_sign_identified, covariance$rho >= 0)

  # The supplied SD is 0.5. Training errors give a constant SD of 1 and
  # a multiplier of 2. The invalid held-out SE is excluded for all models.
  training <- data.frame(replicate = rep(1:2, each = 2L), obs = 2,
    se = 2 * sqrt(expm1(0.25)), centered = c(-1, 1, -1, 1))
  heldout <- transform(training, replicate = replicate + 2L,
    centered = c(-0.5, 1, -1.5, 0.75))
  heldout$se[4L] <- 0
  precision <- oe_diag_rv_precision_scores(training, heldout, bootstrap = 0L)
  scaled <- precision[precision$model == "supplied_scaled", ]
  stopifnot(all(precision$status == "ok"), all(precision$n_calibration == 4L),
    all(precision$n_scored == 3L), abs(scaled$multiplier - 2) < 1e-12,
    abs(scaled$delta) < 1e-12,
    precision$coverage[precision$model == "supplied_fixed"] == 0.25)

  # Two fish per tow, ten standard tows and theta = 0.5 give variance 0.4.
  # Changing a held-out estimate must not change its predicted SD.
  standard <- data.frame(effort = "standard", year = 1L, age = 1L, truth = 2,
    n_equivalent = 1000, n_primary = 10L,
    arith_centered = sqrt(0.4) * c(-1, 1, -1, 1))
  standard$obs <- standard$truth + standard$arith_centered
  training <- rbind(standard, transform(standard, effort = "half", arith_centered = 100))
  heldout <- data.frame(year = 1L, age = 1L, truth = 2, n_equivalent = 1000,
    n_primary = c(5L, 10L, 20L), obs = c(0, 2, 999))
  law <- oe_diag_cadigan(training, heldout)
  stopifnot(law$n_standard_calibration == 4L, all(law$test_ok),
    abs(law$oracle$theta - 0.5) < 1e-6,
    abs(law$calibration$theta - 0.5) < 1e-6,
    max(abs(law$calibration$original$sigma^2 - c(0.8, 0.4, 0.2))) < 1e-6,
    max(abs(law$calibration$log$sigma^2 - log1p(c(0.8, 0.4, 0.2) / 4))) < 1e-6)
  stopifnot(identical(oe_diag_age_pair(c(1L, 2L, 12L, 13L, 14L)),
    c("pair_0", "pair_0", "pair_5", "single_13", "single_14")))

  # Repeated draws may have different truth across population IDs.
  another <- transform(sparse, population_id = "seed19", truth = 20)
  sparse$population_id <- "seed17"
  out <- oe_diagnostics(rbind(sparse, another), calibration = 1:2,
    validation = 3:4, bootstrap = 0L)
  stopifnot(setequal(out$age_year$population_id, c("seed17", "seed19")),
    all(out$age_year$mean_truth[out$age_year$population_id == "seed19"] == 20))
  invisible(TRUE)
}

oe_fitted_diagnostics <- function(predictions) {
  if (!nrow(predictions)) return(data.frame())
  if (!"population_id" %in% names(predictions)) predictions$population_id <- "seed17"
  id <- c("population_id", "component", "treatment", "replicate", "age")
  required <- c(id, "log_error", "log_truth_error", "log_prediction_error", "std_res")
  stopifnot(all(required %in% names(predictions)))
  d <- as.data.frame(predictions)[required]
  groups <- split(d, interaction(d[id], drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(x) {
    keep <- is.finite(x$log_error) & is.finite(x$log_truth_error) &
      is.finite(x$log_prediction_error) & is.finite(x$std_res)
    cbind(x[1L, id], data.frame(n_years = sum(keep),
      truth_error_sd = oe_diag_sd(x$log_truth_error[keep]),
      fitted_error_sd = oe_diag_sd(x$log_error[keep]),
      standardized_sd = oe_diag_sd(x$std_res[keep]),
      mean_prediction_error = oe_diag_mean(x$log_prediction_error[keep])))
  }))
  rownames(out) <- NULL
  out
}
