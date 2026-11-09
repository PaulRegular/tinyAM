# Optional assessment stage after the fixed-population sampling comparison.

.oe_fit_series <- function() c(catch = "Gillnet", index = "RV")

.oe_fit_key <- function(data) paste(data$year, data$age, sep = "/")

oe_assessment_inputs <- function(snapshots) {
  if (!length(snapshots) || !all(vapply(snapshots, function(x) {
    is.data.frame(x$population) && is.data.frame(x$annual)
  }, logical(1)))) {
    stop("Supply the Gillnet truth snapshots with population and annual tables.")
  }
  stack <- function(name) {
    do.call(rbind, lapply(snapshots, function(x) as.data.frame(x[[name]])))
  }
  list(population = stack("population"), annual = stack("annual"))
}

.oe_fit_data <- function(truth) {
  pop <- as.data.frame(truth$population)
  keys <- c("year", "age")
  quantities <- c("N", "Z", "F", "M", "weight", "maturity")
  if (!all(c(keys, quantities) %in% names(pop)) || nrow(pop) != 400L ||
      !setequal(pop$year, 1:20) || !setequal(pop$age, 1:20) ||
      anyDuplicated(pop[keys]) || any(!is.finite(as.matrix(pop[quantities])))) {
    stop("Assessment truth must cover the complete 20-year, 20-age population grid.")
  }
  annual <- as.data.frame(truth$annual)
  if (!all(c("year", "SSB", "recruitment") %in% names(annual)) || nrow(annual) != 20L ||
      !setequal(annual$year, 1:20) || any(!is.finite(annual$SSB) | !is.finite(annual$recruitment))) {
    stop("Supply annual generating SSB and recruitment for all 20 years.")
  }
  tolerance <- 32 * .Machine$double.eps * max(1, max(abs(pop$Z)))
  if (any(abs(pop$M + pop$F - pop$Z) > tolerance)) {
    stop("Realized M must match Z minus Gillnet F.")
  }
  # The assessment shares Gillnet F over 14+. Preserve generating background
  # survival at the true age mix inside that group, including retained ages.
  older <- pop$age >= 14L
  plus <- stats::aggregate(cbind(survivors = pop$N[older] * exp(-pop$M[older]),
    abundance = pop$N[older]), list(year = pop$year[older]), sum)
  if (any(plus$abundance <= 0)) stop("Realized plus-group M requires fish aged 14+ in every year.")
  plus_m <- -log(plus$survivors / plus$abundance)
  supplied <- c(pop$M[!older], plus_m)
  if (any(!is.finite(supplied) | supplied < -tolerance)) {
    stop("Realized M is materially negative or unavailable in a modeled cell.")
  }
  n_floored <- sum(supplied < tolerance)
  m <- pmax(pop$M, tolerance)
  m[older] <- pmax(plus_m, tolerance)[match(pop$year[older], plus$year)]
  grid <- pop[keys]
  paired_age <- cut(pmin(grid$age, 14L), c(0, seq(2, 12, by = 2), 13, 14))
  data <- tinyAM::as_tam_data(list(
    catch = transform(grid, obs = NA_real_, unit = "thousand fish",
      age_band = cut(age, c(0, 3, 5, 10, 13, Inf)),
      paired_age = paired_age),
    index = transform(grid[pop$age <= 14L, ], obs = NA_real_, survey = "RV",
      samp_time = 0, unit = "thousand fish", paired_age = paired_age[pop$age <= 14L]),
    weight = transform(grid, obs = pop$weight, unit = "kg"),
    maturity = transform(grid, obs = pop$maturity, spawn_time = 0),
    pop_covariates = transform(grid, M_realized = m)
  ))
  list(data = data, n_mortality_floored = n_floored, mortality_floor = tolerance)
}

.oe_fit_check_sample <- function(sample, type, expected) {
  group <- if (type == "catch") "fleet" else "survey"
  series <- .oe_fit_series()[[type]]
  if (!is.data.frame(sample) || !all(c("year", "age", "obs", "truth") %in% names(sample)) ||
      anyDuplicated(.oe_fit_key(sample)) ||
      !setequal(.oe_fit_key(sample), .oe_fit_key(expected)) ||
      (group %in% names(sample) && !all(sample[[group]] %in% series)) ||
      ("component" %in% names(sample) && !all(sample$component %in% series))) {
    stop("Modeled ", series, " sample keys must match the input grid.")
  }
}

.oe_fit_precision <- function(candidate) {
  candidate %in% c("supplied_fixed", "supplied_scaled", "cadigan_calibration")
}

.oe_fit_inputs <- function(data, training, validation, rv_fallback, candidate) {
  for (type in names(.oe_fit_series())) {
    input <- data[[type]]
    modeled <- input[input$age <= 14L, ]
    .oe_fit_check_sample(training[[type]], type, modeled)
    .oe_fit_check_sample(validation[[type]], type, modeled)
    row <- match(.oe_fit_key(input), .oe_fit_key(training[[type]]))
    # The direct estimators already combined ages 14+. Keep the complete raw
    # grid for preparation, but put the entire observation at age 14.
    data[[type]]$obs <- training[[type]]$obs[row]
  }

  if (!.oe_fit_precision(candidate)) return(list(data = data, n_rv_fallback = 0L))
  if (length(rv_fallback) != 1L || !is.finite(rv_fallback) || rv_fallback <= 0) {
    stop("RV precision treatment requires one positive fallback log SD from calibration draws.")
  }
  rv <- training$index
  observed <- is.finite(rv$obs) & rv$obs > 0
  sd <- rep(1, nrow(rv))
  if (candidate == "cadigan_calibration") {
    if (!is.numeric(rv$cadigan_log_sd) || length(rv$cadigan_log_sd) != nrow(rv)) {
      stop("Cadigan calibration requires an independently predicted numeric log SD per RV row.")
    }
    sd[observed] <- rv$cadigan_log_sd[observed]
  } else {
    if (!is.numeric(rv$se) || length(rv$se) != nrow(rv)) {
      stop("The modeled RV sample must include a numeric SE per row.")
    }
    known <- observed & is.finite(rv$se) & rv$se > 0
    log_cv <- log(rv$se[known]) - log(rv$obs[known])
    sd[known] <- sqrt(pmax(2 * log_cv, 0) + log1p(exp(-abs(2 * log_cv))))
    sd[observed & !known] <- NA_real_
  }
  unsupported <- observed & (!is.finite(sd) | sd <= 0)
  sd[unsupported] <- rv_fallback
  # Supplied offsets are evaluated on every row. Excluded observations use 1.
  data$index$rv_log_sd <- sd[match(.oe_fit_key(data$index), .oe_fit_key(rv))]
  data$index$rv_sd_scale <- 1
  list(data = data, n_rv_fallback = sum(unsupported))
}

.oe_fit_settings <- function(component, candidate) {
  sd_formula <- switch(candidate,
    paired_age = ~ paired_age,
    iid_age = ~ iid(age),
    free_age = ~ factor(age),
    cadigan_calibration = ~ 0,
    supplied_fixed = ~ 0,
    supplied_scaled = ~ 0 + rv_sd_scale,
    ~ 1)
  catch_settings <- list(sd_form = if (component == "Gillnet") sd_formula else ~ 1)
  index_settings <- list(q_form = ~ logistic(age),
    sd_form = if (component == "RV") sd_formula else ~ 1)
  if (.oe_fit_precision(candidate)) index_settings$sd_supplied <- ~ rv_log_sd
  list(
    ages = 1:14,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "rw"),
    F_settings = list(process = "off", mu_form = ~ factor(year) + age_band, mean_ages = 5:10),
    M_settings = list(process = "off", mu_supplied = ~ M_realized, mean_ages = 5:10),
    catch_settings = catch_settings,
    index_settings = index_settings
  )
}

.oe_fit_join <- function(estimate, truth, keys) {
  d <- merge(estimate, truth, by = keys, all.x = TRUE)
  valid <- is.finite(d$lwr) & is.finite(d$upr) & is.finite(d$truth)
  d$covered <- ifelse(valid, d$lwr <= d$truth & d$upr >= d$truth, NA)
  d$error <- d$est - d$truth
  d$relative_error <- ifelse(d$truth > 0, d$error / d$truth, NA_real_)
  d
}

.oe_fit_fbar <- function(year, F, N) {
  d <- stats::aggregate(cbind(FN = ifelse(N == 0, 0, F * N), N = N), list(year = year), sum)
  data.frame(year = d$year, value = d$FN / d$N)
}

.oe_fit_population <- function(fit, truth) {
  pop <- fit$pop
  adult <- truth$population[truth$population$age %in% 5:10, ]
  fbar <- with(adult, .oe_fit_fbar(year, F, N))
  annual <- list(ssb = pop$ssb, recruitment = pop$recruitment, F_bar = pop$F_bar)
  targets <- list(
    ssb = data.frame(year = truth$annual$year, truth = truth$annual$SSB),
    recruitment = data.frame(year = truth$annual$year, truth = truth$annual$recruitment),
    F_bar = data.frame(year = fbar$year, truth = fbar$value)
  )
  annual <- do.call(rbind, lapply(names(annual), function(quantity) {
    estimate <- annual[[quantity]][c("year", "est", "lwr", "upr", "unit")]
    d <- .oe_fit_join(estimate, targets[[quantity]], "year")
    d$quantity <- quantity
    d
  }))
  # The age-14+ fitted rate spans multiple biological ages. Compare single
  # ages on like support, and show the adult fleet average separately.
  gillnet_est <- pop$F[pop$F$age < 14L, c("year", "age", "est", "lwr", "upr", "unit")]
  gillnet_truth <- data.frame(year = truth$population$year, age = truth$population$age,
    truth = truth$population$F)
  list(annual = annual, gillnet_F = .oe_fit_join(gillnet_est, gillnet_truth, c("year", "age")))
}

.oe_log_error <- function(value, reference) {
  valid <- is.finite(value) & value > 0 & is.finite(reference) & reference > 0
  error <- rep(NA_real_, length(value))
  error[valid] <- log(value[valid]) - log(reference[valid])
  error
}

.oe_fit_predictions <- function(fit, training, validation, truth) {
  series <- .oe_fit_series()
  abundance <- truth$population
  abundance$age <- pmin(abundance$age, 14L)
  abundance <- stats::aggregate(N ~ year + age, abundance, sum)
  do.call(rbind, lapply(names(series), function(type) {
    fitted <- fit$obs_pred[[type]]
    train <- training[[type]]
    test <- validation[[type]]
    trained <- match(.oe_fit_key(fitted), .oe_fit_key(train))
    tested <- match(.oe_fit_key(fitted), .oe_fit_key(test))
    if (anyNA(trained) || anyNA(tested) || nrow(fitted) != nrow(train) ||
        nrow(fitted) != nrow(test) ||
        !isTRUE(all.equal(fitted$obs, train$obs[trained], tolerance = 0))) {
      stop("Fitted and held-out ", series[[type]], " rows must match the study inputs.")
    }
    d <- fitted[c("year", "age", "obs", "pred", "sd", "std_res")]
    d$component <- series[[type]]
    d$truth <- train$truth[trained]
    d$holdout_obs <- test$obs[tested]
    d$log_error <- .oe_log_error(d$obs, d$pred)
    d$log_truth_error <- .oe_log_error(d$obs, d$truth)
    d$log_prediction_error <- .oe_log_error(d$pred, d$truth)
    d$holdout_log_error <- .oe_log_error(d$holdout_obs, d$pred)
    d$q <- d$q_lwr <- d$q_upr <- d$q_truth <- NA_real_
    d$q_covered <- NA
    if (type == "index") {
      d[c("q", "q_lwr", "q_upr")] <- fitted[c("q", "q_lwr", "q_upr")]
      n <- abundance$N[match(.oe_fit_key(d), .oe_fit_key(abundance))]
      d$q_truth <- ifelse(n > 0, d$truth / n, NA_real_)
      available <- is.finite(d$q_lwr) & is.finite(d$q_upr) & is.finite(d$q_truth)
      d$q_covered[available] <- d$q_lwr[available] <= d$q_truth[available] &
        d$q_upr[available] >= d$q_truth[available]
    }
    d
  }))
}

.oe_fit_holdout <- function(predictions) {
  do.call(rbind, lapply(split(predictions, predictions$component), function(d) {
    eligible <- is.finite(d$holdout_obs) & d$holdout_obs > 0
    n_eligible <- sum(eligible)
    n_scored <- sum(is.finite(d$holdout_log_error[eligible]))
    # Every treatment must use the same held-out rows. Keep coverage counts
    # if a prediction is unusable, but do not score a reduced subset.
    error <- if (n_eligible > 0 && n_eligible == n_scored) d$holdout_log_error[eligible] else NA_real_
    data.frame(component = d$component[1L], n_eligible = n_eligible, n_scored = n_scored,
      mean_log_error = mean(error), log_MAE = mean(abs(error)), log_RMSE = sqrt(mean(error^2)))
  }))
}

.oe_fit_pair <- function(data, pair, truth, treatment, component, candidate, rv_fallback) {
  warnings <- character()
  attempt <- data.frame(replicate = pair$replicate, holdout_replicate = pair$holdout_replicate,
    treatment = treatment, candidate = candidate,
    converged = FALSE, scored = FALSE, optimizer_code = NA_integer_, pd_hessian = NA,
    max_gradient = NA_real_, objective = NA_real_, n_rv_fallback = NA_integer_,
    rv_multiplier = NA_real_, gillnet_sd_young = NA_real_,
    gillnet_sd_middle = NA_real_, gillnet_sd_old = NA_real_, warnings = "", error = "")
  out <- list()
  withCallingHandlers(tryCatch({
    input <- .oe_fit_inputs(data, pair$training, pair$validation, rv_fallback, candidate)
    attempt$n_rv_fallback <- input$n_rv_fallback
    fit <- do.call(tinyAM::fit_tam, c(list(data = input$data, silent = TRUE),
      .oe_fit_settings(component, candidate)))
    checks <- tinyAM::check_tam(fit)
    out$checks <- checks
    attempt$converged <- checks$is_converged
    attempt$max_gradient <- checks$max_gradient
    if (!is.null(fit$opt$objective)) attempt$objective <- fit$opt$objective
    if (!is.null(checks$optimizer_code)) attempt$optimizer_code <- checks$optimizer_code
    if (!is.null(checks$pd_hessian)) attempt$pd_hessian <- checks$pd_hessian
    if (checks$is_converged) {
      gillnet <- subset(fit$obs_pred$catch, obs > 0)
      attempt$gillnet_sd_young <- mean(gillnet$sd[gillnet$age <= 3L])
      attempt$gillnet_sd_middle <- mean(gillnet$sd[gillnet$age %in% 4:10])
      attempt$gillnet_sd_old <- mean(gillnet$sd[gillnet$age >= 11L])
      if (.oe_fit_precision(candidate)) {
        d <- fit$obs_pred$index
        row <- which(d$obs > 0)[1L]
        rv <- input$data$index
        attempt$rv_multiplier <- d$sd[row] /
          rv$rv_log_sd[match(.oe_fit_key(d[row, ]), .oe_fit_key(rv))]
      }
      out <- c(out, .oe_fit_population(fit, truth))
      out$predictions <- .oe_fit_predictions(fit, pair$training, pair$validation, truth)
      out$holdout <- .oe_fit_holdout(out$predictions)
      attempt$scored <- all(out$holdout$n_eligible > 0L &
        out$holdout$n_eligible == out$holdout$n_scored & is.finite(out$holdout$log_RMSE))
    }
  }, error = function(e) {
    attempt$error <<- conditionMessage(e)
  }), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  attempt$warnings <- paste(warnings, collapse = " | ")
  out$attempt <- attempt
  out$warnings <- warnings
  for (table in intersect(c("annual", "gillnet_F", "predictions", "holdout"), names(out))) {
    out[[table]]$replicate <- pair$replicate
    out[[table]]$holdout_replicate <- pair$holdout_replicate
    out[[table]]$treatment <- treatment
  }
  out
}

oe_assessment <- function(truth, pairs, candidates, rv_fallback = NULL,
                          checkpoint = function(replicate, treatment, run) run()) {
  allowed <- list(Gillnet = c("paired_age", "iid_age", "free_age"),
    RV = c("paired_age", "iid_age", "free_age", "cadigan_calibration",
      "supplied_fixed", "supplied_scaled"))
  components <- names(candidates)
  if (length(candidates) && (!is.character(candidates) ||
      length(components) != length(candidates) || anyDuplicated(components) ||
      any(!components %in% names(allowed)) ||
      !all(vapply(components, function(name) candidates[[name]] %in% allowed[[name]], logical(1))))) {
    stop("Supply at most one supported candidate per component, named Gillnet or RV.")
  }
  ids <- vapply(pairs, `[[`, numeric(1), "replicate")
  holdout <- vapply(pairs, `[[`, numeric(1), "holdout_replicate")
  if (length(pairs) != 10L || any(!is.finite(c(ids, holdout))) ||
      anyDuplicated(c(ids, holdout))) {
    stop("Supply ten training and ten independent holdout draws with distinct numeric IDs.")
  }
  input <- .oe_fit_data(truth)
  treatments <- c("baseline", if (length(candidates)) paste(tolower(components), candidates, sep = "_"))
  components <- c("", components)
  candidates <- c("constant", unname(candidates))
  results <- lapply(pairs, function(pair) {
    lapply(seq_along(treatments), function(i) {
      checkpoint(pair$replicate, treatments[i], function() {
        .oe_fit_pair(input$data, pair, truth, treatments[i], components[i], candidates[i], rv_fallback)
      })
    })
  })
  results <- unlist(results, recursive = FALSE)
  bind <- function(table) {
    out <- as.data.frame(do.call(rbind, lapply(results, `[[`, table)))
    if (!is.null(truth$population_id)) out$population_id <- rep(truth$population_id, nrow(out))
    out
  }
  tables <- c("attempt", "holdout", "annual", "gillnet_F", "predictions")
  out <- stats::setNames(lapply(tables, bind), c("attempts", tables[-1L]))
  expected <- expand.grid(replicate = ids, treatment = treatments)
  fit_key <- function(x) paste(x$replicate, x$treatment)
  if (anyDuplicated(fit_key(out$attempts)) ||
      !setequal(fit_key(out$attempts), fit_key(expected))) {
    stop("Assessment checkpoints must retain one attempt per training draw and treatment.")
  }
  success <- split(out$attempts$replicate[out$attempts$scored],
    out$attempts$treatment[out$attempts$scored])
  common <- sort(Reduce(intersect, success[treatments]))
  out$common_replicates <- common
  out$comparison_allowed <- length(treatments) > 1L && length(common) >= 8L
  if (out$comparison_allowed) {
    out$comparison <- lapply(out[c("holdout", "annual", "gillnet_F", "predictions")],
      function(x) x[x$replicate %in% common, , drop = FALSE])
    common_intervals <- function(x, keys, column) {
      group <- interaction(x[c("replicate", keys)], drop = TRUE)
      available <- ave(as.integer(!is.na(x[[column]])), group, FUN = function(values) {
        length(values) == length(treatments) && all(values == 1L)
      }) == 1L
      x[[column]][!available] <- NA
      x
    }
    out$comparison$annual <- common_intervals(out$comparison$annual, c("year", "quantity"), "covered")
    out$comparison$gillnet_F <- common_intervals(out$comparison$gillnet_F, c("year", "age"), "covered")
    out$comparison$predictions <- common_intervals(out$comparison$predictions,
      c("component", "year", "age"), "q_covered")
  }
  out$mortality <- input[c("n_mortality_floored", "mortality_floor")]
  out
}
