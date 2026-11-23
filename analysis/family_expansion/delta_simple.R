# Exact observation-family benchmarks around default SimSurvey mortality.
# Run through simulate_delta.R; SimSurvey remains a study-only dependency.

delta_share_exact_draw <- function(mu, setting, scheme, seed) {
  withr::with_seed(seed, {
    if (scheme == "poisson-link") {
      lambda <- log(2) * mu / setting$y50
      positive <- -expm1(-lambda)
      seen <- matrix(stats::runif(length(mu)), nrow(mu)) < positive
      error <- matrix(stats::rnorm(length(mu), -setting$sd^2 / 2, setting$sd), nrow(mu))
      obs <- (mu / positive) * exp(error) * seen
    } else {
      # Independent compound-Poisson/Gamma construction, not tinyAM's simulator.
      power <- setting$tweedie_power
      phi <- setting$tweedie_phi
      lambda <- mu^(2 - power) / (phi * (2 - power))
      shape <- (2 - power) / (power - 1)
      scale <- phi * (power - 1) * mu^(power - 1)
      encounters <- stats::rpois(length(mu), as.vector(lambda))
      obs <- matrix(stats::rgamma(length(mu), shape = encounters * shape,
        scale = as.vector(scale)), nrow(mu), ncol(mu))
    }
    stopifnot(all(is.finite(obs) & obs >= 0))
    list(obs = obs, p_zero = exp(-lambda))
  })
}

delta_share_ablation_support <- function(truth, calibration) {
  rows <- lapply(c("catch", "index"), function(component) {
    s <- calibration$values[calibration$values$component == component, ]
    lambda <- log(2) * as.vector(truth$mean[[component]]) / s$y50
    slope <- 1 - lambda * exp(-lambda) / (-expm1(-lambda))
    small <- lambda < 1e-3
    slope[small] <- lambda[small] / 2 - lambda[small]^2 / 12 + lambda[small]^4 / 720
    data.frame(component = component, truth$reference[c("year", "age")],
      intensity = lambda, positive_mean_sensitivity = slope)
  })
  list(applicable_scheme = "poisson-link", fixed_y50 = TRUE, fixed_positive_sd = TRUE,
    warning = paste("Positive-only ablation must retain known y50 and positive SD.",
      "Local mean sensitivity approaches zero at low intensity;",
      "this diagnostic does not replace the fitted Hessian/gradient checks."),
    rows = do.call(rbind, rows))
}

# Small exact-family benchmarks use the same default SimSurvey population
# law as the spatial study, without generating a spatial grid or fisheries.
delta_share_truth <- function(id, truth = NULL) {
  stopifnot(length(id) == 1L, is.finite(id), id > 0, id == as.integer(id))
  seed <- 104000L + 100L * id
  years <- 1:30
  ages <- 1:20
  share <- 0.10
  keys <- expand.grid(year = years, age = ages)
  key <- function(d) paste(d$year, d$age, sep = "/")
  if (is.null(truth)) {
    if (!requireNamespace("SimSurvey", quietly = TRUE)) {
      stop("Load the study's pinned SimSurvey or supply a population/reference bundle.")
    }
    pop <- withr::with_seed(seed, SimSurvey::sim_abundance(years = years))
    terminal <- length(ages)
    for (t in seq_along(years)[-1L]) {
      pop$N[terminal, t] <- pop$N[terminal, t] +
        pop$N[terminal, t - 1L] * exp(-pop$Z[terminal, t - 1L])
    }
    # Adding the 20+ carryover changes abundance. Refresh weight and mature
    # biomass using the pinned biological functions before computing truth.
    pop <- getFromNamespace(".population_biology", "SimSurvey")(pop, required = TRUE)
    reference <- data.frame(keys, N = as.vector(t(pop$N)) / 1000,
      Z = as.vector(t(pop$Z)), ssb = as.vector(t(pop$S)) / 1000)
    data <- list(
      catch = transform(keys, obs = NA_real_, unit = "thousand fish"),
      index = transform(keys, obs = NA_real_, unit = "thousand fish",
        survey = "RV", samp_time = 0),
      weight = transform(keys, obs = as.vector(t(pop$B / pop$N)), unit = "kg"),
      maturity = transform(keys, obs = as.vector(t(pop$S / pop$B)), spawn_time = 0))
  } else {
    stopifnot(is.list(truth), is.data.frame(truth$reference),
      is.list(truth$data), all(c("weight", "maturity") %in% names(truth$data)))
    reference <- truth$reference
    stopifnot(!anyDuplicated(key(reference)), setequal(key(reference), key(keys)),
      all(c("N", "Z", "ssb") %in% names(reference)))
    reference <- reference[match(key(keys), key(reference)), , drop = FALSE]
    data <- truth$data
  }
  stopifnot(all(is.finite(reference$N) & reference$N > 0),
    all(is.finite(reference$Z) & reference$Z > 0),
    all(is.finite(reference$ssb) & reference$ssb >= 0))
  reference$F_realized <- share * reference$Z
  reference$M_realized <- (1 - share) * reference$Z
  reference$catch_expected <- with(reference, N * share * (-expm1(-Z)))
  reference$index_expected <- reference$N
  reference$F_reference <- reference$F_realized
  reference$M_reference <- reference$M_realized
  reference$deaths_expected <- with(reference, N * (-expm1(-Z)))
  data$pop_covariates <- transform(reference[c("year", "age", "M_realized")],
    known_m = M_realized)
  matrices <- lapply(c("N", "Z", "F_realized", "M_realized"), function(name) {
    matrix(reference[[name]], length(years), length(ages),
      dimnames = list(year = years, age = ages))
  })
  names(matrices) <- c("N", "Z", "F", "M")
  means <- list(catch = matrix(reference$catch_expected, length(years), length(ages),
    dimnames = dimnames(matrices$N)), index = matrices$N)
  predicted <- matrices$N[-length(years), ] * exp(-matrices$Z[-length(years), ])
  expected <- predicted[, -length(ages), drop = FALSE]
  expected[, ncol(expected)] <- rowSums(predicted[, tail(seq_along(ages), 2L)])
  actual <- matrices$N[-1L, -1L, drop = FALSE]
  # Supplied spatial references can be integer balanced, unlike this light
  # benchmark. Do not silently accept that rounding as an exact cohort law.
  stopifnot(max(abs(actual - expected) / pmax(1, abs(expected))) < 1e-12)
  audit <- data.frame(id = id, driver = "low_f", valid = TRUE,
    target_share = share, realized_share = sum(reference$catch_expected) /
      sum(reference$deaths_expected), nonpositive_m = 0L,
    minimum_m = min(reference$M_realized),
    max_abs_balance = max(abs(reference$Z - reference$F_realized - reference$M_realized)),
    max_abs_removal_error = max(abs(with(reference,
      N * F_realized * (-expm1(-Z)) / Z) - reference$catch_expected)))
  c(list(id = as.integer(id), driver = "low_f", effort = "standard", share = share,
    years = years, ages = ages, seed = seed, mean = means, reference = reference,
    data = data, ssb = rowSums(matrix(reference$ssb, length(years))),
    recruitment = matrices$N[, 1L], audit = audit), matrices)
}

delta_share_calibration <- function(truths = NULL) {
  ids <- 101:106
  if (is.null(truths)) {
    truths <- lapply(ids, delta_share_truth)
  } else {
    stopifnot(is.list(truths), length(truths) == length(ids))
    truths <- lapply(seq_along(ids), function(i) delta_share_truth(ids[i], truths[[i]]))
  }
  target <- c(catch = .10, index = .15)
  positive_sd <- c(catch = .25, index = .35)
  power <- 1.5
  values <- lapply(names(target), function(component) {
    mu <- do.call(rbind, lapply(truths, function(x) x$mean[[component]]))
    older <- as.vector(mu[, -1L])
    y50 <- exp(stats::uniroot(function(log_scale) {
      mean(exp(-log(2) * older / exp(log_scale))) - target[[component]]
    }, range(log(older)) + c(-20, 20), tol = 1e-10)$root)
    phi <- exp(stats::uniroot(function(log_phi) {
      mean(exp(-older^(2 - power) / (exp(log_phi) * (2 - power)))) - target[[component]]
    }, log(mean(older^(2 - power)) / (2 - power)) + c(-20, 20), tol = 1e-10)$root)
    data.frame(component = component, target = unname(target[[component]]),
      y50 = y50, sd = unname(positive_sd[[component]]),
      tweedie_power = power, tweedie_phi = phi, tweedie_sd = sqrt(phi),
      poisson_zero = mean(exp(-log(2) * older / y50)),
      tweedie_zero = mean(exp(-older^(2 - power) / (phi * (2 - power)))))
  })
  list(ids = ids, values = do.call(rbind, values),
    basis = "Six independent default-Z/R truth histories; ages 2:20+; no monitoring draws",
    population_seeds = 104000L + 100L * ids, share = .10,
    years = 1:30, ages = 1:20, q = 1,
    simsurvey_revision = "29a1b03350c274f88dac108fb95f9446eac79233")
}

delta_share_exact <- function(id, scheme = c("poisson-link", "tweedie"), calibration,
                              truth = NULL) {
  scheme <- match.arg(scheme)
  stopifnot(is.data.frame(calibration$values),
    setequal(calibration$values$component, c("catch", "index")),
    !anyDuplicated(calibration$values$component))
  truth <- delta_share_truth(id, truth)
  seed_base <- c(`poisson-link` = 184000L, tweedie = 188000L)[[scheme]]
  sampling <- probabilities <- list()
  key <- function(d) paste(d$year, d$age, sep = "/")
  make_draw <- function(draw) {
    data <- truth$data
    for (j in seq_along(truth$mean)) {
      component <- names(truth$mean)[j]
      setting <- calibration$values[calibration$values$component == component, ]
      seed <- seed_base + id + 10000L * j + 100000L * (draw == "holdout")
      result <- delta_share_exact_draw(truth$mean[[component]], setting, scheme, seed)
      row <- match(key(data[[component]]), key(truth$reference))
      stopifnot(!anyNA(row))
      data[[component]]$obs <- as.vector(result$obs)[row]
      if (component == "index") {
        data[[component]]$survey <- "RV"
        data[[component]]$samp_time <- 0
      }
      name <- paste(draw, component, sep = "_")
      probabilities[[name]] <<- data.frame(component = component, draw = draw,
        truth$reference[c("year", "age")], p_zero = as.vector(result$p_zero))
      sampling[[name]] <<- data.frame(component = component, draw = draw,
        year = truth$years, n_aged = NA_integer_, planned_n = NA_integer_,
        reported_total = rowSums(result$obs), truth_total = rowSums(truth$mean[[component]]),
        status = "exact_family_draw", seed = seed)
    }
    data <- tinyAM::as_tam_data(data)
    tinyAM::check_tam_data(data, quiet = TRUE)
    data
  }
  data <- make_draw("training")
  holdout <- make_draw("holdout")
  exact_settings <- calibration$values
  exact_settings$scheme <- scheme
  exact_settings$family <- if (scheme == "poisson-link") "delta_lognormal" else "tweedie"
  exact_settings$delta_type <- if (scheme == "poisson-link") "poisson-link" else NA_character_
  exact_settings$Q_truth <- if (scheme == "poisson-link") 0 else NA_real_
  exact_settings$power_truth <- if (scheme == "tweedie") exact_settings$tweedie_power else NA_real_
  list(id = truth$id, driver = "low_f", effort = "standard", scheme = scheme,
    reference = truth$reference, data = data, holdout = holdout,
    sampling = do.call(rbind, sampling), truth = truth, audit = truth$audit,
    exact_settings = exact_settings, probabilities = do.call(rbind, probabilities),
    ablation = if (scheme == "poisson-link") delta_share_ablation_support(truth, calibration) else NULL,
    provenance = list(population_seed = truth$seed, monitoring_seed_base = seed_base,
      calibration_ids = calibration$ids, simsurvey_revision = calibration$simsurvey_revision,
      population_law = "Default SimSurvey Z/R/biology with exact 20+ carryover",
      reference_fishing_share = .10, q = 1, units = "thousand fish"))
}
