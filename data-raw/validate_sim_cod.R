# Rscript data-raw/validate_sim_cod.R
# Deterministic checks of the monitoring expansion; no SimSurvey installation,
# population simulation or assessment fit is required.

validate_sim_cod_monitoring <- function(repo = ".") {
  code <- new.env(parent = globalenv())
  sys.source(file.path(repo, "data-raw", "sim_cod_obs.R"), envir = code)
  raise <- code$raise_sim_cod_catch
  if (!is.function(raise)) stop("The generator must define raise_sim_cod_catch().")

  operations <- data.frame(
    fishery = "Gillnet", sim = 1L, year = 1L, operation = 1:3,
    reported_n = c(100, 400, 300), sampling_probability = c(.5, 1, .5),
    sampled = c(TRUE, TRUE, FALSE)
  )
  samples <- data.frame(
    fishery = "Gillnet", sim = 1L, year = 1L,
    operation = rep(1:2, each = 4),
    length = c(30, 30, 30, 60, 30, 60, 60, 60),
    age = c(2, NA, NA, 3, 2, 3, NA, NA),
    aged = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    measured = TRUE
  )
  original <- list(operations = operations, samples = samples)
  evaluate <- function(op = operations, sa = samples) {
    raise(op, sa, ages = 1:4, years = 1, fleets = "Gillnet",
      age_length_group = 3)
  }
  equal <- function(actual, expected, label) {
    if (!isTRUE(all.equal(actual, expected, check.attributes = FALSE,
                         tolerance = 1e-12))) stop(label)
  }

  result <- evaluate()
  catch <- result$catch[order(result$catch$age), ]
  # Operation 1 expands to (75,25), then its 1/.5 operation weight
  # gives (150,50). Operation 2 expands to (100,300). H=(250,350).
  # Calibrating these to all three reported totals (800) gives:
  equal(catch$obs, c(0, 1000 / 3, 1400 / 3, 0),
    "Conditional length-age expansion or operation weights are incorrect.")
  equal(sum(catch$obs), 800, "Estimated catch must equal reported annual catch.")
  stopifnot(identical(original, list(operations = operations, samples = samples)))

  shuffled <- evaluate(operations[c(3, 1, 2), ], samples[c(8, 2, 5, 1, 4, 7, 3, 6), ])
  equal(shuffled$catch$obs, result$catch$obs,
    "Monitoring expansion depends on input row order.")
  scaled <- operations
  scaled$reported_n <- scaled$reported_n * 10
  equal(evaluate(scaled)$catch$obs, 10 * result$catch$obs,
    "Reported catch units are not propagated linearly.")

  # True operation catches and unobserved ages are unavailable to the estimator.
  hidden_operations <- operations
  hidden_operations$n <- c(1e8, 2e8, 3e8)
  hidden_samples <- samples
  hidden_samples$true_age <- rep(99, nrow(hidden_samples))
  equal(evaluate(hidden_operations, hidden_samples)$catch$obs, result$catch$obs,
    "The catch estimator must not use hidden catch or age truth.")

  unaged_bin <- samples
  unaged_bin$age[unaged_bin$operation == 1 & unaged_bin$length == 60] <- NA
  unaged_bin$aged[unaged_bin$operation == 1 & unaged_bin$length == 60] <- FALSE
  stopifnot(all(is.na(evaluate(sa = unaged_bin)$catch$obs)))
  stopifnot(all(is.na(evaluate(sa = samples[FALSE, ])$catch$obs)))
  zero <- operations
  zero$reported_n <- 0
  equal(evaluate(zero, samples[FALSE, ])$catch$obs, rep(0, 4),
    "A genuinely zero reported annual catch should remain zero.")
  cat("Simulated-cod monitoring: deterministic expansion checks passed.\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) validate_sim_cod_monitoring()
