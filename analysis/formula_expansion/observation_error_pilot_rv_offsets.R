# Independent seed 21 survey calibration, restricted to the pilot design.
oe_pilot_cadigan_calibration <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  source <- readRDS(path)
  manifest <- source$manifest
  stopifnot(identical(manifest$phase, "study"), manifest$population_seed == 21L,
    manifest$grid_n == 40L, identical(manifest$simsurvey_revision, .observation_error_revision),
    identical(manifest$field_seeds, c(17L, 438L, 72L)))

  offsets <- as.data.frame(source$diagnostics$cadigan_offsets)
  d <- subset(offsets, component == "RV" & om == "default" & bin == "modeled" &
    effort == "standard")
  fields <- unique(d[c("population_id", "theta")])
  expected <- paste0("seed", manifest$field_seeds)
  n_primary <- unique(d$n_primary)
  n_equivalent <- unique(d$n_equivalent)
  stopifnot(nrow(fields) == length(expected), setequal(fields$population_id, expected),
    all(is.finite(fields$theta) & fields$theta >= 0), identical(as.numeric(n_primary), 160),
    length(n_equivalent) == 1L, is.finite(n_equivalent), n_equivalent > 0)
  fields <- fields[match(expected, fields$population_id), , drop = FALSE]

  list(theta = stats::median(fields$theta), field_theta = fields, source = path,
    md5 = unname(tools::md5sum(path)), grid = manifest$grid_n,
    population_seed = manifest$population_seed, years = manifest$years,
    field_seeds = manifest$field_seeds, calibration_draws = manifest$calibration,
    n_primary = n_primary, n_equivalent = n_equivalent)
}

oe_pilot_rv_offset <- function(training_rv, model, baseline_index = NULL, calibration = NULL) {
  model <- match.arg(model, c("rv_design_annual", "rv_cadigan"))
  keys <- c("year", "age")
  stopifnot(all(c(keys, "obs", "se") %in% names(training_rv)), nrow(training_rv) > 0L,
    !anyNA(training_rv[keys]), !anyDuplicated(training_rv[keys]))
  out <- training_rv[keys]

  if (model == "rv_design_annual") {
    positive <- is.finite(training_rv$obs) & training_rv$obs > 0 &
      is.finite(training_rv$se) & training_rv$se > 0
    out$rv_log_sd <- rep(NA_real_, nrow(out))
    out$rv_log_sd[positive] <- oe_diag_reported_log_sd(training_rv$obs[positive], training_rv$se[positive])
    out$offset_source <- "annual"
    needs_fallback <- !is.finite(out$rv_log_sd) | out$rv_log_sd <= 0
    for (year in unique(out$year[needs_fallback])) {
      rows <- out$year == year & needs_fallback
      pooled <- oe_pooled_rv_sd(training_rv[training_rv$year != year, , drop = FALSE],
        out$age[rows])
      out$rv_log_sd[rows] <- pooled
      out$offset_source[rows] <- ifelse(attr(pooled, "fallback"), "pooled_median", "pooled_age")
    }
  } else {
    stopifnot(all(c(keys, "pred", "sd") %in% names(baseline_index)),
      !anyNA(baseline_index[keys]), !anyDuplicated(baseline_index[keys]),
      setequal(.oe_fit_key(baseline_index), .oe_fit_key(training_rv)),
      length(calibration$theta) == 1L, is.finite(calibration$theta), calibration$theta >= 0,
      length(training_rv$n_primary) == nrow(out), length(training_rv$n_equivalent) == nrow(out),
      all(training_rv$n_primary == calibration$n_primary),
      isTRUE(all.equal(training_rv$n_equivalent,
        rep(calibration$n_equivalent, nrow(out)), tolerance = 1e-10)))
    rows <- match(.oe_fit_key(training_rv), .oe_fit_key(baseline_index))
    median_index <- baseline_index$pred[rows]
    fitted_sd <- baseline_index$sd[rows]
    stopifnot(all(is.finite(median_index) & median_index > 0),
      all(is.finite(fitted_sd) & fitted_sd > 0))
    # Convert the baseline's lognormal median in thousand fish to the
    # arithmetic mean in fish per tow, then freeze it for the SD fit.
    log_mean_tow <- log(median_index) + fitted_sd^2 / 2 + log(1000) -
      log(training_rv$n_equivalent)
    out$rv_log_sd <- sqrt(log1p((exp(-log_mean_tow) + calibration$theta) / training_rv$n_primary))
    out$offset_source <- "cadigan"
  }

  stopifnot(all(is.finite(out$rv_log_sd) & out$rv_log_sd > 0))
  out
}

oe_pilot_apply_rv_offset <- function(input, settings, offsets) {
  keys <- .oe_fit_key(offsets)
  index_keys <- .oe_fit_key(input$index)
  stopifnot(!anyDuplicated(keys), !anyDuplicated(index_keys), setequal(keys, index_keys))
  rows <- match(index_keys, keys)
  input$index$rv_log_sd <- offsets$rv_log_sd[rows]
  input$index$rv_sd_scale <- 1
  settings$index_settings$sd_supplied <- ~ rv_log_sd
  settings$index_settings$sd_form <- ~ 0 + rv_sd_scale
  obs <- input$index$obs
  fallback <- sum(is.finite(obs) & obs > 0 & startsWith(offsets$offset_source[rows], "pooled_"))
  list(input = input, settings = settings, fallback = fallback)
}
