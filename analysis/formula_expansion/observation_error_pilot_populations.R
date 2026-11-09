# The pilot uses the demographic settings from oe_population() with the
# default spatial field only. Source observation_error_sim.R first.
.oe_pilot_history <- function(population_seed) {
  .observation_error_with_seed(population_seed, SimSurvey::sim_abundance(
    Z = function(ages, years) outer(0.3 + 0.06 * .oe_gillnet_q(ages), rep(1, length(years))),
    maturity = SimSurvey::sim_maturity(k = 1, x50 = 5), maturity_by = "age"))
}

oe_pilot_population <- function(population_seed, field_seed, grid_n) {
  .observation_error_check_revision()
  stopifnot(grid_n == 40L)
  history <- .oe_pilot_history(population_seed)
  grid <- SimSurvey::make_grid(res = rep(280 / grid_n, 2))
  pop <- .observation_error_with_seed(field_seed,
    SimSurvey::sim_distribution(history, grid = grid, ays_covar = SimSurvey::sim_ays_covar()))
  totals <- round(history$N)
  pop$sp_N[, N := .oe_balance(N, totals[.BY$age, .BY$year]), by = .(year, age)]
  pop <- SimSurvey::round_sim(pop)
  stopifnot(isTRUE(all.equal(pop$N, totals, tolerance = 0)))
  pop
}

.oe_pilot_regime <- function(ssb) {
  stopifnot(length(ssb) == 20L, all(is.finite(ssb) & ssb > 0))
  smoothed <- as.numeric(stats::filter(ssb, rep(1 / 3, 3)))
  initial <- smoothed[2L]
  middle <- smoothed[10L]
  terminal <- smoothed[19L]
  peak_year <- which.max(smoothed)
  peak <- smoothed[peak_year]

  if (terminal >= 0.8 * initial && terminal <= 1.2 * initial &&
      peak <= 1.5 * min(smoothed, na.rm = TRUE)) {
    return("stable")
  }
  if (peak_year %in% 5:16 && peak >= 1.3 * initial && terminal <= 0.75 * peak) {
    return("rise_fall")
  }
  if (terminal <= 0.65 * initial && terminal < middle && middle < initial &&
      peak < 1.3 * initial) {
    return("decline")
  }
  NA_character_
}

oe_pilot_population_manifest <- function() {
  .observation_error_check_revision()
  regimes <- c("stable", "decline", "rise_fall")
  selected <- stats::setNames(vector("list", length(regimes)), regimes)
  for (seed in 10001:11000) {
    history <- .oe_pilot_history(seed)
    ssb <- colSums(history$S) / 1000
    regime <- .oe_pilot_regime(ssb)
    if (!is.na(regime) && length(selected[[regime]]) < 8L) {
      selected[[regime]][[length(selected[[regime]]) + 1L]] <- list(seed = seed, ssb = ssb)
    }
    if (all(lengths(selected) == 8L)) break
  }
  stopifnot(all(lengths(selected) == 8L))

  # Interleave trajectory classes. Repeat monitoring for stocks two and six
  # in each class, assigned before monitoring and fitting.
  manifest <- expand.grid(regime = regimes, within_regime = seq_len(8L),
    stringsAsFactors = FALSE)
  manifest$population_id <- sprintf("pop%02d", seq_len(nrow(manifest)))
  manifest$repeat_monitoring <- manifest$within_regime %in% c(2L, 6L)
  stocks <- Map(function(regime, i) selected[[regime]][[i]],
    manifest$regime, manifest$within_regime)
  manifest$population_seed <- vapply(stocks, `[[`, integer(1), "seed")

  # These stream bases are spaced far enough apart to add year - 1 to each
  # removals or monitoring seed without sharing seeds across RNG roles.
  roles <- c("field", "removals", "gillnet_fit", "rv_fit", "gillnet_holdout",
    "rv_holdout", "gillnet_repeat_fit", "rv_repeat_fit", "gillnet_repeat_holdout",
    "rv_repeat_holdout")
  seed_base <- 900000000L + 1000L * seq_len(nrow(manifest))
  for (i in seq_along(roles)) {
    manifest[[paste0(roles[i], "_seed")]] <- seed_base + 50L * i
  }

  trajectories <- do.call(rbind, lapply(seq_along(stocks), function(i) {
    data.frame(population_id = manifest$population_id[i], regime = manifest$regime[i],
      year = seq_along(stocks[[i]]$ssb), ssb = as.numeric(stocks[[i]]$ssb))
  }))
  list(populations = manifest, ssb = trajectories, last_screened_seed = seed)
}
