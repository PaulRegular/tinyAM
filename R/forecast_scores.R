# Scoring uses original held-out observations, never a subsequent fitted state.
.bind_rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), rows)
  if (!length(rows)) return(data.frame())
  keys <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    for (key in setdiff(keys, names(x))) x[[key]] <- NA
    x[keys]
  })
  rownames_out <- do.call(rbind, rows)
  rownames(rownames_out) <- NULL
  rownames_out
}

.forecast_pairs <- function(predictions, observations) {
  keys <- c("year", "age", intersect(c("survey", "fleet"), names(observations)))
  observations <- observations[!observations$is_proj, , drop = FALSE]
  obs <- unique(observations[c(keys, "obs", intersect("unit", names(observations)))])
  if (anyDuplicated(obs[keys])) cli::cli_abort("Forecast scoring requires unique observations by year, age and series.")
  p <- predictions[predictions$is_proj, , drop = FALSE]
  p$obs <- p$unit <- NULL
  if (!nrow(p)) return(data.frame())
  d <- merge(p, obs, by = keys, all.x = TRUE, sort = FALSE)
  .add_persistence(d, obs, setdiff(keys, "year"), positive = "survey" %in% keys)
}

.add_persistence <- function(pairs, history, keys, positive = FALSE) {
  pairs$baseline <- NA_real_
  pairs$baseline_year <- NA_real_
  for (i in seq_len(nrow(pairs))) {
    keep <- is.finite(history$obs) & history$year <= pairs$fold[i]
    if (positive) keep <- keep & history$obs > 0
    for (key in keys) keep <- keep & as.character(history[[key]]) == as.character(pairs[[key]][i])
    d <- history[which(keep), , drop = FALSE]
    if (nrow(d)) {
      j <- which.max(d$year)
      pairs$baseline[i] <- d$obs[j]
      pairs$baseline_year[i] <- d$year[j]
    }
  }
  pairs$baseline_gap <- pairs$year - pairs$baseline_year
  pairs
}

.score_forecasts <- function(d, log = TRUE) {
  valid <- is.finite(d$obs) & is.finite(d$pred)
  if (log) valid <- valid & d$obs > 0 & d$pred > 0
  baseline_valid <- valid & is.finite(d$baseline)
  if (log) baseline_valid <- baseline_valid & d$baseline > 0
  transform <- if (log) base::log else identity
  error <- transform(d$obs[valid]) - transform(d$pred[valid])
  e <- abs(transform(d$obs[baseline_valid]) - transform(d$pred[baseline_valid]))
  b <- abs(transform(d$obs[baseline_valid]) - transform(d$baseline[baseline_valid]))
  denominator <- sum(b)
  data.frame(RMSE = if (length(error)) sqrt(mean(error^2)) else NA_real_,
    MAE = if (length(error)) mean(abs(error)) else NA_real_,
    MASE = if (length(b) && denominator > 0) sum(e) / denominator else NA_real_,
    n = sum(valid), n_mase = sum(baseline_valid),
    mase_reason = if (!length(b)) "No valid forecast/persistence pairs" else if (denominator == 0) "Zero persistence error" else NA_character_)
}

#' Persistence-relative hindcast forecast error
#'
#' Compare forecasts with holding the last observation available at the forecast
#' origin constant. Values below one indicate smaller absolute errors than
#' persistence. This uses the evaluation-period convention of Carvalho et al.,
#' rather than a training-period scaling denominator.
#' @param data Matched forecast rows containing `obs`, `pred` and `baseline`.
#'   The baseline must use observations available at the forecast origin.
#' @param log Use log-scale errors (default); observations and baselines must
#'   then be positive. Use `FALSE` for proportions or mean modeled age.
#' @return A number, or `NA` when there are no valid pairs or persistence has
#'   zero error. See the score tables from [fit_hindcast()] for counts and reasons.
#' @details Numerator and denominator use the same forecast cases:
#'   \deqn{\mathrm{MASE}=\sum_i|z_i-\widehat z_i|/\sum_i|z_i-z_{0,i}|.}
#'   For multiple horizons the baseline remains fixed at the origin. This is
#'   an extension of the paper's one-year calculation. Predictions are conditional
#'   medians, as in existing log-RMSE calculations; this is not a probability score.
#' @references Carvalho et al. (2021). A cookbook for using model diagnostics in
#'   integrated stock assessments. Fisheries Research 240, 105959.
#'   \doi{10.1016/j.fishres.2021.105959}.
#' @examples
#' compute_hindcast_mase(data.frame(obs = c(2, 4), pred = c(2.1, 3.8),
#'                                 baseline = c(1, 2)))
#' @export
compute_hindcast_mase <- function(data, log = TRUE) {
  if (!all(c("obs", "pred", "baseline") %in% names(data))) {
    cli::cli_abort("MASE requires matched {.field obs}, {.field pred} and {.field baseline} columns.")
  }
  .score_forecasts(data, log)$MASE
}

.aggregate_indices <- function(data, supports, predicted = FALSE) {
  if (!nrow(data)) return(data.frame())
  groups <- intersect(c("model", "fold", "horizon", "year", "survey"), names(data))
  .bind_rows(lapply(split(data, interaction(data[groups], drop = TRUE)), function(d) {
    ages <- supports[[as.character(d$survey[1])]]
    if (anyDuplicated(d$age) || !setequal(d$age, ages)) return(NULL)
    values <- d$obs
    if (any(!is.finite(values) | values <= 0)) return(NULL)
    out <- d[1, groups, drop = FALSE]
    if ("unit" %in% names(d)) out$unit <- d$unit[1]
    out$obs <- sum(values)
    if (predicted) out$pred <- if (all(is.finite(d$pred) & d$pred > 0)) sum(d$pred) else NA_real_
    out
  }))
}

.index_score_tables <- function(pairs, observations) {
  if (!nrow(pairs)) return(list(pairs = pairs, aggregate = data.frame(), scores = data.frame()))
  supports <- lapply(split(observations$age, observations$survey), function(x) sort(unique(x)))
  aggregate <- .aggregate_indices(pairs, supports, TRUE)
  history <- .aggregate_indices(observations[!observations$is_proj, ], supports)
  if (nrow(aggregate)) aggregate <- .add_persistence(aggregate, history, "survey")
  tables <- list()
  for (level in c("age", "aggregate")) {
    d <- if (level == "age") pairs else aggregate
    if (!nrow(d)) next
    if (level == "aggregate") d$age <- NA_real_
    for (h in c(0, sort(unique(d$horizon)))) {
      x <- if (h == 0) d else d[d$horizon == h, ]
      groups <- if (level == "age") c("survey", "age") else "survey"
      tables <- c(tables, lapply(split(x, interaction(x[groups], drop = TRUE)), function(z) {
        cbind(data.frame(level = level, survey = as.character(z$survey[1]), age = z$age[1],
          horizon = if (h == 0) NA_real_ else h, scale = "log"), .score_forecasts(z))
      }))
    }
  }
  scores <- .bind_rows(tables)
  for (h in c(0, sort(unique(pairs$horizon)))) {
    x <- if (h == 0) pairs else pairs[pairs$horizon == h, ]
    s <- scores[scores$level == "age" & !is.na(scores$horizon), ]
    if (h != 0) s <- s[s$horizon == h, ]
    # Equal weights for ages within survey, surveys within horizon, then horizons.
    means <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
    survey <- if (nrow(s)) aggregate(s$MASE, s[c("survey", "horizon")], means) else data.frame()
    hs <- if (nrow(survey)) vapply(split(survey$x, survey$horizon), means, numeric(1)) else numeric()
    score <- .score_forecasts(x)
    score$MASE <- means(hs)
    score$mase_reason <- if (is.finite(score$MASE)) NA_character_ else "No eligible series"
    scores <- rbind(scores, cbind(data.frame(level = "overall", survey = NA_character_, age = NA_real_,
      horizon = if (h == 0) NA_real_ else h, scale = "log"), score))
  }
  list(pairs = pairs, aggregate = aggregate, scores = scores)
}

.catch_compositions <- function(pairs, observations, ages) {
  if (!nrow(pairs)) return(list(pairs = data.frame(), mean_age = data.frame(), scores = data.frame()))
  groups <- intersect(c("model", "fold", "horizon", "year", "fleet"), names(pairs))
  proportions <- .bind_rows(lapply(split(pairs, interaction(pairs[groups], drop = TRUE)), function(d) {
    if (anyDuplicated(d$age) || !setequal(d$age, ages) ||
      any(!is.finite(d$obs) | d$obs < 0 | !is.finite(d$pred) | d$pred < 0) ||
      sum(d$obs) <= 0 || sum(d$pred) <= 0) return(NULL)
    d$obs <- d$obs / sum(d$obs); d$pred <- d$pred / sum(d$pred)
    d$unit <- "proportion"
    d
  }))
  if (!nrow(proportions)) return(list(pairs = proportions, mean_age = data.frame(), scores = data.frame()))
  history_groups <- c("year", intersect("fleet", names(observations)))
  history <- .bind_rows(lapply(split(observations[!observations$is_proj, ],
    interaction(observations[!observations$is_proj, history_groups, drop = FALSE], drop = TRUE)), function(d) {
    if (!setequal(d$age, ages) || anyDuplicated(d$age) || any(!is.finite(d$obs) | d$obs < 0) || sum(d$obs) <= 0) return(NULL)
    d$obs <- d$obs / sum(d$obs); d
  }))
  proportions <- .add_persistence(proportions, history, c("age", intersect("fleet", names(history))))
  means <- function(d, groups, prediction) .bind_rows(lapply(split(d, interaction(d[groups], drop = TRUE)), function(z) {
    x <- z[1, groups, drop = FALSE]; x$obs <- sum(z$age * z$obs)
    if (prediction) x$pred <- sum(z$age * z$pred)
    x
  }))
  mean_age <- means(proportions, groups, TRUE)
  mean_age$unit <- "years"
  mean_history <- if (nrow(history)) means(history, history_groups, FALSE) else data.frame()
  mean_age <- .add_persistence(mean_age, mean_history, intersect("fleet", names(mean_history)))
  scores <- list()
  for (h in c(0, sort(unique(proportions$horizon)))) {
    p <- if (h == 0) proportions else proportions[proportions$horizon == h, ]
    m <- if (h == 0) mean_age else mean_age[mean_age$horizon == h, ]
    scores <- c(scores, list(cbind(data.frame(metric = "composition", age = NA_real_, horizon = if (h == 0) NA_real_ else h), .score_forecasts(p, FALSE)),
      cbind(data.frame(metric = "mean_age", age = NA_real_, horizon = if (h == 0) NA_real_ else h), .score_forecasts(m, FALSE))))
    scores <- c(scores, lapply(split(p, p$age), function(z) cbind(data.frame(metric = "proportion", age = z$age[1], horizon = if (h == 0) NA_real_ else h), .score_forecasts(z, FALSE))))
  }
  scores <- .bind_rows(scores)
  if ("fleet" %in% names(pairs)) {
    scores$fleet <- NA_character_
    by_fleet <- lapply(unique(as.character(pairs$fleet)), function(fleet) {
      p <- pairs[as.character(pairs$fleet) == fleet, , drop = FALSE]
      o <- observations[as.character(observations$fleet) == fleet, , drop = FALSE]
      p$fleet <- o$fleet <- NULL
      d <- .catch_compositions(p, o, ages)$scores
      if (nrow(d)) d$fleet <- fleet
      d
    })
    scores <- .bind_rows(c(list(scores), by_fleet))
  }
  scores$scale <- ifelse(scores$metric == "mean_age", "years", "proportion")
  list(pairs = proportions, mean_age = mean_age, scores = scores)
}
