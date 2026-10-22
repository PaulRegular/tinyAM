#' Compute Mohn's rho for retrospective analyses
#'
#' @description
#' Calculates Mohn's rho (the average proportional retrospective bias)
#' from a long-format data frame containing terminal and peeled estimates
#' across retrospective runs.
#'
#' @param data A data frame with columns (e.g., data frames in `pop` output from [fit_retro()]):
#'
#' - **year**: integer or numeric assessment year of the estimate.
#' - **fold**: integer or numeric label for the retrospective peel
#'   (typically the terminal year of the dataset used for that run).
#' - **est**: numeric estimate (e.g., SSB, F, recruitment) on which
#'   Mohn's rho is calculated.
#' - **is_proj**: logical; `TRUE` for projection rows, `FALSE` otherwise.
#'
#' @details
#' Projection rows are excluded. The largest retained `fold` is the reference
#' assessment. Each earlier peel contributes its estimate at `year == fold`,
#' compared with the reference estimate for that same year:
#' \deqn{\rho=\frac{1}{K}\sum_{k=1}^K
#' \frac{\widehat X^{(k)}_{t_k}-\widehat X^{(ref)}_{t_k}}
#' {\widehat X^{(ref)}_{t_k}}.}
#' Positive values mean the peeled fits tend to estimate higher values than
#' the reference fit. Supply a single quantity and, for age-specific output,
#' one age at a time. `fold` must already identify terminal historical years;
#' no year shifting or relabeling is performed. Each matched peel has equal weight.
#' If failed fits were dropped, the reference is the latest retained fit.
#' Missing differences are omitted; zero reference estimates can yield infinite
#' or undefined differences.
#'
#' @return A single numeric proportional difference. Returns `NaN` (a missing
#'   numeric value) if there are no non-missing comparisons.
#'
#' @examples
#' df <- data.frame(
#'   year = rep(2010:2015, each = 3),
#'   fold = rep(c(2013, 2014, 2015), times = 6),
#'   est = c(
#'     100, 105, 110, 90, 95, 100, 85, 92, 98,
#'     80, 90, 95, 78, 88, 93, 76, 86, 91
#'   ),
#'   is_proj = FALSE
#' )
#' compute_mohns_rho(df)
#'
#' @references
#' Mohn, R. (1999). The retrospective problem in sequential population analysis:
#' An investigation using cod fishery and simulated data. *ICES Journal of
#' Marine Science*, 56(4), 473–488.
#'
#' @export
compute_mohns_rho <- function(data) {
  d <- data[!data$is_proj, ]
  td <- d[d$fold == max(d$fold), c("year", "est")] # terminal data
  rd <- d[
    d$fold != max(d$fold) &
      d$year == d$fold,
    c("year", "est", "fold")
  ] # retrospective data
  cd <- merge(rd, td, by = "year", suffixes = c("_r", "_t")) # combined data (r = retro, t = terminal)
  cd$pdiff <- (cd$est_r - cd$est_t) / cd$est_t
  mean(cd$pdiff, na.rm = TRUE)
}



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
  scores$scale <- ifelse(scores$metric == "mean_age", "years", "proportion")
  list(pairs = proportions, mean_age = mean_age, scores = scores)
}


#' Compute hindcast RMSE (observed vs projected)
#'
#' @description
#' Calculates the root–mean–squared error (RMSE) between observed values and
#' one–step–ahead projections from a hindcast run (e.g., a data frame from
#' `hindcasts$obs_pred$catch` or `hindcasts$obs_pred$index`).
#'
#' @param data A long-format data frame with columns (e.g., from `obs_pred` in
#'   [fit_hindcast()] output):
#'
#' - **year**: integer assessment year.
#' - **age**: integer model age.
#' - **obs**: numeric observed value for the given `year` × `age`.
#' - **pred**: numeric projected value for the given `year` × `age`.
#' - **fold**: integer or numeric label for the hindcast peel
#'   (typically the terminal year used in that run).
#' - **is_proj**: logical; `TRUE` for projection rows (the one–step–ahead
#'   predictions), `FALSE` otherwise.
#'
#' @param log Logical; if `TRUE` (default), compute RMSE on the log scale.
#'   Zeros in `obs`/`pred` are converted to `NA` before logging.
#'
#' @details
#' Matches projected rows (`is_proj == TRUE`) with available historical
#' observations by year and age, also using `type` and `survey` when present.
#' Include `type` when combining catch and index tables, and `survey` when
#' combining surveys. Conflicting observations for the same key cause an error.
#' Repeated copies of the same observation across folds count only once per
#' prediction. Each matched forecast receives equal weight.
#' Matched rows from `result$forecast_pairs$index` are also accepted directly.
#'
#' On the natural scale,
#' \deqn{\mathrm{RMSE}=\sqrt{\frac{1}{K}\sum_i(Y_i-\widetilde Y_i)^2}.}
#' With `log = TRUE`, both values are logged first. Log RMSE measures departures
#' on a relative scale and is often more useful when catch and survey units differ.
#' It still combines errors from different series and is not a likelihood score.
#' Zeros and missing pairs are omitted for log RMSE; zeros are retained on the
#' natural scale. Predictions beyond the available observations do not contribute.
#' All projection rows are eligible; use [fit_hindcast()] for one-year horizons
#' or subset the input to the desired forecast horizon.
#'
#' @return A single numeric RMSE. Returns `NaN` (a missing numeric value) if
#'   there are no non-missing matched forecast errors.
#'
#' @examples
#' \dontrun{
#' # Suppose `hindcasts <- fit_hindcast(fit, folds = 3)`
#' # Overall RMSE for the index series on the log scale
#' compute_hindcast_rmse(hindcasts$obs_pred$index, log = TRUE)
#'
#' # RMSE for catch on the natural scale
#' compute_hindcast_rmse(hindcasts$obs_pred$catch, log = FALSE)
#' }
#'
#' @export
compute_hindcast_rmse <- function(data, log = TRUE) {
  if (!"is_proj" %in% names(data) || "baseline" %in% names(data)) {
    if (!"baseline" %in% names(data)) data$baseline <- rep(NA_real_, nrow(data))
    return(.score_forecasts(data, log)$RMSE)
  }
  keys <- c("year", "age", intersect(c("model", "type", "survey", "fleet"), names(data)))
  proj_d <- data[data$is_proj, c(keys, "pred"), drop = FALSE]
  obs_d <- unique(data[!data$is_proj & !is.na(data$obs), c(keys, "obs"), drop = FALSE])
  if (anyDuplicated(obs_d[keys])) {
    cli::cli_abort("Hindcast scoring requires one observed value per year, age, and series. Separate distinct series with {.field type} or {.field survey}.")
  }
  d <- merge(obs_d, proj_d, by = keys)
  if (log) {
    d$obs[d$obs == 0] <- NA
    d$pred[d$pred == 0] <- NA
    d$obs <- log(d$obs)
    d$pred <- log(d$pred)
  }
  d$sq_error <- (d$obs - d$pred)^2
  sqrt(mean(d$sq_error, na.rm = TRUE))
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

