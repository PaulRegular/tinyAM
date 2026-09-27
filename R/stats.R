
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
#'   est = c(100, 105, 110, 90, 95, 100, 85, 92, 98,
#'           80, 90, 95, 78, 88, 93, 76, 86, 91),
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
  td <- d[d$fold == max(d$fold), c("year", "est")]   # terminal data
  rd <- d[d$fold != max(d$fold) &
            d$year == d$fold,
          c("year", "est", "fold")]                         # retrospective data
  cd <- merge(rd, td, by = "year", suffixes = c("_r", "_t"))      # combined data (r = retro, t = terminal)
  cd$pdiff <- (cd$est_r - cd$est_t) / cd$est_t
  mean(cd$pdiff, na.rm = TRUE)
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
  keys <- c("year", "age", intersect(c("type", "survey"), names(data)))
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
  d$sq_error <- (d$obs - d$pred) ^ 2
  sqrt(mean(d$sq_error, na.rm = TRUE))
}


