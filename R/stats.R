
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
#' The function:
#'
#' 1. Coerces `fold` to numeric.
#' 2. If `max(year) > max(fold)`, adds 1 to `fold` to handle
#'    conventions where each peel's estimate is reported for the
#'    subsequent calendar year (i.e., “stepped forward by one year”).
#' 3. Identifies terminal estimates (`fold == max(fold)`) as the
#'    reference series.
#' 4. For each peel, extracts the estimate at the peel year
#'    (`year == fold`) and merges with the terminal estimate for the
#'    same year.
#' 5. Computes proportional difference
#'    \deqn{\rho = \frac{ \text{est}_{\text{retro}} - \text{est}_{\text{terminal}} }
#'                { \text{est}_{\text{terminal}} }}
#'    and returns the mean across peels.
#'
#' Observations with zero terminal estimates will produce `Inf`/`NaN`
#' in the proportional difference; consider filtering or transforming
#' input accordingly if this is a concern.
#'
#' @return
#' A single numeric value: the mean proportional difference (Mohn's rho).
#' Returns `NA` if no valid pairs are available.
#'
#' @examples
#' df <- data.frame(
#'   year = rep(2010:2015, each = 3),
#'   fold = rep(c(2013, 2014, 2015), times = 6),
#'   est = c(100, 105, 110, 90, 95, 100, 85, 92, 98,
#'           80, 90, 95, 78, 88, 93, 76, 86, 91)
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


  ## YOU ARE HERE FIXING UP THIS FUNCTION AND THE HINDCAST_RMSE FUNCTION TO MAKE SURE IT HANDELS IS_PROJ CORRECTLY


  d <- data
  d$fold <- as.numeric(as.character(d$fold))
  if (max(d$year) > max(d$fold)) {
    d$fold <- d$fold + 1 # plus one if the estimates are stepped forward by one year
  }
  td <- d[d$fold == max(d$fold), c("year", "est")]   # terminal data
  rd <- d[d$fold != max(d$fold) &
            d$year == d$fold,
          c("year", "est", "fold")]                         # retrospective data
  cd <- merge(rd, td, by = "year", suffixes = c("_r", "_t"))      # combined data (r = retro, t = terminal)
  cd$pdiff <- (cd$est_r - cd$est_t) / cd$est_t
  mean(cd$pdiff)

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
#' The function:
#'
#' 1. Selects one–step–ahead projections via `is_proj == TRUE` (projected rows).
#' 2. Extracts observed values at the peel year via `year == fold`
#'    (observed rows).
#' 3. Merges observed and projected values by `(year, age)`.
#' 4. Optionally transforms to log scale (after replacing exact zeros with `NA`).
#' 5. Computes squared errors \eqn{(\mathrm{obs} - \mathrm{pred})^2} and returns
#'    \deqn{\mathrm{RMSE} = \sqrt{\mathrm{mean}\big[(\mathrm{obs} - \mathrm{pred})^2\big]}}
#'    with `na.rm = TRUE`.
#'
#' Notes:
#' - If `log = TRUE`, exact zeros are dropped (set to `NA`) before `log()`.
#'   Consider adding a small constant beforehand if you prefer to retain zeros.
#' - Only year–age pairs present in **both** the observed-at-fold and projected
#'   sets contribute to the RMSE (via the merge).
#'
#' @return
#' A single numeric value: the RMSE between observed and one–step–ahead projected
#' values for matched `(year, age)` pairs. Returns `NA` if no valid pairs exist.
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
  proj_d <- data[data$is_proj, c("year", "age", "pred")]
  obs_d <- data[data$year == data$fold, c("year", "age", "obs")]
  d <- merge(obs_d, proj_d, by = c("year", "age"))
  if (log) {
    d$obs[d$obs == 0] <- NA
    d$pred[d$pred == 0] <- NA
    d$obs <- log(d$obs)
    d$pred <- log(d$pred)
  }
  d$sq_error <- (d$obs - d$pred) ^ 2
  sqrt(mean(d$sq_error, na.rm = TRUE))
}


