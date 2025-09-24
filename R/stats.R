
#' Compute Mohn's rho for retrospective analyses
#'
#' @description
#' Calculates Mohn's rho (the average proportional retrospective bias)
#' from a long-format data frame containing terminal and peeled estimates
#' across retrospective runs.
#'
#' @param data A data frame with columns:
#' \describe{
#'   \item{year}{Integer or numeric assessment year of the estimate.}
#'   \item{retro_year}{Integer or numeric label for the retrospective peel
#'     (typically the terminal year of the dataset used for that run).}
#'   \item{est}{Numeric estimate (e.g., SSB, F, recruitment) on which
#'     Mohn's rho is calculated.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Coerces `retro_year` to numeric.
#'   \item If `max(year) > max(retro_year)`, adds 1 to `retro_year` to handle
#'         conventions where each peel's estimate is reported for the
#'         subsequent calendar year (i.e., “stepped forward by one year”).
#'   \item Identifies terminal estimates (`retro_year == max(retro_year)`) as the
#'         reference series.
#'   \item For each peel, extracts the estimate at the peel year
#'         (`year == retro_year`) and merges with the terminal estimate for the
#'         same year.
#'   \item Computes proportional difference
#'         \eqn{(est_{retro} - est_{terminal}) / est_{terminal}} and returns
#'         the mean across peels.
#' }
#'
#' Observations with zero terminal estimates will produce `Inf`/`NaN`
#' in the proportional difference; consider filtering or transforming
#' input accordingly if this is a concern.
#'
#' @return A single numeric value: the mean proportional difference
#' (Mohn's rho). `NA` if no valid pairs are available.
#'
#' @examples
#' df <- data.frame(
#'   year = rep(2010:2015, each = 3),
#'   retro_year = rep(c(2013, 2014, 2015), times = 6),
#'   est = c(100, 105, 110, 90, 95, 100, 85, 92, 98, 80, 90, 95, 78, 88, 93, 76, 86, 91)
#' )
#' compute_mohns_rho(df)
#'
#' @references
#' Mohn, R. (1999). The retrospective problem in sequential population analysis:
#' An investigation using cod fishery and simulated data. \emph{ICES Journal of
#' Marine Science}, 56(4), 473–488.
#'
#' @export
compute_mohns_rho <- function(data) {

  d <- data
  d$retro_year <- as.numeric(as.character(d$retro_year))
  if (max(d$year) > max(d$retro_year)) {
    d$retro_year <- d$retro_year + 1 # plus one if the estimates are stepped forward by one year
  }
  td <- d[d$retro_year == max(d$retro_year), c("year", "est")]   # terminal data
  rd <- d[d$retro_year != max(d$retro_year) &
            d$year == d$retro_year,
          c("year", "est", "retro_year")]                         # retrospective data
  cd <- merge(rd, td, by = "year", suffixes = c("_r", "_t"))      # combined data (r = retro, t = terminal)
  cd$pdiff <- (cd$est_r - cd$est_t) / cd$est_t
  mean(cd$pdiff)

}

