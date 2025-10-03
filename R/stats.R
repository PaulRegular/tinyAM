#' Compute Mohn's rho for retrospective analyses
#'
#' @description
#' Calculates Mohn's rho (the average proportional retrospective bias)
#' from a long-format data frame containing terminal and peeled estimates
#' across retrospective runs.
#'
#' @param data A data frame with columns:
#'
#' - **year**: integer or numeric assessment year of the estimate.
#' - **fold**: integer or numeric label for the retrospective peel
#'   (typically the terminal year of the dataset used for that run).
#' - **est**: numeric estimate (e.g., SSB, F, recruitment) on which
#'   Mohn's rho is calculated.
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

