#' Double-logistic selectivity curves
#'
#' Use this additive formula marker in `F_settings$mu_form` or
#' `index_settings$q_form`, including inside [for_fleet()] and [for_survey()].
#' Ordinary formula coefficients provide the F or q scale; the curve multiplies
#' it. This curve does not select an age at which F must equal its scale.
#'
#' @param x Name of a numeric age/size column. F curves require age or a catch
#'   covariate constant over years at each age within its applicable fleet.
#' @param by Optional categorical column giving separate curves.
#' @details
#' The four-parameter curve follows WHAM's raw product:
#' \deqn{S(x)=\frac{1}{1+\exp\{-(x-a_{50,up})/k_{up}\}}
#'             \frac{1}{1+\exp\{(x-a_{50,down})/k_{down}\}},}
#' with positive \eqn{k_{up},k_{down}}. There is no peak normalization or
#' midpoint ordering. The accompanying F/q coefficient is a scale multiplier,
#' not necessarily maximum realized F/q. Reversed midpoints are permitted;
#' sparse tails and alternative parameter combinations can yield weak or
#' multimodal estimates. Check uncertainty and recovery before routine use.
#'
#' WHAM's age counter runs from 1 to the number of ages; tinyAM instead uses
#' the numeric coordinates in `x`. Account for this offset in comparisons.
#' F curves support no residual process, IID or AR1 residuals; static curves
#' cancel from intrinsic RW/correlated-RW increments and are rejected there.
#' Catchability curves multiply q after its inverse link, including logit q.
#'
#' @references Miller, T. J. and Stock, B. C., WHAM, `src/selectivity.hpp`,
#'   `get_selAA()`, double-logistic model 3, revision
#'   `ffcbfd093161e8210aa20716933975191b661807` (GPL-3).
#'   <https://github.com/timjmiller/wham/blob/ffcbfd093161e8210aa20716933975191b661807/src/selectivity.hpp>
#' @return Formula marker; direct calls raise an error.
#' @examples
#' ~ double_logistic(age)
#' ~ double_logistic(age, by = fleet)
#' ~ 0 + for_fleet("Trawl", ~ factor(year) + logistic(age)) +
#'   for_fleet("Gillnet", ~ factor(year) + double_logistic(age))
#' @seealso [catchability_curves], [for_fleet()], [for_survey()], [fit_tam()]
#' @export
double_logistic <- function(x, by = NULL) {
  cli::cli_abort("double_logistic() is a formula marker for F_settings$mu_form or index_settings$q_form.")
}

.is_selectivity_curve <- function(term) term$type %in% c("logistic", "double_logistic")

# Parameter names keep q/F roles explicit, including inside scoped formulas.
.curve_parameter_names <- function(term, component = "q") {
  id <- term$id
  fields <- switch(term$type,
    logistic = c(a50 = "a50", slope = "log_slope"),
    double_logistic = c(a50_up = "a50_up", scale_up = "log_scale_up",
      a50_down = "a50_down", scale_down = "log_scale_down"))
  for (field in names(fields)) {
    type <- fields[[field]]
    prefix <- if (startsWith(type, "log_")) paste0("log_", component, "_", substring(type, 5L)) else
      paste0(component, "_", type)
    term[[paste0(field, "_parameter")]] <- paste0(prefix, "_", id)
  }
  term
}

.curve_parameters <- function(term) {
  labels <- vapply(term$groups, `[[`, character(1), "group")
  mid <- vapply(term$groups, `[[`, numeric(1), "midpoint")
  width <- vapply(term$groups, `[[`, numeric(1), "width")
  values <- switch(term$type,
    logistic = list(a50 = mid, slope = log(4 / width)),
    double_logistic = list(a50_up = mid - width / 4, scale_up = log(width / 8),
      a50_down = mid + width / 4, scale_down = log(width / 8)))
  stats::setNames(lapply(values, stats::setNames, nm = labels),
    vapply(names(values), function(x) term[[paste0(x, "_parameter")]], character(1)))
}

.curve_log_selectivity <- function(x, par, term, group) {
  get <- function(field) par[[term[[paste0(field, "_parameter")]]]][group]
  if (term$type == "logistic") {
    predictor <- exp(get("slope")) * (x - get("a50"))
    return(-RTMB::logspace_add(0, -predictor))
  }
  # Adapted from WHAM get_selAA(), pinned in the help above. GPL-3.
  up <- (x - get("a50_up")) / exp(get("scale_up"))
  down <- (x - get("a50_down")) / exp(get("scale_down"))
  -RTMB::logspace_add(0, -up) - RTMB::logspace_add(0, down)
}

.curve_jacobian <- function(term, data) {
  parameters <- .curve_parameters(term)
  lengths <- lengths(parameters)
  offsets <- cumsum(c(0L, lengths))
  value <- unlist(parameters, use.names = FALSE)
  tape <- RTMB::MakeTape(function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    parameters <- stats::setNames(lapply(seq_along(lengths), function(i)
      p[seq.int(offsets[i] + 1L, offsets[i + 1L])]), names(parameters))
    out <- numeric(nrow(data))
    for (i in seq_along(term$groups)) {
      rows <- term$groups[[i]]$rows
      out[rows] <- .curve_log_selectivity(data[[term$variable]][rows], parameters, term, i)
    }
    out
  }, value)
  tape$jacobian(value)
}
