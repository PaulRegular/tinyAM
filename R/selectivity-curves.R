#' Asymptotic and dome-shaped selectivity curves
#'
#' Use these additive formula markers in `F_settings$mu_form` or
#' `index_settings$q_form`, including inside [for_fleet()] and [for_survey()].
#' Ordinary formula coefficients provide the F or q scale; the curve multiplies
#' it. These curves do not select an age at which F must equal its scale.
#'
#' @param x Name of a numeric age/size column. F curves require age or a catch
#'   covariate constant over years at each age within its applicable fleet.
#' @param by Optional categorical column giving separate curves.
#' @param b Fixed nonnegative shape, including `Inf`, or `NULL` to estimate it.
#' @details
#' ## Warped normal
#'
#' With \eqn{u=(x-a_{50})/s}, positive \eqn{s}, and \eqn{b\geq0},
#' \deqn{g_b(u)=(1+b)\exp\{-u/(1+b)\}-b,\qquad S(x)=2^{-g_b(u)^2}.}
#' The ascending half-selection coordinate is \eqn{a_{50}} and the slope there
#' is \eqn{\log(2)/s}. At `b = 0` the curve rises towards one. Finite positive
#' b gives a dome with upper tail \eqn{2^{-b^2}}. At `b = Inf`,
#' \eqn{g(u)=1-u}: the curve is Gaussian, centered at \eqn{a_{50}+s}.
#' The peak, asymmetry and tail are coupled; this is not a double-normal curve.
#'
#' The fitted shape is \eqn{c=1/(1+b)} with bounds zero and one. A series around
#' \eqn{cu=0} preserves the Gaussian limit and its first and second derivatives.
#' [tidy_par()] includes derived peak and tail estimates and b where finite.
#' Infinite b or peak coordinates have no finite Wald SE/interval. Symmetric
#' Wald inference is unreliable at an active shape boundary; inspect the shape,
#' data coverage and simpler fixed-shape alternatives.
#'
#' ## Double logistic
#'
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
#' @return Formula markers; direct calls raise an error.
#' @examples
#' ~ logistic(age)
#' ~ warped_normal(age)                 # estimate the shape
#' ~ warped_normal(age, b = 0)          # asymptotic
#' ~ warped_normal(age, b = Inf)        # Gaussian limit
#' ~ double_logistic(age, by = fleet)
#' ~ 0 + for_fleet("Trawl", ~ factor(year) + logistic(age)) +
#'   for_fleet("Gillnet", ~ factor(year) + warped_normal(age))
#' @seealso [catchability_curves], [for_fleet()], [fit_tam()], [tidy_par()]
#' @export
warped_normal <- function(x, by = NULL, b = NULL) {
  cli::cli_abort("warped_normal() is a formula marker for F_settings$mu_form or index_settings$q_form.")
}

#' @rdname warped_normal
#' @export
double_logistic <- function(x, by = NULL) {
  cli::cli_abort("double_logistic() is a formula marker for F_settings$mu_form or index_settings$q_form.")
}

.is_selectivity_curve <- function(term) term$type %in% c("logistic", "warped_normal", "double_logistic")

# Parameter names keep q/F roles explicit, including inside scoped formulas.
.curve_parameter_names <- function(term, component = "q") {
  id <- term$id
  fields <- switch(term$type,
    logistic = c(a50 = "a50", slope = "log_slope"),
    warped_normal = c(a50 = "a50", scale = "log_scale",
      if (is.null(term$b)) c(shape = "shape")),
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
    warped_normal = c(list(a50 = mid, scale = log(width / 4)),
      if (is.null(term$b)) list(shape = rep(.5, length(mid)))),
    double_logistic = list(a50_up = mid - width / 4, scale_up = log(width / 8),
      a50_down = mid + width / 4, scale_down = log(width / 8)))
  stats::setNames(lapply(values, stats::setNames, nm = labels),
    vapply(names(values), function(x) term[[paste0(x, "_parameter")]], character(1)))
}

# Taped branching (RTMB >= 2.0), not an R branch fixed at the starting shape.
# The series has the exact continuous first/second derivatives at c = 0.
.warped_g <- function(u, shape) {
  scalar <- function(x, shape) {
    "if" <- RTMB::ADoverload("if")
    z <- shape * x
    if (abs(z) < .05) {
      1 - x * (1 - z/2 + z^2/6 - z^3/24 + z^4/120 - z^5/720 +
        z^6/5040 - z^7/40320 + z^8/362880)
    } else 1 + expm1(-z) / shape
  }
  RTMB::Vectorize(scalar)(u, shape)
}

.curve_log_selectivity <- function(x, par, term, group) {
  get <- function(field) par[[term[[paste0(field, "_parameter")]]]][group]
  if (term$type == "logistic") {
    predictor <- exp(get("slope")) * (x - get("a50"))
    return(-RTMB::logspace_add(0, -predictor))
  }
  if (term$type == "warped_normal") {
    u <- (x - get("a50")) / exp(get("scale"))
    shape <- if (is.null(term$b)) get("shape") else 1 / (1 + term$b)
    return(-log(2) * .warped_g(u, shape)^2)
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

.warped_peak_multiplier <- function(shape) {
  "if" <- RTMB::ADoverload("if")
  if (abs(shape) < .01) {
    1 + shape/2 + shape^2/3 + shape^3/4 + shape^4/5 + shape^5/6 +
      shape^6/7 + shape^7/8 + shape^8/9
  } else -log1p(-shape)/shape
}

.warped_tail <- function(shape) {
  "if" <- RTMB::ADoverload("if")
  if (shape == 0) 0*shape else exp(-log(2)*((1-shape)/shape)^2)
}

# Derived curve quantities use the joint fitted fixed-parameter covariance,
# including covariances between midpoint, scale and shape.
.tidy_selectivity_curves <- function(fit, interval) {
  terms <- Filter(function(x) x$type == "warped_normal", c(fit$dat$q_terms, fit$dat$F_terms))
  if (!length(terms)) return(data.frame())
  p <- .tam_parameter_summary(fit, "Estimate")
  sdr <- fit[["sdrep"]]
  has_cov <- !is.null(sdr) && isTRUE(sdr$pdHess) && !is.null(sdr$cov.fixed)
  z <- stats::qnorm(.5 + interval/2)
  out <- list()
  for (term in terms) for (i in seq_along(term$groups)) {
    names <- c(term$a50_parameter, term$scale_parameter, term$shape_parameter)
    theta <- vapply(names, function(nm) p[[nm]][i], numeric(1))
    shape <- if (is.null(term$b)) theta[3L] else 1/(1+term$b)
    # Infinite b/peak cannot have finite Wald uncertainty.
    metrics <- c("shape", "b", "peak", "tail")
    if (is.null(term$b)) metrics <- setdiff(metrics, "shape")
    finite <- metrics[!(metrics == "b" & shape == 0) & !(metrics == "peak" & shape == 1)]
    evaluate <- function(x) {
      "c" <- RTMB::ADoverload("c")
      shape_value <- if (is.null(term$b)) x[3L] else 1/(1+term$b)
      values <- list(shape = shape_value, tail = .warped_tail(shape_value))
      if ("b" %in% finite) values$b <- (1-shape_value)/shape_value
      if ("peak" %in% finite) values$peak <- x[1L] + exp(x[2L])*.warped_peak_multiplier(shape_value)
      RTMB::AD(do.call(c, values[finite]))
    }
    value <- as.numeric(evaluate(theta))
    se <- rep(NA_real_, length(value))
    fixed <- finite %in% c("shape", "b", "tail") & !is.null(term$b)
    se[fixed] <- 0
    if (has_cov) {
      positions <- vapply(names, function(nm) which(names(sdr$par.fixed) == nm)[i], integer(1))
      if (all(!is.na(positions))) {
        jac <- RTMB::MakeTape(evaluate, theta)$jacobian(theta)
        variance <- rowSums((jac %*% sdr$cov.fixed[positions, positions, drop = FALSE])*jac)
        se <- sqrt(pmax(0, variance))
      }
    }
    component <- if (startsWith(term$id, "F_")) "F" else "q"
    tab <- data.frame(par = paste0(component,"_curve_",finite), term = term$id,
      coef = term$groups[[i]]$group, est = value, se = se,
      lwr = value-z*se, upr = value+z*se, se_scale = "reported", derived = TRUE,
      is_fixed = fixed, unit = ifelse(finite %in% c("shape", "b", "tail"), "unitless", NA_character_))
    bounded <- finite %in% c("shape", "tail")
    tab$lwr[bounded] <- pmax(0,tab$lwr[bounded]); tab$upr[bounded] <- pmin(1,tab$upr[bounded])
    if ("b" %in% finite) tab$lwr[finite == "b"] <- pmax(0,tab$lwr[finite == "b"])
    infinite <- setdiff(metrics, finite)
    if (length(infinite)) tab <- rbind(tab, data.frame(par = paste0(component,"_curve_",infinite),
      term = term$id, coef = term$groups[[i]]$group, est = Inf, se = NA_real_,
      lwr = NA_real_, upr = NA_real_, se_scale = "reported", derived = TRUE,
      is_fixed = !is.null(term$b), unit = ifelse(infinite == "b","unitless",NA_character_)))
    if (!is.null(term$scope)) tab[[term$scope$type]] <- term$scope$name
    out[[length(out)+1L]] <- tab
  }
  .bind_rows(out)
}
