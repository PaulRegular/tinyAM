#' @rdname catchability_curves
#' @export
double_logistic <- function(x, by = NULL) {
  cli::cli_abort("double_logistic() is a formula marker for F_settings$mu_form or index_settings$q_form.")
}

#' @rdname catchability_curves
#' @export
normal <- function(x, by = NULL) {
  cli::cli_abort("normal() is a formula marker for F_settings$mu_form or index_settings$q_form.")
}

.is_selectivity_curve <- function(term) term$type %in% c("logistic", "normal", "double_logistic")

# Parameter names keep q/F roles explicit, including inside scoped formulas.
.curve_parameter_names <- function(term, component = "q") {
  id <- term$id
  fields <- switch(term$type,
    logistic = c(a50 = "a50", slope = "log_slope"),
    normal = c(peak = "peak", scale = "log_scale"),
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
    normal = list(peak = mid, scale = log(width / 4)),
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
  if (term$type == "normal") {
    distance <- (x - get("peak")) / exp(get("scale"))
    return(-0.5 * distance^2)
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
