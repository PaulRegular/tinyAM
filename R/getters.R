#' Extract common population estimates
#'
#' Return a ready-to-plot table of spawning biomass, recruitment or mortality,
#' without navigating the nested tidy output. These functions select the same
#' estimates, uncertainty and units used in [tidy_tam()] and the dashboard.
#'
#' @param ... One or more fitted `tam_fit` models or `tam_ref` reporting
#'   references, or one named list of these objects.
#' @param model_list A named list of models, instead of `...`.
#' @param interval Confidence level in `(0, 1)`, default `0.95`. Fitted summaries
#'   are rebuilt when their cached interval differs. Reference uncertainty is
#'   retained as supplied, because it cannot be recalculated from a reporting
#'   reference.
#' @details
#' `get_ssb()` and `get_recruitment()` return annual estimates. `get_fbar()`
#' and `get_mbar()` return the model's abundance-weighted means over the ages
#' chosen by `F_settings$mean_ages` and `M_settings$mean_ages`, respectively.
#' `get_f()` and `get_m()` return year-by-age surfaces; F is total fishing
#' mortality, summed over fleets. They do not return process deviations.
#'
#' A single fit returns its table without a model column. Multiple fits and
#' named-list inputs add `model`, using the supplied names. Missing reference
#' quantities are omitted rather than reconstructed. Explicit reference
#' `comparison_scales` follow [tidy_tam()]'s conversion rules; otherwise known
#' incompatible units must be converted before combining models.
#'
#' These are lightweight selectors. They never refit models or calculate
#' observation OSA residuals or joint-sample process diagnostics. Supplied
#' objects are not modified.
#'
#' @return A data frame retaining available keys, `est`, confidence limits,
#'   `se`, `se_scale`, `unit` and `is_proj`. Available metadata, such as SSB
#'   timing, is retained. Annual summaries have `year`; surfaces also have `age`.
#' @examples
#' \dontrun{
#' get_ssb(fit)
#' get_recruitment(fit, interval = 0.8)
#' get_fbar(Base = fit, Alternative = other_fit)
#' get_f(model_list = list(Base = fit, Alternative = other_fit))
#' }
#' @seealso [tidy_pop()], [tidy_tam()], [plot_trend()], [plot_matrix()]
#' @name get_population
NULL

#' @rdname get_population
#' @export
get_ssb <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("ssb", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

#' @rdname get_population
#' @export
get_recruitment <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("recruitment", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

#' @rdname get_population
#' @export
get_fbar <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("F_bar", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

#' @rdname get_population
#' @export
get_mbar <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("M_bar", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

#' @rdname get_population
#' @export
get_f <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("F", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

#' @rdname get_population
#' @export
get_m <- function(..., model_list = NULL, interval = 0.95) {
  .get_population_quantity("M", list(...), as.list(substitute(list(...)))[-1L],
    model_list, interval)
}

.get_population_quantity <- function(quantity, dots, expressions, model_list, interval) {
  if (!is.numeric(interval) || length(interval) != 1L || !is.finite(interval) ||
      interval <= 0 || interval >= 1) {
    cli::cli_abort("{.arg interval} must be one number strictly between zero and one.")
  }
  info <- .dots_or_list(dots, expressions, model_list, "model_list")
  named_list <- length(dots) == 1L && is.list(dots[[1L]]) &&
    !.is_reporting_fit(dots[[1L]])
  add_label <- !info$using_dots || named_list || length(info$fits) > 1L
  populations <- .population_tables(info$fits, interval)
  selected <- lapply(populations, function(x) x[intersect(names(x), quantity)])
  tab <- stack_nested(selected, label = if (add_label) "model" else NULL)[[quantity]]
  if (!is.data.frame(tab)) {
    cli::cli_abort("The requested quantity {.val {quantity}} is unavailable in the supplied models.")
  }
  .unit_label(tab, "")
  tab
}
