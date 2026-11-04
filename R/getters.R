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

# Population-only collection shared by tidy_tam() and the quantity getters.
# Keep comparison conversions here so lightweight consumers use the same units
# and reported-versus-log SE conventions as dashboards.
.population_tables <- function(model_list, interval) {
  pop_list <- lapply(model_list, function(fit) {
    if (inherits(fit, "tam_ref")) return(fit$pop)
    stored <- attr(fit$pop, "interval", exact = TRUE)
    if (!is.null(fit$pop) && !is.null(stored) &&
        isTRUE(all.equal(stored, interval))) {
      fit$pop
    } else {
      tidy_pop(fit, interval = interval)
    }
  })
  references <- Filter(function(fit) inherits(fit, "tam_ref"), model_list)
  comparison_scales <- lapply(references, `[[`, "comparison_scales")
  comparison_scales <- Filter(function(x) length(x) > 0L, comparison_scales)
  if (length(comparison_scales)) {
    comparison_scales <- lapply(comparison_scales, function(x) {
      if (is.list(x)) x <- unlist(x, use.names = TRUE)
      if (!is.numeric(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
          anyDuplicated(names(x)) || any(!is.finite(x)) || any(x <= 0)) {
        cli::cli_abort("A {.cls tam_ref} {.field comparison_scales} value must be a named vector of positive finite numbers.")
      }
      x
    })
    if (length(comparison_scales) > 1L &&
        !all(vapply(comparison_scales[-1L], identical, logical(1), comparison_scales[[1L]]))) {
      cli::cli_abort("Assessment references use different {.field comparison_scales}; the fitted outputs cannot be shown on one scale.")
    }
    comparison_scales <- comparison_scales[[1L]]
    reference_pop <- references[[1L]]$pop
    for (i in seq_along(model_list)) {
      if (inherits(model_list[[i]], "tam_ref")) next
      for (metric in intersect(names(comparison_scales), names(pop_list[[i]]))) {
        tab <- pop_list[[i]][[metric]]
        scale <- comparison_scales[[metric]]
        reference <- reference_pop[[metric]]
        from <- unique(tab$unit[!is.na(tab$unit)])
        to <- unique(reference$unit[!is.na(reference$unit)])
        if (length(from) == 1L && length(to) == 1L) .check_unit_conversion(from, to, scale)
        for (field in intersect(c("est", "lwr", "upr"), names(tab))) {
          tab[[field]] <- tab[[field]] * scale
        }
        if ("se" %in% names(tab)) {
          natural <- if ("se_scale" %in% names(tab)) tab$se_scale == "reported" else rep(TRUE, nrow(tab))
          tab$se[which(natural)] <- tab$se[which(natural)] * scale
        }
        if (is.data.frame(reference) && "unit" %in% names(reference)) {
          units <- unique(as.character(reference$unit[!is.na(reference$unit) &
            nzchar(as.character(reference$unit))]))
          if (length(units) == 1L) tab$unit <- units[[1L]]
        }
        pop_list[[i]][[metric]] <- tab
      }
    }
  }
  pop_list
}
