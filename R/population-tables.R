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
