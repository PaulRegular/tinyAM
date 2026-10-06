.assessment_age_summary <- function(table, ages, weights = NULL) {
  if (is.null(table) || !all(c("year", "age", "est") %in% names(table))) return(NULL)
  table <- table[table$age %in% ages, c("year", "age", "est"), drop = FALSE]
  if (!nrow(table) || anyDuplicated(paste(table$year, table$age))) return(NULL)
  rows <- lapply(split(table, table$year), function(x) {
    if (!setequal(x$age, ages) || any(!is.finite(x$est))) return(NULL)
    value <- sum(x$est)
    if (!is.null(weights)) {
      w <- weights$est[match(paste(x$year, x$age), paste(weights$year, weights$age))]
      if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) return(NULL)
      value <- stats::weighted.mean(x$est, w)
    }
    data.frame(year = x$year[1L], est = value)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) do.call(rbind, rows) else NULL
}

.assessment_comparison_unit <- function(metric, scale) {
  if (metric %in% c("F", "M", "F_bar", "M_bar", "Z_bar")) return("per year")
  numbers <- metric %in% c("N", "recruitment", "abundance")
  units <- if (numbers) c("1" = "fish", "0.001" = "thousand fish",
                          "1e-06" = "million fish", "1e-09" = "billion fish") else
    c("1" = "kg", "0.001" = "t", "1e-06" = "kt")
  unit <- unname(units[as.character(scale)])
  if (is.na(unit)) paste(if (numbers) "fish" else "kg", "x", scale) else unit
}

.assessment_percent_differences <- function(fit, reference,
                                             scales = c(ssb = 1e-3,
                                               recruitment = 1e-3, N = 1e-3,
                                               F = 1, M = 1,
                                               abundance = 1e-3, biomass = 1e-3,
                                               F_bar = 1, M_bar = 1),
                                             assumptions = NULL,
                                             comparison_aggregates = character()) {
  if (!is.null(reference$comparison_scales)) {
    scales[names(reference$comparison_scales)] <- reference$comparison_scales
  }
  native <- attr(reference, "source_pop")
  if (is.null(native)) native <- reference$pop
  # Convert source quantities to fish/kg before constructing common definitions.
  for (name in names(native)) {
    x <- native[[name]]
    if (is.null(x$unit)) next
    convert <- if (name %in% c("N", "recruitment", "abundance")) {
      function(value, unit) value * .translation_number_multiplier(unit)
    } else if (name %in% c("biomass", "biomass_at_age", "ssb", "ssb_mat")) {
      .translation_biomass_to_kg
    } else NULL
    if (!is.null(convert)) {
      converted <- tryCatch(mapply(convert, x$est, x$unit), error = function(e) NULL)
      if (is.null(converted)) native[[name]] <- NULL else native[[name]]$est <- converted
    }
  }
  ages <- fit$dat$ages
  rows <- lapply(names(scales), function(metric) {
    source <- native[[metric]]
    tiny <- fit$pop[[metric]]
    if (metric == "biomass_at_age" && is.null(tiny)) tiny <- fit$pop$biomass_mat
    unit <- .assessment_comparison_unit(metric, scales[[metric]])
    definition <- "Matching annual age-specific states"
    reason <- "No numerical output is available on the matching definition."
    status <- "unavailable"
    if (metric == "recruitment") {
      source_age <- if (!is.null(source$age)) unique(na.omit(source$age)) else integer()
      if (!length(source_age) && !is.null(assumptions)) {
        values <- unique(assumptions$value[assumptions$setting == "recruitment_age"])
        if (length(values) == 1L && grepl("^[0-9]+$", values)) source_age <- as.integer(values)
      }
      if (length(source_age) != 1L || source_age != min(ages)) {
        source <- NULL
        status <- "non_equivalent"
        reason <- "Recruitment ages differ or the accepted recruitment age is unresolved."
      } else definition <- paste("Recruitment at age", source_age)
    }
    if (metric %in% c("abundance", "biomass", "ssb")) {
      age_metric <- switch(metric, abundance = "N", biomass = "biomass_at_age", ssb = "ssb_mat")
      source <- .assessment_age_summary(native[[age_metric]], ages)
      definition <- paste("Sum of reported", age_metric, "over ages", paste(range(ages), collapse = "-"))
      # When no biomass surface is reported, compare N using the same biological
      # inputs in both models. This is a common-definition calculation, not a
      # replacement for the accepted model's native biomass or spawning biomass.
      if (is.null(source) && metric %in% c("biomass", "ssb") &&
          !is.null(native$N) && !is.null(fit$dat$W)) {
        n <- native$N
        index <- cbind(match(n$year, fit$dat$years), match(n$age, ages))
        biology <- fit$dat$W[index]
        if (metric == "ssb") biology <- biology * fit$dat$P[index]
        n$est <- n$est * biology
        source <- .assessment_age_summary(n, ages)
        definition <- paste("Accepted N with shared translated weights",
                            if (metric == "ssb") "and maturity" else "",
                            "over ages", paste(range(ages), collapse = "-"))
      }
      if (is.null(source) && metric %in% comparison_aggregates &&
          metric == "ssb" && !is.null(native$ssb) &&
          (!"age" %in% names(native$ssb) || all(is.na(native$ssb$age)))) {
        source <- native$ssb[c("year", "est")]
        definition <- "Reported aggregate female SSB; source age-specific contributions are unavailable"
      }
      if (is.null(source) && !is.null(native[[metric]])) {
        status <- "non_equivalent"
        reason <- "Native aggregate cannot be reconciled to the modeled ages and biological definition."
      }
    }
    if (metric %in% c("F_bar", "M_bar")) {
      rate <- if (metric == "F_bar") "F" else "M"
      mean_ages <- fit$dat[[paste0(rate, "_settings")]]$mean_ages
      source <- if (is.null(native$N)) NULL else
        .assessment_age_summary(native[[rate]], mean_ages, native$N)
      definition <- paste("N-weighted", rate, "over ages", paste(mean_ages, collapse = ","))
      if (is.null(source) && !is.null(native[[metric]])) {
        status <- "non_equivalent"
        reason <- "Native mortality mean lacks a recoverable matching age range and population weighting."
      }
    }
    keys <- if (metric %in% c("N", "F", "M", "biomass_at_age")) c("year", "age") else "year"
    if (!is.null(source) && !is.null(tiny) &&
        all(c(keys, "est") %in% names(source)) && all(c(keys, "est") %in% names(tiny))) {
      source <- source[c(keys, "est")]
      tiny <- tiny[c(keys, "est")]
      if (!anyDuplicated(source[keys]) && !anyDuplicated(tiny[keys])) {
        names(source)[ncol(source)] <- "source"
        names(tiny)[ncol(tiny)] <- "tinyAM"
        common <- merge(source, tiny, by = keys)
        if (nrow(common) && any(is.finite(common$source) & is.finite(common$tinyAM))) {
          if (!"age" %in% names(common)) common$age <- NA_integer_
          common$source <- common$source * scales[[metric]]
          common$tinyAM <- common$tinyAM * scales[[metric]]
          common$metric <- metric
          common$unit <- unit
          common$comparison_status <- "matched"
          common$definition <- definition
          common$reason <- ""
          common$absolute_difference <- abs(common$tinyAM - common$source)
          common$percent_difference <- ifelse(is.finite(common$source) & common$source != 0,
            100 * (common$tinyAM - common$source) / common$source, NA_real_)
          return(common[c("metric", "year", "age", "source", "tinyAM", "absolute_difference",
                           "percent_difference", "unit", "comparison_status", "definition", "reason")])
        }
      } else reason <- "Output dimensions do not map uniquely to matching year-age groups."
    }
    data.frame(metric = metric, year = NA_integer_, age = NA_integer_, source = NA_real_,
               tinyAM = NA_real_, absolute_difference = NA_real_, percent_difference = NA_real_,
               unit = unit, comparison_status = status, definition = definition, reason = reason)
  })
  do.call(rbind, rows)
}
