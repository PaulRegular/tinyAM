.assessment_age_summary <- function(table, ages, weights = NULL, model_ages = ages) {
  if (is.null(table) || !all(c("year", "age", "est") %in% names(table))) return(NULL)
  resolve_age <- function(x) {
    age <- suppressWarnings(as.integer(as.character(x$age)))
    if ("age_group" %in% names(x) && length(model_ages)) {
      group <- as.character(x$age_group)
      plus <- !is.na(group) & grepl("^[0-9]+\\+$", group)
      start <- suppressWarnings(as.integer(sub("\\+$", "", group)))
      age[plus & start == max(model_ages)] <- max(model_ages)
      age[plus & start != max(model_ages)] <- NA_integer_
    }
    age
  }
  table$.age <- resolve_age(table)
  table <- table[table$.age %in% ages, c("year", ".age", "est"), drop = FALSE]
  names(table)[names(table) == ".age"] <- "age"
  if (!nrow(table) || anyDuplicated(paste(table$year, table$age))) return(NULL)
  if (!is.null(weights)) {
    if (!all(c("year", "age", "est") %in% names(weights))) return(NULL)
    weights$.age <- resolve_age(weights)
    weights <- weights[weights$.age %in% ages, c("year", ".age", "est"), drop = FALSE]
    names(weights)[names(weights) == ".age"] <- "age"
    if (anyDuplicated(paste(weights$year, weights$age))) return(NULL)
  }
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

.assessment_age_comparison_cells <- function(source, tiny, ages, groups = NULL) {
  if (is.null(source) || is.null(tiny) ||
      !all(c("year", "est") %in% names(source)) ||
      !all(c("year", "age", "est") %in% names(tiny))) {
    return(list(rows = NULL, reason = "An age-specific source or tinyAM surface is unavailable."))
  }
  if (is.null(groups)) groups <- list()
  if (!is.list(groups) || (length(groups) &&
      (is.null(names(groups)) || anyNA(names(groups)) || any(!nzchar(names(groups))) ||
       anyDuplicated(names(groups))))) {
    stop("comparison_age_groups must be a named list of reported age labels and tinyAM age vectors.")
  }
  if (length(groups)) {
    valid_groups <- vapply(groups, function(x) {
      is.numeric(x) && length(x) > 0L && !anyNA(x) && all(is.finite(x)) &&
        !anyDuplicated(x) && all(x %in% ages)
    }, logical(1))
    mapped <- unlist(groups, use.names = FALSE)
    if (!all(valid_groups) || anyDuplicated(mapped)) {
      stop("comparison_age_groups must map each reported group to distinct modeled ages.")
    }
  }
  source_label <- if ("age" %in% names(source)) as.character(source$age) else
    rep(NA_character_, nrow(source))
  if ("age_group" %in% names(source)) {
    group <- as.character(source$age_group)
    has_plus <- !is.na(group) & grepl("^[0-9]+\\+$", group)
    source_label[has_plus] <- group[has_plus]
  }
  plus <- !is.na(source_label) & grepl("^[0-9]+\\+$", source_label)
  plus_age <- suppressWarnings(as.integer(sub("\\+$", "", source_label)))
  ordinary_plus <- plus & !source_label %in% names(groups) & plus_age == max(ages)
  source_label[ordinary_plus] <- as.character(max(ages))
  unknown <- !is.na(source_label) & !grepl("^[0-9]+$", source_label) &
    !source_label %in% names(groups)
  if (any(unknown)) {
    labels <- paste(unique(source_label[unknown]), collapse = ", ")
    return(list(rows = NULL, reason = paste0(
      "Reported age group(s) ", labels,
      " have no explicit comparison_age_groups mapping."
    )))
  }
  if (anyDuplicated(paste(source$year, source_label))) {
    return(list(rows = NULL, reason = "Source age groups do not map uniquely within year."))
  }
  if (anyDuplicated(paste(tiny$year, tiny$age))) {
    return(list(rows = NULL, reason = "tinyAM age rows do not map uniquely within year."))
  }

  rows <- list()
  mapped_ages <- if (length(groups)) unlist(groups, use.names = FALSE) else numeric()
  direct_age <- suppressWarnings(as.integer(source_label))
  direct <- !is.na(direct_age) & direct_age %in% ages & !direct_age %in% mapped_ages
  if (any(direct)) {
    selected <- which(direct)
    index <- match(paste(source$year[selected], direct_age[selected]),
                   paste(tiny$year, tiny$age))
    keep <- which(!is.na(index) & is.finite(source$est[selected]) &
                    is.finite(tiny$est[index]))
    if (length(keep)) {
      rows[[length(rows) + 1L]] <- data.frame(
        year = source$year[selected[keep]], age = as.character(direct_age[selected[keep]]),
        source = source$est[selected[keep]], tinyAM = tiny$est[index[keep]],
        definition = paste0("Reported age ", direct_age[selected[keep]],
                            " compared with tinyAM age ", direct_age[selected[keep]]),
        stringsAsFactors = FALSE
      )
    }
  }
  for (label in names(groups)) {
    source_group <- which(source_label == label)
    if (!length(source_group)) next
    for (year in unique(source$year[source_group])) {
      selected <- source_group[source$year[source_group] == year]
      if (length(selected) != 1L || !is.finite(source$est[selected])) next
      tiny_group <- tiny[tiny$year == year & tiny$age %in% groups[[label]], , drop = FALSE]
      if (nrow(tiny_group) != length(groups[[label]]) ||
          !setequal(tiny_group$age, groups[[label]]) || any(!is.finite(tiny_group$est))) next
      age_text <- if (length(groups[[label]]) > 1L &&
                      all(diff(sort(groups[[label]])) == 1L)) {
        paste(range(groups[[label]]), collapse = "-")
      } else {
        paste(groups[[label]], collapse = ", ")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        year = year, age = label, source = source$est[selected],
        tinyAM = sum(tiny_group$est),
        definition = paste0("Reported age group ", label,
                            " compared with the sum of tinyAM ages ",
                            age_text),
        stringsAsFactors = FALSE
      )
    }
  }
  list(rows = if (length(rows)) do.call(rbind, rows) else NULL,
       reason = if (length(rows)) "" else "No complete matching age groups were available.")
}

.assessment_age_group_total <- function(cells, ages, groups = NULL) {
  if (is.null(cells) || !nrow(cells)) return(NULL)
  if (is.null(groups)) groups <- list()
  mapped <- if (length(groups)) unlist(groups, use.names = FALSE) else numeric()
  expected <- c(as.character(ages[!ages %in% mapped]), names(groups))
  rows <- lapply(split(cells, cells$year), function(x) {
    if (anyDuplicated(x$age) || !setequal(x$age, expected) ||
        any(!is.finite(x$source)) || any(!is.finite(x$tinyAM))) return(NULL)
    data.frame(year = x$year[1L], source = sum(x$source), tinyAM = sum(x$tinyAM))
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

.tam_reference_comparisons <- function(fit, reference,
                                             scales = c(ssb = 1e-3,
                                               recruitment = 1e-3, N = 1e-3,
                                               F = 1, M = 1,
                                               abundance = 1e-3, biomass = 1e-3,
                                               biomass_at_age = 1e-3,
                                               F_bar = 1, M_bar = 1),
                                             assumptions = NULL,
                                             comparison_aggregates = character(),
                                             comparison_age_groups = list(),
                                             comparison_definitions = list()) {
  if (!is.null(reference$comparison_scales)) {
    scales[names(reference$comparison_scales)] <- reference$comparison_scales
  }
  native <- attr(reference, "source_pop")
  if (is.null(native)) native <- reference$pop
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
    comparison_definition <- comparison_definitions[[metric]]
    if (metric %in% c("N", "F", "M", "biomass_at_age", "ssb_mat") &&
        !is.null(source) && !is.null(tiny)) {
      cells <- .assessment_age_comparison_cells(
        source, tiny, ages, comparison_age_groups[[metric]]
      )
      if (!is.null(cells$rows) && nrow(cells$rows)) {
        common <- cells$rows
        common$source <- common$source * scales[[metric]]
        common$tinyAM <- common$tinyAM * scales[[metric]]
        common$metric <- metric
        common$unit <- unit
        common$comparison_status <- if (is.null(comparison_definition$status)) "matched" else
          comparison_definition$status
        common$reason <- if (is.null(comparison_definition$reason)) "" else
          comparison_definition$reason
        common$definition <- if (is.null(comparison_definition$definition))
          common$definition else comparison_definition$definition
        return(common[c("metric", "year", "age", "source", "tinyAM",
                         "unit",
                         "comparison_status", "definition", "reason")])
      }
      status <- if (grepl("no explicit comparison_age_groups", cells$reason, fixed = TRUE))
        "non_equivalent" else "unavailable"
      reason <- cells$reason
      source <- tiny <- NULL
    }
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
      if (metric == "abundance") {
        cells <- .assessment_age_comparison_cells(
          native$N, fit$pop$N, ages, comparison_age_groups[["N"]]
        )
        totals <- .assessment_age_group_total(cells$rows, ages,
                                              comparison_age_groups[["N"]])
        if (!is.null(totals)) {
          source <- data.frame(year = totals$year, est = totals$source)
          tiny <- data.frame(year = totals$year, est = totals$tinyAM)
          definition <- paste0("Total abundance across modeled ages ",
                               paste(range(ages), collapse = "-"),
                               "; explicit reported age groups are matched to sums of tinyAM ages.")
        } else {
          source <- tiny <- NULL
          status <- "non_equivalent"
          reason <- if (!is.null(cells$reason) && nzchar(cells$reason)) cells$reason else
            "Accepted age groups do not cover the modeled ages completely for total abundance."
        }
      } else {
        age_metric <- if (metric == "biomass") "biomass_at_age" else "ssb_mat"
        source <- .assessment_age_summary(native[[age_metric]], ages,
                                           model_ages = fit$dat$ages)
        definition <- paste("Sum of reported", age_metric, "over ages", paste(range(ages), collapse = "-"))
        spawn_time <- fit$dat$ssb_settings$spawn_time
        if (is.null(spawn_time)) spawn_time <- 0
        if (is.null(source) && metric %in% c("biomass", "ssb") &&
            (metric != "ssb" || spawn_time == 0) &&
            !is.null(native$N) && !is.null(fit$dat$W)) {
          n <- native$N
          index <- cbind(match(n$year, fit$dat$years), match(n$age, ages))
          biology <- fit$dat$W[index]
          if (metric == "ssb") biology <- biology * fit$dat$P[index]
          n$est <- n$est * biology
          source <- .assessment_age_summary(n, ages, model_ages = fit$dat$ages)
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
          status <- if (is.null(comparison_definition$status)) {
            "non_equivalent"
          } else {
            comparison_definition$status
          }
          definition <- if (is.null(comparison_definition$definition)) {
            definition
          } else {
            comparison_definition$definition
          }
          reason <- if (is.null(comparison_definition$reason)) {
            "Native aggregate cannot be reconciled to the modeled ages and biological definition."
          } else {
            comparison_definition$reason
          }
        }
      }
    }
    if (metric %in% c("F_bar", "M_bar")) {
      rate <- if (metric == "F_bar") "F" else "M"
      mean_ages <- fit$dat[[paste0(rate, "_settings")]]$mean_ages
      source <- if (is.null(native$N)) NULL else
        .assessment_age_summary(native[[rate]], mean_ages, native$N,
                                model_ages = fit$dat$ages)
      definition <- paste("N-weighted", rate, "over ages", paste(mean_ages, collapse = ","))
      if (is.null(source) && !is.null(native[[metric]])) {
        status <- "non_equivalent"
        reason <- "Native mortality mean lacks a recoverable matching age range and population weighting."
      }
    }
    keys <- if (metric %in% c("N", "F", "M", "biomass_at_age", "ssb_mat")) c("year", "age") else "year"
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
          common$comparison_status <- if (is.null(comparison_definition$status)) "matched" else
            comparison_definition$status
          common$definition <- if (is.null(comparison_definition$definition)) definition else
            comparison_definition$definition
          common$reason <- if (is.null(comparison_definition$reason)) "" else
            comparison_definition$reason
          return(common[c("metric", "year", "age", "source", "tinyAM",
                           "unit", "comparison_status", "definition", "reason")])
        }
      } else reason <- "Output dimensions do not map uniquely to matching year-age groups."
    }
    data.frame(metric = metric, year = NA_integer_, age = NA_character_, source = NA_real_,
               tinyAM = NA_real_,
               unit = unit, comparison_status = status, definition = definition, reason = reason,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "scales") <- scales
  out
}

.assessment_percent_differences <- function(reference) {
  common <- reference$comparisons
  if (is.null(common)) {
    cli::cli_abort("Translate the accepted reference with a fitted template before calculating differences.")
  }
  common$absolute_difference <- abs(common$tinyAM - common$source)
  common$percent_difference <- ifelse(is.finite(common$source) & common$source != 0,
    100 * (common$tinyAM - common$source) / common$source, NA_real_)
  common[c("metric", "year", "age", "source", "tinyAM", "absolute_difference",
           "percent_difference", "unit", "comparison_status", "definition", "reason")]
}
