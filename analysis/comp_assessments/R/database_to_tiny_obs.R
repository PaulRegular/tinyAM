.translation_abort <- function(message) {
  cli::cli_abort(message, call = NULL)
}

.translation_rows <- function(inputs, assessment_id) {
  required <- c("assessment_id", "type", "measure", "basis", "survey", "fleet",
                "sex", "region", "season", "year", "year_basis", "age",
                "value", "unit", "sampling_time", "source_type",
                "source_reference", "transformation", "notes")
  missing <- setdiff(required, names(inputs))
  if (length(missing)) {
    .translation_abort(paste("inputs is missing required columns:",
                             paste(missing, collapse = ", ")))
  }
  if (length(assessment_id) != 1L || is.na(assessment_id) ||
      !nzchar(assessment_id)) {
    .translation_abort("assessment_id must be one non-empty value.")
  }

  x <- inputs[!is.na(inputs$assessment_id) &
                inputs$assessment_id == assessment_id, , drop = FALSE]
  if (!nrow(x)) {
    .translation_abort(paste("No input rows found for assessment_id:", assessment_id))
  }
  x$year <- suppressWarnings(as.numeric(as.character(x$year)))
  x$age <- suppressWarnings(as.numeric(as.character(x$age)))
  x$value <- suppressWarnings(as.numeric(as.character(x$value)))
  x
}

.translation_measure <- function(x, type, measure) {
  x[!is.na(x$type) & x$type == type & !is.na(x$measure) &
      x$measure == measure, , drop = FALSE]
}

.translation_one_group <- function(x, columns, label) {
  for (column in columns) {
    values <- unique(as.character(x[[column]][!is.na(x[[column]]) &
                                               nzchar(as.character(x[[column]]))]))
    if (length(values) > 1L) {
      .translation_abort(paste(label, "contains multiple", column,
                               "groups; select or aggregate them explicitly."))
    }
  }
}

.translation_key <- function(x, columns) {
  do.call(paste, c(lapply(x[columns], function(z) {
    z <- as.character(z)
    z[is.na(z)] <- ""
    z
  }), sep = "\r"))
}

.translation_calendar_rows <- function(x, label) {
  if ("year_basis" %in% names(x)) {
    x <- x[!is.na(x$year_basis) & x$year_basis == "calendar_year", , drop = FALSE]
  }
  if (!nrow(x)) .translation_abort(paste("No calendar-year", label, "inputs are recorded."))
  if (anyNA(x$year) || anyNA(x$age) || anyNA(x$value) ||
      any(!is.finite(x$year)) || any(!is.finite(x$age)) ||
      any(x$year != as.integer(x$year)) || any(x$age != as.integer(x$age)) ||
      any(!is.finite(x$value)) || any(x$value < 0)) {
    .translation_abort(paste(label, "rows need finite, non-negative values and whole years and ages."))
  }
  x
}

.translation_grid <- function(years, ages) {
  expand.grid(year = years, age = ages)
}

.translation_source_provenance <- function(x, component, method) {
  collapse <- function(column) {
    values <- unique(as.character(x[[column]][!is.na(x[[column]]) &
                                                  nzchar(as.character(x[[column]]))]))
    paste(values, collapse = "; ")
  }
  data.frame(component = component, survey = NA_character_,
             source_type = collapse("source_type"), method = method,
             source_reference = collapse("source_reference"),
             notes = collapse("notes"), stringsAsFactors = FALSE)
}

.translation_surface <- function(x, years, ages, label) {
  x <- x[x$year %in% years & x$age %in% ages, , drop = FALSE]
  if (anyDuplicated(.translation_key(x, c("year", "age")))) {
    .translation_abort(paste(label, "has multiple values for a year-age row."))
  }
  grid <- .translation_grid(years, ages)
  i <- match(.translation_key(grid, c("year", "age")),
             .translation_key(x, c("year", "age")))
  if (anyNA(i)) {
    .translation_abort(paste(label, "does not cover the requested year-age grid."))
  }
  data.frame(year = grid$year, age = grid$age, obs = x$value[i])
}

.translation_biomass_to_kg <- function(value, unit) {
  unit <- tolower(trimws(as.character(unit)))
  if (grepl("kt", unit, fixed = TRUE)) return(value * 1e6)
  if (grepl("kg", unit, fixed = TRUE)) return(value)
  if (grepl("tonne", unit, fixed = TRUE) || grepl("(^|[^a-z])t($|[^a-z])", unit)) {
    return(value * 1000)
  }
  .translation_abort(paste("Cannot convert biomass unit to kg:", unit))
}

.translation_number_multiplier <- function(unit) {
  unit <- tolower(trimws(as.character(unit)))
  count_unit <- "(fish|individuals?|numbers?|counts?)"
  if (grepl(paste0("(million|10\\^?6|1,?000,?000)\\s*", count_unit), unit)) return(1e6)
  if (grepl(paste0("(billion|10\\^?9|1,?000,?000,?000)\\s*", count_unit), unit)) return(1e9)
  if (grepl(paste0("(thousand|10\\^?3|1,?000)\\s*", count_unit), unit)) return(1e3)
  if (grepl(count_unit, unit)) return(1)
  .translation_abort(paste("Cannot convert number unit to individual fish:", unit))
}

.translation_index_multiplier <- function(unit) {
  unit <- tolower(trimws(as.character(unit)))
  if (grepl("^(native[ _]+)?survey[ _]?index$", unit)) return(1)
  .translation_number_multiplier(unit)
}

.translation_catch_at_age <- function(x, years, ages) {
  direct <- .translation_measure(x, "catch", "numbers_at_age")
  composition <- .translation_measure(x, "catch", "proportion_at_age")
  totals <- .translation_measure(x, "catch", "total_numbers")
  if (!nrow(direct) && !nrow(composition)) {
    .translation_abort("No source catch-at-age or number-proportion inputs are recorded.")
  }
  if (nrow(direct)) direct <- .translation_calendar_rows(direct, "catch numbers-at-age")
  if (nrow(composition)) composition <- .translation_calendar_rows(composition, "catch age composition")
  source_rows <- rbind(direct, composition)
  .translation_one_group(source_rows, c("fleet", "sex", "region", "season"),
                         "Catch-at-age input")
  direct <- direct[direct$year %in% years & direct$age >= min(ages), , drop = FALSE]
  composition <- composition[composition$year %in% years &
                               composition$age >= min(ages), , drop = FALSE]
  totals <- totals[totals$year %in% years, , drop = FALSE]

  converted <- list()
  provenance <- list()
  units <- list()
  k <- 0L

  if (nrow(direct)) {
    multiplier <- vapply(direct$unit, .translation_number_multiplier, numeric(1))
    direct$value <- direct$value * multiplier
    if (anyDuplicated(.translation_key(direct, c("year", "age")))) {
      .translation_abort("Catch numbers-at-age has duplicate year-age rows before plus-group aggregation.")
    }
    converted[[length(converted) + 1L]] <- direct[c("year", "age", "value")]
    units[[length(units) + 1L]] <- data.frame(
      source_unit = direct$unit, multiplier_to_fish = multiplier,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    provenance[[k]] <- .translation_source_provenance(
      direct, "catch",
      "Source numbers-at-age retained and converted to individual fish; ages above the selected maximum summed into the plus age.")
  }

  if (nrow(composition)) {
    basis <- as.character(composition$basis)
    if (anyNA(basis) || !all(basis == "proportion_numbers")) {
      .translation_abort("Catch proportions can be reconstructed only when basis = 'proportion_numbers'.")
    }
    if (any(composition$value > 1)) {
      .translation_abort("Catch number proportions must be in [0, 1].")
    }
    if (!nrow(totals)) {
      .translation_abort("Catch number proportions need matching annual total_numbers inputs.")
    }
    totals <- totals[!is.na(totals$year_basis) & totals$year_basis == "calendar_year", , drop = FALSE]
    totals$year <- suppressWarnings(as.numeric(as.character(totals$year)))
    totals$value <- suppressWarnings(as.numeric(as.character(totals$value)))
    if (!nrow(totals) || anyNA(totals$year) || anyNA(totals$value) ||
        any(!is.finite(totals$year)) || any(!is.finite(totals$value)) ||
        any(totals$year != as.integer(totals$year)) || any(totals$value < 0) ||
        any(!is.na(totals$age))) {
      .translation_abort("Catch total_numbers rows need finite calendar years, non-negative values, and no age dimension.")
    }
    group_columns <- c("year", "fleet", "sex", "region", "season")
    composition_key <- .translation_key(composition, group_columns)
    total_key <- .translation_key(totals, group_columns)
    if (anyDuplicated(total_key)) {
      .translation_abort("Catch total_numbers has multiple rows for a year and catch stream.")
    }
    total_index <- match(composition_key, total_key)
    if (anyNA(total_index)) {
      .translation_abort("Some catch number-proportion rows have no matching total_numbers value for the same year and catch stream.")
    }
    composition_age_key <- .translation_key(composition, c(group_columns, "age"))
    if (anyDuplicated(composition_age_key)) {
      .translation_abort("Catch proportions contain duplicate year-age rows before plus-group aggregation.")
    }
    multiplier <- vapply(totals$unit[total_index], .translation_number_multiplier, numeric(1))
    composition$value <- composition$value * totals$value[total_index] * multiplier
    converted[[length(converted) + 1L]] <- composition[c("year", "age", "value")]
    units[[length(units) + 1L]] <- data.frame(
      source_unit = totals$unit[total_index], multiplier_to_fish = multiplier,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    source <- rbind(composition, totals[unique(total_index), , drop = FALSE])
    provenance[[k]] <- .translation_source_provenance(
      source, "catch",
      "Number proportions multiplied by matching annual total_numbers without renormalizing the source proportions; ages above the selected maximum summed into the plus age.")
  }

  catch_source <- if (length(converted)) {
    do.call(rbind, converted)
  } else {
    data.frame(year = numeric(), age = numeric(), value = numeric())
  }
  catch_key <- .translation_key(catch_source, c("year", "age"))
  if (anyDuplicated(catch_key)) {
    .translation_abort("Catch inputs contain overlapping year-age values from multiple representations.")
  }
  catch_source$age <- pmin(catch_source$age, max(ages))
  catch_source <- stats::aggregate(value ~ year + age, catch_source, sum)
  grid <- .translation_grid(years, ages)
  catch <- data.frame(year = grid$year, age = grid$age,
                      obs = catch_source$value[match(
                        .translation_key(grid, c("year", "age")),
                        .translation_key(catch_source, c("year", "age")))])

  list(catch = catch,
       units = unique(do.call(rbind, units)),
       provenance = do.call(rbind, provenance))
}

.translation_sampling_time <- function(rows, survey, sampling_times) {
  values <- suppressWarnings(as.numeric(as.character(rows$sampling_time)))
  values <- unique(values[is.finite(values)])
  if (length(values) > 1L) {
    .translation_abort(paste("Conflicting sampling_time values for survey:", survey))
  }
  if (length(values)) return(values[[1]])
  if (!is.null(sampling_times) && survey %in% names(sampling_times)) {
    value <- as.numeric(sampling_times[[survey]])
    if (length(value) == 1L && is.finite(value) && value >= 0 && value <= 1) {
      return(value)
    }
  }
  .translation_abort(paste("Survey timing is unknown for", survey,
                           "; supply a named sampling_times value from 0 to 1."))
}

.translation_weights_for_index <- function(survey, all_weights, stock_weights,
                                            years, ages) {
  matching <- all_weights[!is.na(all_weights$survey) &
                            all_weights$survey == survey, , drop = FALSE]
  if (nrow(matching)) return(matching)
  stock_weights
}

.translation_index <- function(x, all_weights, stock_weights, years, ages,
                               sampling_times, surveys = NULL,
                               exclude_index_years = NULL) {
  supported_measures <- c("numbers_at_age", "proportion_at_age",
                          "total_numbers", "total_biomass", "log_index_sd")
  unsupported <- x[!is.na(x$type) & x$type == "index" &
                     (is.na(x$measure) | !x$measure %in% supported_measures), ,
                   drop = FALSE]
  if (!is.null(surveys)) unsupported <- unsupported[unsupported$survey %in% surveys, , drop = FALSE]
  if (nrow(unsupported)) {
    descriptions <- unique(paste0(unsupported$survey, " (", unsupported$measure, ")"))
    .translation_abort(paste(
      "Unsupported index measure(s) would be omitted from tinyAM obs:",
      paste(descriptions, collapse = "; "),
      ". Select other surveys explicitly or add a documented translation."
    ))
  }

  direct <- .translation_measure(x, "index", "numbers_at_age")
  composition <- x[!is.na(x$type) & x$type == "index" &
                     !is.na(x$measure) & x$measure == "proportion_at_age", , drop = FALSE]
  sampled_age <- direct[grepl("sampled", tolower(direct$unit)), , drop = FALSE]
  direct <- direct[!grepl("sampled", tolower(direct$unit)), , drop = FALSE]
  if (nrow(sampled_age)) composition <- rbind(composition, sampled_age)

  totals <- x[!is.na(x$type) & x$type == "index" &
                !is.na(x$measure) & x$measure %in% c("total_numbers", "total_biomass"), , drop = FALSE]
  known_surveys <- unique(c(as.character(direct$survey), as.character(composition$survey)))
  known_surveys <- known_surveys[!is.na(known_surveys) & nzchar(known_surveys)]
  if (!is.null(surveys)) {
    unknown <- setdiff(surveys, known_surveys)
    if (length(unknown)) .translation_abort(paste("Unknown survey(s):", paste(unknown, collapse = ", ")))
    direct <- direct[direct$survey %in% surveys, , drop = FALSE]
    composition <- composition[composition$survey %in% surveys, , drop = FALSE]
    totals <- totals[totals$survey %in% surveys, , drop = FALSE]
    known_surveys <- surveys
  }
  if (!is.null(exclude_index_years)) {
    if (is.null(names(exclude_index_years)) || any(!nzchar(names(exclude_index_years)))) {
      .translation_abort("exclude_index_years must be a named list of survey years.")
    }
    drop_years <- function(d) {
      for (s in names(exclude_index_years)) {
        if (s %in% d$survey) {
          yr <- as.numeric(exclude_index_years[[s]])
          d <- d[!(d$survey == s & d$year %in% yr), , drop = FALSE]
        }
      }
      d
    }
    direct <- drop_years(direct)
    composition <- drop_years(composition)
    totals <- drop_years(totals)
  }

  result <- list()
  provenance <- list()
  k <- 0L
  if (nrow(direct)) {
    direct <- .translation_calendar_rows(direct, "direct index-at-age")
    if (anyNA(direct$survey) || any(!nzchar(as.character(direct$survey)))) {
      .translation_abort("Every direct index row needs its original survey name.")
    }
    for (survey in unique(as.character(direct$survey))) {
      .translation_one_group(direct[direct$survey == survey, , drop = FALSE],
                             c("sex", "region", "season"),
                             paste("Direct index for", survey))
    }
    unit_counts <- vapply(split(direct$unit, direct$survey), function(z) length(unique(z)), integer(1))
    if (any(unit_counts > 1L)) {
      .translation_abort("A survey has multiple direct index units; convert them before combining rows.")
    }
    direct <- direct[direct$year %in% years & direct$age >= min(ages), , drop = FALSE]
    if (nrow(direct)) {
      direct_source <- direct
      direct_multiplier <- vapply(direct$unit, .translation_index_multiplier, numeric(1))
      direct$value <- direct$value * direct_multiplier
      direct$age <- pmin(direct$age, max(ages))
      direct <- stats::aggregate(value ~ year + age + survey, direct, sum)
      direct$samp_time <- vapply(seq_len(nrow(direct)), function(i) {
        s <- direct$survey[[i]]
        rows <- x[!is.na(x$survey) & x$survey == s & x$year == direct$year[[i]], , drop = FALSE]
        .translation_sampling_time(rows, s, sampling_times)
      }, numeric(1))
      direct$obs <- direct$value
      result[[length(result) + 1L]] <- direct[c("year", "age", "obs", "survey", "samp_time")]
      for (survey in unique(as.character(direct_source$survey))) {
        source <- direct_source[direct_source$survey == survey, , drop = FALSE]
        method <- if (grepl("^(native[ _]+)?survey[ _]?index$",
                            tolower(trimws(source$unit[[1]])))) {
          "native survey index values retained on their source scale"
        } else {
          "source numbers-at-age converted to individual fish"
        }
        k <- k + 1L
        provenance[[k]] <- data.frame(
          survey = survey,
          method = paste(method,
                         "ages above the selected maximum summed into the plus age"),
          source_reference = paste(unique(source$source_reference), collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (nrow(composition)) {
    composition <- .translation_calendar_rows(composition, "index age-composition")
    for (survey in unique(as.character(composition$survey))) {
      comp <- composition[composition$survey == survey &
                            composition$year %in% years, , drop = FALSE]
      if (!nrow(comp)) next
      .translation_one_group(comp, c("sex", "region", "season"),
                             paste("Index composition for", survey))
      totals_survey <- totals[!is.na(totals$survey) & totals$survey == survey &
                                totals$year %in% years, , drop = FALSE]
      if (!nrow(totals_survey)) {
        .translation_abort(paste("Index age composition for", survey,
                                 "has no matching total index values."))
      }
      w <- .translation_weights_for_index(survey, all_weights, stock_weights,
                                           years, ages)
      for (year in unique(comp$year)) {
        cp <- comp[comp$year == year, , drop = FALSE]
        total <- totals_survey[totals_survey$year == year, , drop = FALSE]
        if (!nrow(total)) {
          .translation_abort(paste("Index age composition for", survey, year,
                                   "has no matching total index value."))
        }
        if (nrow(total) != 1L) .translation_abort(paste("Multiple aggregate index values for", survey, year))
        if (anyDuplicated(cp$age)) .translation_abort(paste("Duplicate composition ages for", survey, year))
        if (nrow(cp) > 1L && length(unique(cp$basis)) > 1L) {
          .translation_abort(paste("Mixed composition bases for", survey, year))
        }
        basis <- if (grepl("sampled", tolower(paste(cp$unit, collapse = " ")))) {
          "proportion_numbers"
        } else unique(as.character(cp$basis))
        p <- cp$value
        if (basis == "proportion_numbers" && grepl("sampled", tolower(paste(cp$unit, collapse = " ")))) {
          p <- p / sum(p)
        }
        if (length(basis) != 1L || !basis %in% c("proportion_numbers", "proportion_biomass")) {
          .translation_abort(paste("Unsupported age-composition basis for", survey, year))
        }
        if (any(!is.finite(p)) || any(p < 0) || sum(p) <= 0) {
          .translation_abort(paste("Invalid age composition for", survey, year))
        }
        p <- p / sum(p)
        ww <- w[w$year == year & w$age %in% cp$age, , drop = FALSE]
        if (anyDuplicated(ww$age) || !setequal(ww$age, cp$age)) {
          .translation_abort(paste("Weight-at-age is incomplete for", survey, year,
                                   "age-composition conversion."))
        }
        ww <- ww$value[match(cp$age, ww$age)]
        if (any(!is.finite(ww)) || any(ww <= 0)) {
          .translation_abort(paste("Weights must be positive for", survey, year,
                                   "age-composition conversion."))
        }
        total_value <- total$value[[1]]
        total_measure <- total$measure[[1]]
        if (total_measure == "total_numbers") {
          total_numbers <- total_value * .translation_number_multiplier(total$unit[[1]])
          if (basis == "proportion_numbers") {
            age_numbers <- total_numbers * p
          } else {
            p_numbers <- p / ww
            age_numbers <- total_numbers * p_numbers / sum(p_numbers)
          }
        } else {
          total_kg <- .translation_biomass_to_kg(total_value, total$unit[[1]])
          if (basis == "proportion_numbers") {
            age_numbers <- total_kg * p / sum(p * ww)
          } else {
            age_numbers <- total_kg * p / ww
          }
        }
        translated <- data.frame(year = year, age = cp$age, obs = age_numbers,
                                 survey = survey,
                                 samp_time = .translation_sampling_time(
                                   rbind(cp, total), survey, sampling_times))
        translated <- translated[translated$age >= min(ages), , drop = FALSE]
        translated$age <- pmin(translated$age, max(ages))
        translated <- stats::aggregate(obs ~ year + age + survey + samp_time,
                                       translated, sum)
        k <- k + 1L
        result[[length(result) + 1L]] <- translated
      provenance[[k]] <- data.frame(
        survey = survey,
          method = paste0("age composition x aggregate index reconstructed as numbers-at-age using ",
                          if (any(w$survey == survey, na.rm = TRUE)) {
                            "survey-specific weights"
                          } else {
                            "selected biological weight-at-age series"
                          }, "; ages above the selected maximum summed into the plus age"),
          source_reference = paste(unique(c(cp$source_reference, total$source_reference,
                                             w$source_reference)), collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(result)) .translation_abort("No translatable age-structured index rows were found.")
  out <- do.call(rbind, result)
  out <- out[out$year %in% years & out$age %in% ages, , drop = FALSE]
  key <- .translation_key(out, c("year", "age", "survey"))
  if (anyDuplicated(key)) .translation_abort("Translated index rows are duplicated by year, age, and survey.")
  if (any(!is.finite(out$obs)) || any(out$obs < 0) || anyNA(out$samp_time) ||
      any(out$samp_time < 0 | out$samp_time > 1)) {
    .translation_abort("Translated index values or sampling times are invalid.")
  }
  attr(out, "translation") <- if (length(provenance)) unique(do.call(rbind, provenance)) else data.frame()
  out
}

#' Translate canonical assessment inputs into tinyAM observation tables
#'
#' This analysis-local helper preserves source values where they are already
#' numbers-at-age and makes required transformations explicit. When an
#' annual total number of catch removals and number proportions-at-age are
#' recorded, it multiplies the matching values without renormalizing them.
#' When both an aggregate biomass index and an age composition are available,
#' it derives numbers-at-age using the matching survey weight-at-age when
#' available, or the selected weight series otherwise. Direct indices recorded in a native
#' survey-index scale are preserved without converting them to fish counts.
#' Survey timing must be in the source rows or supplied explicitly; unknown
#' timing is never replaced silently.
#' Index measures without an age-abundance translation are reported as errors
#' unless their surveys are explicitly excluded with `surveys`.
#' The translated index also carries `q_block` (age) and `q_key`
#' (survey-by-age) factors so the fit can state catchability sharing explicitly.
#'
#' @param assessment_id One assessment identifier from `assessments.csv`.
#' @param inputs Canonical `inputs.csv` data.
#' @param years Optional consecutive calendar years to retain. Defaults to the
#'   overlap supported by the selected weight and maturity series.
#' @param ages Optional consecutive model ages. Older catch/index values are
#'   summed into the maximum selected age as a plus group.
#' @param weight_survey Optional survey label selecting the weight series used
#'   for tinyAM biology when multiple weight series overlap.
#' @param sampling_times Optional named numeric vector supplying timing for
#'   surveys whose source rows do not contain a value. Values must be in [0, 1].
#' @param surveys Optional survey names to include. Omitting a source survey is
#'   an analysis choice and should be recorded in the assumption audit.
#' @param exclude_index_years Optional named list of years excluded by the
#'   accepted source assessment, one vector per survey.
#'
#' @return A `tinyAM` observation list with a `translation` attribute recording
#'   catch/index reconstruction methods, unit conversions, and source references.
database_to_tiny_obs <- function(assessment_id, inputs, years = NULL, ages = NULL,
                                 weight_survey = NULL, sampling_times = NULL,
                                 surveys = NULL, exclude_index_years = NULL) {
  x <- .translation_rows(inputs, assessment_id)
  if (!is.null(sampling_times)) {
    if (is.null(names(sampling_times)) || any(!nzchar(names(sampling_times))) ||
        any(!is.finite(as.numeric(sampling_times))) ||
        any(as.numeric(sampling_times) < 0 | as.numeric(sampling_times) > 1)) {
      .translation_abort("sampling_times must be a named numeric vector with values in [0, 1].")
    }
  }

  source_weights <- .translation_calendar_rows(
    .translation_measure(x, "weight", "weight_at_age"), "weight-at-age")
  all_weights <- source_weights
  if (!is.null(weight_survey)) {
    all_weights <- all_weights[!is.na(all_weights$survey) &
                                 all_weights$survey == weight_survey, , drop = FALSE]
    if (!nrow(all_weights)) .translation_abort(paste("No weight series for survey:", weight_survey))
  }
  selected_weight_source <- all_weights
  maturity_source <- .translation_measure(x, "maturity", "maturity_at_age")
  if (any(!is.na(maturity_source$year_basis) & maturity_source$year_basis == "birth_cohort") &&
      !any(!is.na(maturity_source$year_basis) & maturity_source$year_basis == "calendar_year")) {
    .translation_abort("Cohort-specific maturity needs an explicit cohort-to-year mapping before tinyAM use.")
  }
  maturity_source <- .translation_calendar_rows(maturity_source, "maturity-at-age")
  if (any(maturity_source$value > 1)) .translation_abort("Maturity values must be proportions in [0, 1].")

  if (is.null(years)) {
    supported_years <- intersect(unique(all_weights$year), unique(maturity_source$year))
    if (!length(supported_years)) .translation_abort("Weight and calendar-year maturity have no overlapping years.")
    years <- seq.int(min(supported_years), max(supported_years))
  }
  if (is.null(ages)) {
    supported_ages <- intersect(unique(all_weights$age), unique(maturity_source$age))
    if (!length(supported_ages)) .translation_abort("Weight and maturity have no overlapping ages.")
    ages <- seq.int(min(supported_ages), max(supported_ages))
  }
  years <- as.numeric(years)
  ages <- as.numeric(ages)
  if (!length(years) || anyNA(years) || any(!is.finite(years)) ||
      any(years != as.integer(years)) ||
      !isTRUE(all.equal(years, as.numeric(seq.int(min(years), max(years)))))) {
    .translation_abort("years must be a non-empty consecutive sequence of whole years.")
  }
  if (!length(ages) || anyNA(ages) || any(!is.finite(ages)) ||
      any(ages != as.integer(ages)) ||
      !isTRUE(all.equal(ages, as.numeric(seq.int(min(ages), max(ages)))))) {
    .translation_abort("ages must be a non-empty consecutive sequence of whole ages.")
  }

  weight <- .translation_surface(all_weights, years, ages, "Weight-at-age")
  maturity <- .translation_surface(maturity_source, years, ages, "Maturity-at-age")
  catch_translation <- .translation_catch_at_age(x, years, ages)
  catch <- catch_translation$catch

  index <- .translation_index(x, source_weights, selected_weight_source, years, ages,
                              sampling_times, surveys, exclude_index_years)
  index$q_block <- factor(index$age, levels = ages)
  index$q_key <- interaction(index$survey, index$q_block, drop = TRUE,
                             lex.order = TRUE)
  obs <- list(catch = catch, index = index, weight = weight, maturity = maturity)
  index_provenance <- attr(index, "translation")
  source_provenance <- rbind(
    catch_translation$provenance,
    do.call(rbind, lapply(seq_len(nrow(index_provenance)), function(i) {
      survey_rows <- x[!is.na(x$survey) &
                         x$survey == index_provenance$survey[[i]] &
                         x$type == "index", , drop = FALSE]
      data.frame(component = "index", survey = index_provenance$survey[[i]],
                 source_type = paste(unique(survey_rows$source_type), collapse = "; "),
                 method = index_provenance$method[[i]],
                 source_reference = index_provenance$source_reference[[i]],
                 notes = "See translation decisions for survey-specific timing and weight choices.",
                 stringsAsFactors = FALSE)
    })),
    .translation_source_provenance(
      all_weights, "weight",
      paste("Selected weight series:", weight_survey)),
    .translation_source_provenance(
      maturity_source, "maturity",
      "Calendar-year maturity retained on the requested year-age grid.")
  )
  attr(obs, "translation") <- list(
    assessment_id = assessment_id,
    years = years,
    ages = ages,
    weight_survey = weight_survey,
    sampling_times_override = sampling_times,
    excluded_index_years = exclude_index_years,
    catch_method = paste(unique(catch_translation$provenance$method),
                         collapse = "; "),
    catch_units = catch_translation$units,
    source_provenance = source_provenance,
    index = attr(index, "translation")
  )
  attr(obs$index, "translation") <- NULL

  if (requireNamespace("tinyAM", quietly = TRUE)) tinyAM::check_obs(obs)
  obs
}
