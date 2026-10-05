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

.translation_biology_rows <- function(x, label) {
  if (!nrow(x)) .translation_abort(paste("No", label, "inputs are recorded."))
  x$year <- suppressWarnings(as.numeric(as.character(x$year)))
  x$age <- suppressWarnings(as.numeric(as.character(x$age)))
  x$value <- suppressWarnings(as.numeric(as.character(x$value)))
  static <- is.na(x$year)
  if (any(static)) {
    static_rows <- x[static, , drop = FALSE]
    if (any(!is.na(static_rows$year_basis) & nzchar(as.character(static_rows$year_basis))) ||
        anyNA(static_rows$age) || any(!is.finite(static_rows$age)) ||
        any(static_rows$age != as.integer(static_rows$age)) || anyNA(static_rows$value) ||
        any(!is.finite(static_rows$value)) || any(static_rows$value < 0)) {
      .translation_abort(paste(label, "time-invariant rows need finite, non-negative values for whole ages and a blank year basis."))
    }
    if (any(!static)) {
      x <- rbind(.translation_calendar_rows(x[!static, , drop = FALSE], label),
                 static_rows)
    } else {
      x <- static_rows
    }
  } else {
    x <- .translation_calendar_rows(x, label)
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

.translation_repeat_biology <- function(x, years, reference_year, label) {
  static <- x[is.na(x$year), , drop = FALSE]
  if (nrow(static)) {
    annual <- x[!is.na(x$year), , drop = FALSE]
    group_columns <- intersect(c("fleet", "survey", "sex", "region", "season"), names(x))
    group_key <- function(d) .translation_key(d, group_columns)
    static_groups <- unique(group_key(static))
    annual_groups <- unique(group_key(annual))
    if (any(static_groups %in% annual_groups)) {
      .translation_abort(paste(label, "cannot mix time-invariant and annual values within a series."))
    }
    expanded <- lapply(static_groups, function(group) {
      values <- static[group_key(static) == group, , drop = FALSE]
      if (anyDuplicated(values$age)) {
        .translation_abort(paste(label, "time-invariant series need one value per age."))
      }
      do.call(rbind, lapply(years, function(year) {
        values$year <- year
        values$year_basis <- "calendar_year"
        values
      }))
    })
    x <- do.call(rbind, c(list(annual), expanded))
    rownames(x) <- NULL
  }
  if (is.null(reference_year)) return(x)
  if (length(reference_year) != 1L || !is.finite(reference_year) ||
      reference_year != as.integer(reference_year)) {
    .translation_abort(paste(label, "reference year must be one whole year."))
  }
  x <- x[x$year == reference_year, , drop = FALSE]
  if (!nrow(x) || anyDuplicated(x$age)) {
    .translation_abort(paste(label, "reference year needs one source value per age."))
  }
  do.call(rbind, lapply(years, function(year) {
    x$year <- year
    x
  }))
}

.translation_biomass_to_kg <- function(value, unit) {
  unit <- tolower(trimws(as.character(unit)))
  if (grepl("^(thousand|1000|1,000) (t|tonnes?)( |$)", unit)) return(value * 1e6)
  if (grepl("kt", unit, fixed = TRUE)) return(value * 1e6)
  if (grepl("kg", unit, fixed = TRUE)) return(value)
  if (grepl("tonne", unit, fixed = TRUE) || grepl("(^|[^a-z])t($|[^a-z])", unit)) {
    return(value * 1000)
  }
  .translation_abort(paste("Cannot convert biomass unit to kg:", unit))
}

.translation_weight_to_kg <- function(value, unit) {
  unit <- tolower(trimws(as.character(unit)))
  if (grepl("^kg($|[/ _])", unit)) return(value)
  if (grepl("^g($|[/ _])", unit)) return(value / 1000)
  .translation_abort(paste("Cannot convert weight-at-age unit to kg per fish:", unit))
}

.translation_catch_weights <- function(x, composition) {
  weights <- .translation_measure(x, "catch_weight", "weight_at_age")
  columns <- c("year", "fleet", "sex", "region", "season")
  weights <- weights[.translation_key(weights, columns) %in%
                       .translation_key(composition, columns), , drop = FALSE]
  columns <- c(columns, "age")
  if (anyDuplicated(.translation_key(weights, columns))) {
    .translation_abort("Catch weights contain duplicate year-age rows for the catch stream.")
  }
  i <- match(.translation_key(composition, columns), .translation_key(weights, columns))
  if (anyNA(i)) {
    .translation_abort("Biomass-based catch reconstruction needs matching catch weights for every source composition age; stock weights are not substituted.")
  }
  weights <- weights[i, , drop = FALSE]
  weights$value <- mapply(.translation_weight_to_kg, weights$value, weights$unit)
  if (any(!is.finite(weights$value)) || any(weights$value <= 0)) {
    .translation_abort("Catch weights must be positive and finite for biomass reconstruction.")
  }
  weights
}

.translation_number_multiplier <- function(unit) {
  unit <- tolower(trimws(as.character(unit)))
  unit <- gsub("_", " ", unit, fixed = TRUE)
  unit <- gsub("\\s+", " ", unit)
  count_unit <- "(fish|individuals?|numbers?|counts?)"
  if (grepl(paste0("(million|10\\^?6|1,?000,?000)\\s*", count_unit), unit)) return(1e6)
  if (grepl(paste0("(billion|10\\^?9|1,?000,?000,?000)\\s*", count_unit), unit)) return(1e9)
  if (grepl(paste0("(thousand|10\\^?3|1,?000)\\s*", count_unit), unit)) return(1e3)
  if (grepl(count_unit, unit)) return(1)
  .translation_abort(paste("Cannot convert number unit to individual fish:", unit))
}

.translation_index_multiplier <- function(unit) {
  unit <- tolower(trimws(as.character(unit)))
  if (grepl("^(native[ _]+)?survey[ _]?index( \\(unit unresolved\\))?$", unit)) return(1)
  .translation_number_multiplier(unit)
}

.translation_catch_at_age <- function(x, years, ages) {
  direct <- .translation_measure(x, "catch", "numbers_at_age")
  composition <- .translation_measure(x, "catch", "proportion_at_age")
  totals <- x[x$type == "catch" & x$measure %in% c("total_numbers", "total_biomass"), , drop = FALSE]
  if (!nrow(direct) && !nrow(composition)) {
    .translation_abort("No source catch-at-age or number-proportion inputs are recorded.")
  }
  if (nrow(direct)) direct <- .translation_calendar_rows(direct, "catch numbers-at-age")
  if (nrow(composition)) composition <- .translation_calendar_rows(composition, "catch age composition")
  source_rows <- rbind(direct, composition)
  .translation_one_group(source_rows, c("fleet", "sex", "region", "season"),
                         "Catch-at-age input")
  direct <- direct[direct$year %in% years & direct$age >= min(ages), , drop = FALSE]
  composition <- composition[composition$year %in% years, , drop = FALSE]
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
    if (anyNA(basis) || !all(basis %in% c("proportion_numbers", "proportion_biomass"))) {
      .translation_abort("Catch proportions need an explicit number or biomass basis.")
    }
    if (any(composition$value > 1)) {
      .translation_abort("Catch proportions must be in [0, 1].")
    }
    if (!nrow(totals)) {
      .translation_abort("Catch proportions need matching annual total_numbers or total_biomass inputs.")
    }
    totals <- totals[!is.na(totals$year_basis) & totals$year_basis == "calendar_year", , drop = FALSE]
    totals$year <- suppressWarnings(as.numeric(as.character(totals$year)))
    totals$value <- suppressWarnings(as.numeric(as.character(totals$value)))
    if (!nrow(totals) || anyNA(totals$year) || anyNA(totals$value) ||
        any(!is.finite(totals$year)) || any(!is.finite(totals$value)) ||
        any(totals$year != as.integer(totals$year)) || any(totals$value < 0) ||
        any(!is.na(totals$age))) {
      .translation_abort("Catch totals need finite calendar years, non-negative values, and no age dimension.")
    }
    group_columns <- c("year", "fleet", "sex", "region", "season")
    composition_key <- .translation_key(composition, group_columns)
    total_key <- .translation_key(totals, group_columns)
    if (anyDuplicated(total_key)) {
      .translation_abort("Catch totals have multiple representations for a year and catch stream; select one explicitly.")
    }
    total_index <- match(composition_key, total_key)
    if (anyNA(total_index)) {
      .translation_abort("Some catch proportion rows have no matching total_numbers or total_biomass value for the same year and catch stream.")
    }
    composition_age_key <- .translation_key(composition, c(group_columns, "age"))
    if (anyDuplicated(composition_age_key)) {
      .translation_abort("Catch proportions contain duplicate year-age rows before plus-group aggregation.")
    }
    multiplier <- rep(NA_real_, nrow(composition))
    methods <- character()
    weight_sources <- list()
    for (group in unique(composition_key)) {
      rows <- which(composition_key == group)
      cp <- composition[rows, , drop = FALSE]
      total <- totals[total_index[rows[[1]]], , drop = FALSE]
      if (length(unique(cp$basis)) != 1L) {
        .translation_abort("Catch composition mixes number and biomass proportions within a year and stream.")
      }
      p <- cp$value
      if (total$measure == "total_numbers" && cp$basis[[1]] == "proportion_numbers") {
        multiplier[rows] <- .translation_number_multiplier(total$unit)
        numbers <- p * total$value * multiplier[rows]
        method <- "Number proportions multiplied by matching annual total_numbers without renormalizing the source proportions"
      } else {
        weights <- .translation_catch_weights(x, cp)
        w <- weights$value
        weight_sources[[length(weight_sources) + 1L]] <- weights
        if (total$measure == "total_biomass") {
          biomass <- .translation_biomass_to_kg(total$value, total$unit)
          if (cp$basis[[1]] == "proportion_numbers") {
            if (sum(p * w) <= 0) .translation_abort("Catch number proportions need a positive weighted sum.")
            numbers <- biomass * p / sum(p * w)
            method <- "Total catch biomass (kg) times number proportions divided by sum(number proportions times catch weight (kg/fish))"
          } else {
            numbers <- biomass * p / w
            method <- "Total catch biomass (kg) times biomass proportions divided by catch weight (kg/fish), without renormalizing proportions"
          }
        } else {
          multiplier[rows] <- .translation_number_multiplier(total$unit)
          if (sum(p / w) <= 0) .translation_abort("Catch biomass proportions need a positive weight-adjusted sum.")
          numbers <- total$value * multiplier[rows] * (p / w) / sum(p / w)
          method <- "Total catch numbers times (biomass proportions / catch weight) normalized over all source composition ages"
        }
      }
      composition$value[rows] <- numbers
      methods <- c(methods, method)
    }
    converted[[length(converted) + 1L]] <- composition[c("year", "age", "value")]
    units[[length(units) + 1L]] <- data.frame(
      source_unit = totals$unit[total_index], multiplier_to_fish = multiplier,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
    source <- rbind(composition, totals[unique(total_index), , drop = FALSE])
    if (length(weight_sources)) source <- rbind(source, do.call(rbind, weight_sources))
    provenance[[k]] <- .translation_source_provenance(
      source, "catch",
      paste(paste(unique(methods), collapse = "; "),
            "Reconstruction uses all source composition ages before selecting model ages and summing older ages into the plus age."))
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
  catch_source <- catch_source[catch_source$age >= min(ages), , drop = FALSE]
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
                                            years, ages,
                                            index_weight_source = "source") {
  matching <- all_weights[!is.na(all_weights$survey) &
                            all_weights$survey == survey, , drop = FALSE]
  if (index_weight_source == "source" && nrow(matching)) return(matching)
  stock_weights
}

.translation_index <- function(x, all_weights, stock_weights, years, ages,
                               sampling_times, surveys = NULL,
                               exclude_index_years = NULL,
                               index_weight_source = "source") {
  supported_measures <- c("numbers_at_age", "proportion_at_age",
                          "total_numbers", "total_biomass", "log_index_sd",
                          "index_sd", "relative_precision_weight",
                          "relative_standard_error")
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
        method <- if (grepl("^(native[ _]+)?survey[ _]?index( \\(unit unresolved\\))?$",
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
      w <- .translation_weights_for_index(
        survey, all_weights, stock_weights, years, ages,
        index_weight_source = index_weight_source
      )
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
                          } else if (any(!is.na(w$survey) & nzchar(w$survey))) {
                            "selected survey weights as an approximation"
                          } else {
                            "stock weights as an approximation"
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

.translation_index_precision <- function(index, x) {
  for (measure in c("log_index_sd", "relative_precision_weight")) {
    rows <- .translation_measure(x, "index", measure)
    if (!nrow(rows)) next
    rows <- rows[!is.na(rows$year_basis) & rows$year_basis == "calendar_year", , drop = FALSE]
    rows$year <- suppressWarnings(as.numeric(as.character(rows$year)))
    rows$age <- suppressWarnings(as.numeric(as.character(rows$age)))
    rows$value <- suppressWarnings(as.numeric(as.character(rows$value)))
    if (anyNA(rows$year) || any(!is.finite(rows$year)) ||
        any(rows$year != as.integer(rows$year)) || anyNA(rows$value) ||
        any(!is.finite(rows$value)) || any(rows$value <= 0) ||
        anyNA(rows$survey) || any(!nzchar(as.character(rows$survey)))) {
      .translation_abort(paste(measure, "must contain positive finite values for calendar-year survey observations."))
    }
    if (anyNA(rows$age) && !all(is.na(rows$age))) {
      .translation_abort(paste(measure, "cannot mix age-specific and aggregate rows."))
    }
    age_specific <- !all(is.na(rows$age))
    if (age_specific && any(rows$age != as.integer(rows$age))) {
      .translation_abort(paste(measure, "ages must be whole numbers."))
    }
    key_columns <- c("year", if (age_specific) "age", "survey")
    key <- .translation_key(rows, key_columns)
    if (anyDuplicated(key)) .translation_abort(paste("Duplicate", measure, "rows."))
    index_key <- .translation_key(index, key_columns)
    column <- if (measure == "log_index_sd") "relative_sd" else measure
    index[[column]] <- rows$value[match(index_key, key)]
  }
  index
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
#'   for tinyAM biology. By default, an unlabeled stock-weight series is used
#'   when available; otherwise a sole available series is used, and multiple
#'   survey-only series require an explicit choice. Use `""` to select stock
#'   weights explicitly.
#' @param index_weight_source Whether age-composition indices use their matching
#'   survey weight series (`"source"`, the default) or the selected stock-weight
#'   series for every survey (`"stock"`).
#' @param sampling_times Optional named numeric vector supplying timing for
#'   surveys whose source rows do not contain a value. Values must be in [0, 1].
#' @param surveys Optional survey names to include. Omitting a source survey is
#'   an analysis choice and should be recorded in the assumption audit.
#' @param exclude_index_years Optional named list of years excluded by the
#'   accepted source assessment, one vector per survey.
#' @param weight_reference_year,maturity_reference_year Optional source year
#'   whose age vector is explicitly assumed constant over the requested years.
#'   Use only when this constant-vector treatment is documented; these arguments
#'   do not infer constancy from an incomplete annual series.
#' @param maturity_multiplier Optional finite non-negative scalar applied to
#'   translated maturity, for a documented source convention such as a female
#'   fraction.
#' @param assumptions Optional canonical assumptions.csv data passed to
#'   database_to_tam_M() to describe the source M structure.
#'
#' @return A `tinyAM` observation list with a `translation` attribute recording
#'   catch/index reconstruction methods, unit conversions, and source references.
database_to_tam_obs <- function(assessment_id, inputs, years = NULL, ages = NULL,
                                 weight_survey = NULL,
                                 index_weight_source = c("source", "stock"),
                                 sampling_times = NULL,
                                 surveys = NULL, exclude_index_years = NULL,
                                 weight_reference_year = NULL,
                                 maturity_reference_year = NULL,
                                 maturity_multiplier = 1, assumptions = NULL) {
  x <- .translation_rows(inputs, assessment_id)
  index_weight_source <- match.arg(index_weight_source)
  if (!is.null(sampling_times)) {
    if (is.null(names(sampling_times)) || any(!nzchar(names(sampling_times))) ||
        any(!is.finite(as.numeric(sampling_times))) ||
        any(as.numeric(sampling_times) < 0 | as.numeric(sampling_times) > 1)) {
      .translation_abort("sampling_times must be a named numeric vector with values in [0, 1].")
    }
  }
  if (length(maturity_multiplier) != 1L || !is.finite(maturity_multiplier) ||
      maturity_multiplier < 0) {
    .translation_abort("maturity_multiplier must be one finite, non-negative number.")
  }

  source_weights <- .translation_biology_rows(
    .translation_measure(x, "weight", "weight_at_age"), "weight-at-age")
  source_weights$value <- mapply(.translation_weight_to_kg,
                                source_weights$value, source_weights$unit)
  weight_series <- as.character(source_weights$survey)
  weight_series[is.na(weight_series) | !nzchar(weight_series)] <- ""
  if (is.null(weight_survey)) {
    stock_rows <- weight_series == ""
    available_series <- unique(weight_series[!stock_rows])
    if (any(stock_rows)) {
      selected_weight_source <- source_weights[stock_rows, , drop = FALSE]
    } else if (length(available_series) == 1L) {
      selected_weight_source <- source_weights
    } else {
      .translation_abort(paste(
        "Multiple survey-specific weight series are available; select one with weight_survey."
      ))
    }
  } else {
    if (length(weight_survey) != 1L || is.na(weight_survey)) {
      .translation_abort("weight_survey must be NULL or one survey label; use an empty string for stock weights.")
    }
    selected_rows <- if (!nzchar(weight_survey)) {
      weight_series == ""
    } else {
      weight_series == weight_survey
    }
    selected_weight_source <- source_weights[selected_rows, , drop = FALSE]
    if (!nrow(selected_weight_source)) {
      label <- if (!nzchar(weight_survey)) "unlabeled stock" else weight_survey
      .translation_abort(paste("No weight series for:", label))
    }
  }
  selected_weight_label <- if (all(is.na(selected_weight_source$survey) |
                                   !nzchar(as.character(selected_weight_source$survey)))) {
    "unlabeled stock weights"
  } else paste(unique(as.character(selected_weight_source$survey)), collapse = "; ")
  maturity_source <- .translation_measure(x, "maturity", "maturity_at_age")
  if (any(!is.na(maturity_source$year_basis) & maturity_source$year_basis == "birth_cohort") &&
      !any(!is.na(maturity_source$year_basis) & maturity_source$year_basis == "calendar_year")) {
    .translation_abort("Cohort-specific maturity needs an explicit cohort-to-year mapping before tinyAM use.")
  }
  maturity_source <- .translation_biology_rows(maturity_source, "maturity-at-age")
  if (any(maturity_source$value > 1)) .translation_abort("Maturity values must be proportions in [0, 1].")

  if (is.null(years)) {
    weight_years <- unique(selected_weight_source$year[!is.na(selected_weight_source$year)])
    maturity_years <- unique(maturity_source$year[!is.na(maturity_source$year)])
    supported_years <- if (length(weight_years) && length(maturity_years)) {
      intersect(weight_years, maturity_years)
    } else {
      c(weight_years, maturity_years)
    }
    if (!length(supported_years)) .translation_abort("Weight and calendar-year maturity have no overlapping years.")
    years <- seq.int(min(supported_years), max(supported_years))
  }
  if (is.null(ages)) {
    supported_ages <- intersect(unique(selected_weight_source$age), unique(maturity_source$age))
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

  selected_weight_source <- .translation_repeat_biology(
    selected_weight_source, years, weight_reference_year, "Weight-at-age"
  )
  index_weight_source_rows <- .translation_repeat_biology(
    source_weights, years, NULL, "Index weight-at-age"
  )
  weight <- .translation_surface(selected_weight_source, years, ages, "Weight-at-age")
  m_translation <- database_to_tam_M(assessment_id, inputs, assumptions,
                                     years = years, ages = ages)
  if (!is.null(m_translation$surface)) {
    m_key <- .translation_key(weight, c("year", "age"))
    m_surface_key <- .translation_key(m_translation$surface, c("year", "age"))
    m_index <- match(m_key, m_surface_key)
    if (anyNA(m_index)) .translation_abort("Numerical M does not match the translated year-age grid.")
    weight$M_assumption <- m_translation$surface$M_assumption[m_index]
  }
  maturity <- .translation_surface(.translation_repeat_biology(
    maturity_source, years, maturity_reference_year, "Maturity-at-age"),
    years, ages, "Maturity-at-age")
  maturity$obs <- maturity$obs * maturity_multiplier
  if (any(maturity$obs > 1)) {
    .translation_abort("Scaled maturity values must remain in proportions between 0 and 1.")
  }
  catch_translation <- .translation_catch_at_age(x, years, ages)
  catch <- catch_translation$catch

  index <- .translation_index(x, index_weight_source_rows, selected_weight_source, years, ages,
                              sampling_times, surveys, exclude_index_years,
                              index_weight_source)
  index <- .translation_index_precision(index, x)
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
      selected_weight_source, "weight",
      paste("Weight converted to kg per fish. Selected weight series:", selected_weight_label,
            if (!is.null(weight_reference_year)) {
              paste("Constant source vector from", weight_reference_year, "expanded across model years.")
            } else if (all(is.na(selected_weight_source$year))) {
              "Time-invariant source vector expanded across model years."
            } else "Annual source values retained.")),
    .translation_source_provenance(
      maturity_source, "maturity",
      paste(
        if (is.null(maturity_reference_year)) {
          if (all(is.na(maturity_source$year))) {
            "Time-invariant source maturity vector expanded across model years."
          } else {
            "Calendar-year maturity retained on the requested year-age grid."
          }
        } else paste("Constant source maturity vector from", maturity_reference_year,
                     "expanded across model years."),
        if (maturity_multiplier != 1) {
          paste("Translated maturity multiplied by", maturity_multiplier,
                "to match the source SSB convention.")
        } else "No additional maturity scaling."
      ))
  )
  if (!is.null(m_translation$surface)) {
    source_provenance <- rbind(
      source_provenance,
      .translation_source_provenance(
        .translation_measure(x, "M", "natural_mortality_at_age"),
        "M", m_translation$notes
      )
    )
  }
  attr(obs, "translation") <- list(
    assessment_id = assessment_id,
    years = years,
    ages = ages,
    weight_survey = weight_survey,
    index_weight_source = index_weight_source,
    weight_reference_year = weight_reference_year,
    maturity_reference_year = maturity_reference_year,
    maturity_multiplier = maturity_multiplier,
    sampling_times_override = sampling_times,
    selected_surveys = if (is.null(surveys)) {
      unique(na.omit(x$survey[x$type == "index"]))
    } else surveys,
    excluded_surveys = if (is.null(surveys)) character() else {
      setdiff(unique(na.omit(x$survey[x$type == "index"])), surveys)
    },
    excluded_index_years = exclude_index_years,
    catch_method = paste(unique(catch_translation$provenance$method),
                         collapse = "; "),
    catch_units = catch_translation$units,
    source_provenance = source_provenance,
    M = list(status = m_translation$status, notes = m_translation$notes,
             settings = m_translation$M_settings),
    index = attr(index, "translation")
  )
  attr(obs$index, "translation") <- NULL

  if (requireNamespace("tinyAM", quietly = TRUE)) tinyAM::check_obs(obs)
  obs
}


#' Translate a numerical natural-mortality input for tinyAM
#'
#' A fixed numerical M surface can be supplied to tinyAM after joining the
#' returned `surface` to `obs$weight`. When the accepted assessment estimates M,
#' this helper reports that structure but does not substitute fitted M outputs
#' or choose a simplified M process.
#'
#' @param assessment_id One assessment identifier from `assessments.csv`.
#' @param inputs Canonical `inputs.csv` data.
#' @param assumptions Optional canonical `assumptions.csv` data, used only to
#'   describe an estimated or otherwise unavailable source M structure.
#' @param years Calendar years required for a time-invariant M input.
#' @param ages Model ages required for a time-invariant M input.
#'
#' @return A list with `status`, a `surface` data frame when numerical source M
#'   is available, a tinyAM `M_settings` template for fixed M, and explanatory
#'   `notes`.
database_to_tam_M <- function(assessment_id, inputs, assumptions = NULL,
                               years = NULL, ages = NULL) {
  x <- .translation_rows(inputs, assessment_id)
  m <- .translation_measure(x, "M", "natural_mortality_at_age")
  if (!nrow(m)) {
    source_m <- ""
    if (!is.null(assumptions) && all(c("assessment_id", "component", "setting", "value") %in% names(assumptions))) {
      a <- assumptions[!is.na(assumptions$assessment_id) &
                         assumptions$assessment_id == assessment_id &
                         assumptions$component == "M", , drop = FALSE]
      source_m <- paste(a$setting, a$value, collapse = "; ")
    }
    status <- if (grepl("estimate|random walk|time series", source_m, ignore.case = TRUE)) {
      "estimated_in_source"
    } else {
      "no_numerical_input"
    }
    return(list(
      assessment_id = assessment_id,
      status = status,
      surface = NULL,
      M_settings = NULL,
      notes = paste(
        "No numerical fixed M input is recorded. Source-estimated M outputs are not used as inputs.",
        if (nzchar(source_m)) paste("Recorded source structure:", source_m) else ""
      )
    ))
  }

  m$value <- suppressWarnings(as.numeric(as.character(m$value)))
  m$year <- suppressWarnings(as.numeric(as.character(m$year)))
  m$age <- suppressWarnings(as.numeric(as.character(m$age)))
  if (anyNA(m$value) || any(!is.finite(m$value)) || any(m$value <= 0)) {
    .translation_abort("Numerical M inputs must be positive and finite.")
  }
  if (any(!is.na(m$year) & m$year != as.integer(m$year)) ||
      any(!is.na(m$age) & m$age != as.integer(m$age))) {
    .translation_abort("M input years and ages must be whole numbers.")
  }
  if (length(unique(m$unit)) != 1L ||
      !grepl("per.?year|yr.?-?1|year.?-?1", unique(m$unit), ignore.case = TRUE)) {
    .translation_abort("M inputs must use one explicit per-year unit.")
  }

  if (is.null(years)) {
    available_years <- unique(m$year[!is.na(m$year)])
    if (!length(available_years)) {
      .translation_abort("Supply years to expand a time-invariant age-specific M vector.")
    }
    years <- seq.int(min(available_years), max(available_years))
  }
  if (is.null(ages)) {
    available_ages <- unique(m$age[!is.na(m$age)])
    if (!length(available_ages)) {
      .translation_abort("Supply ages to expand a time-invariant M input.")
    }
    ages <- seq.int(min(available_ages), max(available_ages))
  }
  years <- as.numeric(years)
  ages <- as.numeric(ages)
  if (!length(years) || anyNA(years) || any(years != as.integer(years)) ||
      !isTRUE(all.equal(years, as.numeric(seq.int(min(years), max(years)))))) {
    .translation_abort("years must be consecutive whole years.")
  }
  if (!length(ages) || anyNA(ages) || any(ages != as.integer(ages)) ||
      !isTRUE(all.equal(ages, as.numeric(seq.int(min(ages), max(ages)))))) {
    .translation_abort("ages must be consecutive whole ages.")
  }

  has_year <- !all(is.na(m$year))
  has_age <- !all(is.na(m$age))
  if (!has_year && !has_age) {
    .translation_abort("M input needs a year and/or age dimension.")
  }
  if (has_year && anyNA(m$year)) .translation_abort("M input mixes time-varying and time-invariant rows.")
  if (has_age && anyNA(m$age)) .translation_abort("M input mixes age-specific and age-invariant rows.")

  if (has_year && has_age) {
    if (anyDuplicated(.translation_key(m, c("year", "age")))) {
      .translation_abort("M inputs contain duplicate year-age rows.")
    }
    surface <- .translation_surface(m, years, ages, "M")
    names(surface)[[3]] <- "M_assumption"
    notes <- "Source numerical M surface retained on the requested year-age grid."
  } else if (has_age) {
    if (anyDuplicated(m$age) || !setequal(m$age, ages)) {
      .translation_abort("Time-invariant age-specific M must contain one value for every requested age.")
    }
    grid <- .translation_grid(years, ages)
    value <- m$value[match(grid$age, m$age)]
    surface <- data.frame(year = grid$year, age = grid$age, M_assumption = value)
    notes <- "Time-invariant age-specific M vector expanded over requested years."
  } else {
    if (anyDuplicated(m$year) || !setequal(m$year, years)) {
      .translation_abort("Age-invariant M must contain one value for every requested year.")
    }
    grid <- .translation_grid(years, ages)
    value <- m$value[match(grid$year, m$year)]
    surface <- data.frame(year = grid$year, age = grid$age, M_assumption = value)
    notes <- "Age-invariant M time series expanded across requested ages."
  }

  list(
    assessment_id = assessment_id,
    status = "fixed_numerical_input",
    surface = surface,
    M_settings = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption),
    notes = notes
  )
}
