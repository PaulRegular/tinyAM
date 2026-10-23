derive_population_total_numbers <- function(outputs, assessment_ids = NULL) {
  fields <- c(
    "assessment_id", "type", "measure", "fleet", "survey", "sex",
    "region", "season", "year", "age", "age_group", "value", "se",
    "lwr", "upr", "unit", "source_type", "source_reference", "notes"
  )
  required <- setdiff(fields, names(outputs))
  if (length(required)) {
    stop("outputs is missing required columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  clean <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    x
  }

  population <- outputs[
    outputs$type == "population" & outputs$measure == "numbers_at_age",
    fields, drop = FALSE
  ]
  if (!is.null(assessment_ids)) {
    population <- population[population$assessment_id %in% assessment_ids, , drop = FALSE]
  }
  if (!nrow(population)) {
    result <- outputs[FALSE, fields, drop = FALSE]
    attr(result, "skipped") <- character()
    return(result)
  }

  population$year <- clean(population$year)
  population$age <- clean(population$age)
  population$age_group <- clean(population$age_group)
  population$sex <- clean(population$sex)
  if (any(!nzchar(population$year))) {
    stop("Population numbers-at-age rows require a year before totals can be derived.", call. = FALSE)
  }

  base_dims <- c("assessment_id", "fleet", "survey", "region", "season")
  key <- function(data, columns) {
    values <- lapply(data[columns], clean)
    do.call(paste, c(values, sep = "\r"))
  }
  series <- split(population, key(population, base_dims))
  existing <- outputs[
    outputs$type == "population" & outputs$measure == "total_numbers",
    fields, drop = FALSE
  ]
  rows <- list()
  skipped <- character()

  for (x in series) {
    x$age_key <- ifelse(nzchar(x$age_group), paste0("group:", x$age_group), paste0("age:", x$age))
    reason <- NULL
    if (any(!nzchar(x$age_key))) {
      reason <- "missing age labels"
    } else if (any(!is.finite(suppressWarnings(as.numeric(x$value))))) {
      reason <- "non-finite age estimates"
    }
    units <- unique(clean(x$unit))
    if (is.null(reason) && (length(units) != 1L || !nzchar(units))) {
      reason <- "mixed or missing units"
    }

    years <- unique(x$year)
    sexes_by_year <- lapply(years, function(year) sort(unique(x$sex[x$year == year])))
    if (is.null(reason) && !all(vapply(sexes_by_year, identical, logical(1), sexes_by_year[[1L]]))) {
      reason <- "sex coverage changes by year"
    }
    sexes <- sexes_by_year[[1L]]

    if (is.null(reason)) for (sex in sexes) {
      age_sets <- lapply(years, function(year) {
        sort(unique(x$age_key[x$year == year & x$sex == sex]))
      })
      if (!all(vapply(age_sets, identical, logical(1), age_sets[[1L]]))) {
        reason <- paste0("age coverage changes by year for sex '", sex, "'")
        break
      }
      for (year in years) {
        current <- x[x$year == year & x$sex == sex, , drop = FALSE]
        if (anyDuplicated(current$age_key)) {
          reason <- paste0("duplicate age groups for year ", year, " and sex '", sex, "'")
          break
        }
      }
      if (!is.null(reason)) break
    }
    if (!is.null(reason)) {
      assessment_id <- clean(x$assessment_id[[1L]])
      skipped <- c(skipped, paste0(assessment_id, ": ", reason))
      next
    }

    output_sex <- if (length(sexes) > 1L) "combined" else sexes[[1L]]
    dimension_values <- vapply(base_dims, function(name) clean(x[[name]][1L]), character(1))
    for (year in years) {
      match_existing <- existing[
        clean(existing$assessment_id) == dimension_values[["assessment_id"]] &
          clean(existing$fleet) == dimension_values[["fleet"]] &
          clean(existing$survey) == dimension_values[["survey"]] &
          clean(existing$region) == dimension_values[["region"]] &
          clean(existing$season) == dimension_values[["season"]] &
          clean(existing$year) == year &
          clean(existing$sex) %in% if (output_sex == "combined") c("", "combined") else output_sex,
        , drop = FALSE
      ]
      if (nrow(match_existing)) next

      current <- x[x$year == year, , drop = FALSE]
      age_labels <- sort(unique(current$age_key))
      sex_note <- if (length(sexes) > 1L) {
        paste0(" and sexes (", paste(sexes[nzchar(sexes)], collapse = ", "), ")")
      } else {
        ""
      }
      rows[[length(rows) + 1L]] <- data.frame(
        assessment_id = dimension_values[["assessment_id"]],
        type = "population",
        measure = "total_numbers",
        fleet = dimension_values[["fleet"]],
        survey = dimension_values[["survey"]],
        sex = output_sex,
        region = dimension_values[["region"]],
        season = dimension_values[["season"]],
        year = year,
        age = "",
        age_group = "",
        value = sum(as.numeric(current$value)),
        se = NA_real_,
        lwr = NA_real_,
        upr = NA_real_,
        unit = units[[1L]],
        source_type = "derived_source_output",
        source_reference = paste0(
          "Derived from the corresponding official population/numbers_at_age rows; ",
          "see their source references"
        ),
        notes = paste0(
          "Sum across all reported age groups (", paste(age_labels, collapse = ", "), ")",
          sex_note,
          ". No standard error or interval was combined."
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  result <- if (!length(rows)) {
    outputs[FALSE, fields, drop = FALSE]
  } else {
    do.call(rbind, rows)[fields]
  }
  attr(result, "skipped") <- skipped
  result
}
