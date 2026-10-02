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
database_to_tiny_M <- function(assessment_id, inputs, assumptions = NULL,
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
