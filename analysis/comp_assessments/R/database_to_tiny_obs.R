database_to_tiny_obs <- function(assessment_id, inputs) {
  required <- c("assessment_id", "type", "measure", "basis", "fleet", "survey",
                "sex", "region", "season", "year", "year_basis", "age",
                "value", "unit", "sampling_time", "source_type",
                "source_reference", "transformation", "notes")
  missing <- setdiff(required, names(inputs))
  if (length(missing)) {
    stop("inputs is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(assessment_id) != 1L || is.na(assessment_id) || !nzchar(assessment_id)) {
    stop("assessment_id must be one non-empty value.", call. = FALSE)
  }

  x <- inputs[inputs$assessment_id == assessment_id, , drop = FALSE]
  if (!nrow(x)) stop("No input rows found for assessment_id: ", assessment_id, call. = FALSE)

  key <- function(d, cols) do.call(paste, c(d[cols], sep = "\r"))
  one_value <- function(x, what) {
    values <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
    if (length(values) != 1L) {
      stop("This assessment needs exactly one ", what,
           " to fit tinyAM's current observation tables.", call. = FALSE)
    }
    values
  }
  one_group <- function(x, what) {
    values <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
    if (length(values) > 1L) {
      stop("This assessment has multiple ", what,
           " values that tinyAM's current observation tables cannot keep separate.", call. = FALSE)
    }
  }
  numeric_rows <- function(d, label, allow_missing = FALSE) {
    if (!nrow(d)) stop("No ", label, " inputs are recorded for this assessment.", call. = FALSE)
    if (!is.numeric(d$year) || !is.numeric(d$age) || !is.numeric(d$value) ||
        anyNA(d$year) || anyNA(d$age) || (!allow_missing && anyNA(d$value))) {
      stop(label, " rows need numeric year and age fields", if (!allow_missing) " and non-missing values" else " and numeric values", ".", call. = FALSE)
    }
    if (any(!is.finite(d$year)) || any(!is.finite(d$age)) ||
        any(d$year != as.integer(d$year)) || any(d$age != as.integer(d$age))) {
      stop(label, " years and ages must be finite whole numbers.", call. = FALSE)
    }
    values <- d$value[!is.na(d$value)]
    if (any(!is.finite(values)) || any(values < 0)) {
      stop(label, " values must be finite and non-negative.", call. = FALSE)
    }
    if (anyDuplicated(key(d, c("year", "age")))) {
      stop(label, " has more than one value for a year-age row; no fleets, sexes, or seasons were combined.", call. = FALSE)
    }
    d
  }
  pick <- function(type, measure = NULL) {
    keep <- x$type == type
    if (!is.null(measure)) keep <- keep & x$measure == measure
    x[keep, , drop = FALSE]
  }
  as_obs <- function(d) {
    d$obs <- d$value
    d
  }

  catch_source <- numeric_rows(pick("catch", "numbers_at_age"),
                               "catch numbers-at-age", allow_missing = TRUE)
  if (nrow(catch_source) && length(unique(as.character(catch_source$fleet))) != 1L) {
    stop("tinyAM's catch table accepts one fleet; catch fleets were not combined.", call. = FALSE)
  }
  for (field in c("sex", "region", "season")) one_group(catch_source[[field]], paste("catch ", field, " group", sep = ""))

  index_source <- numeric_rows(pick("index", "numbers_at_age"),
                               "index numbers-at-age", allow_missing = TRUE)
  if (anyNA(index_source$survey) || any(!nzchar(as.character(index_source$survey)))) {
    stop("Every index row needs its original survey name.", call. = FALSE)
  }
  if (!is.numeric(index_source$sampling_time) || anyNA(index_source$sampling_time) ||
      any(index_source$sampling_time < 0 | index_source$sampling_time > 1)) {
    stop("Every index row needs a numeric sampling_time from 0 to 1.", call. = FALSE)
  }
  if (anyDuplicated(key(index_source, c("year", "age", "survey")))) {
    stop("Index rows must be unique by year, age, and survey.", call. = FALSE)
  }

  weight_source <- numeric_rows(pick("weight", "weight_at_age"), "weight")
  maturity_source <- x[x$type == "maturity" & x$measure == "maturity_at_age" &
                         x$year_basis == "calendar_year", , drop = FALSE]
  if (!nrow(maturity_source) && any(x$type == "maturity" &
                                     x$year_basis == "birth_cohort")) {
    stop("Cohort-specific maturity needs an explicit cohort-to-year mapping before tinyAM use.",
         call. = FALSE)
  }
  maturity_source <- numeric_rows(maturity_source, "calendar-year maturity")
  if (any(maturity_source$value > 1)) stop("Maturity values must be proportions from 0 to 1.", call. = FALSE)
  if (!setequal(key(weight_source, c("year", "age")),
                key(maturity_source, c("year", "age")))) {
    stop("Weight and maturity inputs must cover the same year-age rows.", call. = FALSE)
  }

  grid <- expand.grid(year = sort(unique(weight_source$year)),
                      age = sort(unique(weight_source$age)))
  if (!setequal(key(grid, c("year", "age")), key(weight_source, c("year", "age")))) {
    stop("Weight inputs must cover every year-age combination in the assessment.", call. = FALSE)
  }
  if (min(index_source$year, na.rm = TRUE) < min(grid$year) ||
      max(index_source$year, na.rm = TRUE) > max(grid$year) ||
      min(catch_source$year, na.rm = TRUE) < min(grid$year) ||
      max(catch_source$year, na.rm = TRUE) > max(grid$year) ||
      min(index_source$age, na.rm = TRUE) < min(grid$age) ||
      max(index_source$age, na.rm = TRUE) > max(grid$age) ||
      min(catch_source$age, na.rm = TRUE) < min(grid$age) ||
      max(catch_source$age, na.rm = TRUE) > max(grid$age)) {
    stop("Catch and survey years and ages must fall within the biological input range.", call. = FALSE)
  }

  catch <- grid
  catch$value <- catch_source$value[match(key(grid, c("year", "age")),
                                          key(catch_source, c("year", "age")))]
  catch$obs <- catch$value
  catch$fleet <- one_value(catch_source$fleet, "catch fleet")
  index_source$samp_time <- index_source$sampling_time
  index_source$sampling_time <- NULL
  index <- as_obs(index_source)
  weight <- as_obs(weight_source)
  maturity <- as_obs(maturity_source)
  obs <- list(catch = catch, index = index, weight = weight, maturity = maturity)

  if (requireNamespace("tinyAM", quietly = TRUE)) tinyAM::check_obs(obs)
  obs
}
