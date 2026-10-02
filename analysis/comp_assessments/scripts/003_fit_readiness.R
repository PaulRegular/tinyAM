root <- file.path("analysis", "comp_assessments")
read_table <- function(name) {
  read.csv(file.path(root, "database", name), stringsAsFactors = FALSE,
           na.strings = c("", "NA"), check.names = FALSE)
}

assessments <- read_table("assessments.csv")
assumptions <- read_table("assumptions.csv")
inputs <- read_table("inputs.csv")
outputs <- read_table("outputs.csv")

setting <- function(x, component, name, default = NA_character_) {
  values <- unique(x$value[x$component == component & x$setting == name])
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) values[[1]] else default
}

model_ages <- function(text) {
  ages <- as.integer(unlist(regmatches(text, gregexpr("[0-9]+", text))))
  if (length(ages) >= 2L) seq(min(ages), max(ages)) else integer()
}

grid_complete <- function(x, years, ages) {
  if (!length(years) || !length(ages) || !nrow(x) || anyNA(x$year) || anyNA(x$age)) {
    return(FALSE)
  }
  expected <- expand.grid(year = years, age = ages)
  key <- function(d) paste(d$year, d$age, sep = "\r")
  !anyDuplicated(key(x)) && setequal(key(x), key(expected))
}

source_converter <- function() {
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "R", "database_to_tiny_obs.R"), envir = env)
  env$database_to_tiny_obs
}

readiness <- lapply(seq_len(nrow(assessments)), function(i) {
  id <- assessments$assessment_id[i]
  a <- assumptions[assumptions$assessment_id == id, , drop = FALSE]
  x <- inputs[inputs$assessment_id == id, , drop = FALSE]
  y <- outputs[outputs$assessment_id == id, , drop = FALSE]
  year_start <- suppressWarnings(as.integer(setting(a, "model", "time_series_start_year")))
  if (is.na(year_start)) year_start <- suppressWarnings(as.integer(min(x$year, na.rm = TRUE)))
  year_end <- suppressWarnings(as.integer(setting(a, "model", "data_terminal_year")))
  if (is.na(year_end)) year_end <- assessments$terminal_year[i]
  age_text <- setting(a, "population", "age_range")
  ages <- model_ages(age_text)
  years <- if (!is.na(year_start) && !is.na(year_end)) seq(year_start, year_end) else integer()

  catch <- x[x$type == "catch", , drop = FALSE]
  catch_age <- x[x$type == "catch_at_age", , drop = FALSE]
  index <- x[x$type == "index", , drop = FALSE]
  weight <- x[x$type == "weight", , drop = FALSE]
  catch_weight <- x[x$type == "catch_weight", , drop = FALSE]
  maturity <- x[x$type == "maturity", , drop = FALSE]
  maturity_cohort <- x[x$type == "maturity_cohort", , drop = FALSE]
  m_input <- x[x$type == "M", , drop = FALSE]
  m_output <- y[y$type == "M" & !is.na(y$age), , drop = FALSE]
  index_surveys <- unique(index$survey[!is.na(index$survey) & nzchar(index$survey)])
  index_times <- if (nrow(index)) unique(paste0(index$survey, "=", ifelse(is.na(index$samp_time), "unknown", index$samp_time))) else character()
  samp_time <- suppressWarnings(as.numeric(index$samp_time))
  timing_recorded <- nrow(index) > 0L && all(is.finite(samp_time) & samp_time >= 0 & samp_time <= 1)
  timing_exact <- timing_recorded && !any(grepl("approx", index$notes, ignore.case = TRUE))
  index_age <- nrow(index) > 0L && all(!is.na(index$age))
  catch_grid <- grid_complete(catch, years, ages)
  catch_age_grid <- grid_complete(catch_age, years, ages)
  catch_age_source_grid <- if (nrow(catch_age)) {
    grid_complete(catch_age, seq(min(catch_age$year), max(catch_age$year)),
                  seq(min(catch_age$age), max(catch_age$age)))
  } else FALSE
  weight_grid <- grid_complete(weight, years, ages)
  catch_weight_grid <- grid_complete(catch_weight, years, ages)
  maturity_grid <- grid_complete(maturity, years, ages)
  maturity_cohort_source_grid <- if (nrow(maturity_cohort)) {
    grid_complete(maturity_cohort,
                  seq(min(maturity_cohort$year), max(maturity_cohort$year)),
                  seq(min(maturity_cohort$age), max(maturity_cohort$age)))
  } else FALSE
  m_grid <- grid_complete(m_input, years, ages) || grid_complete(m_output, years, ages)
  has_expected_obs <- all(c("catch", "index", "weight", "maturity") %in% x$type)
  conversion_ok <- FALSE
  check_obs_ok <- FALSE
  conversion_error <- "Required catch, index, weight, and maturity rows are not all present."
  if (has_expected_obs) {
    conversion <- tryCatch(source_converter()(id, inputs), error = identity)
    conversion_ok <- !inherits(conversion, "error")
    if (conversion_ok) {
      conversion_error <- ""
      if (requireNamespace("tinyAM", quietly = TRUE)) {
        check_result <- tryCatch(tinyAM::check_obs(conversion), error = identity)
        check_obs_ok <- !inherits(check_result, "error")
      }
    } else {
      conversion_error <- conditionMessage(conversion)
    }
  }
  output_types <- unique(y$type)
  has_age_output <- function(type) any(y$type == type & !is.na(y$age))
  missing <- c(
    if (!catch_grid) "complete catch-at-age observations on the full modeled year-age grid",
    if (!nrow(catch_age)) "numerical catch-at-age inputs" else if (!catch_age_source_grid) "gaps within the available catch-age year-age series",
    if (!index_age) "age-structured survey indices",
    if (!length(index_surveys)) "survey identities" else if (length(index_surveys) < 2L && grepl("survey|index", setting(a, "data", "assessment_inputs", ""), ignore.case = TRUE)) "other model survey series are not transcribed",
    if (!timing_exact) "exact survey sampling times (recorded times may be season-level approximations)",
    if (!weight_grid) "complete stock weight-at-age matrix for all modeled years and ages",
    if (!maturity_grid) {
      if (nrow(maturity_cohort) && maturity_cohort_source_grid) {
        "cohort-specific maturity needs an explicit cohort-to-year mapping before tinyAM use"
      } else if (nrow(maturity_cohort)) {
        "gaps within the available cohort-specific maturity series"
      } else {
        "complete maturity-at-age matrix for all modeled years and ages"
      }
    },
    if (!m_grid) "numerical M-at-age values over the modeled years",
    if (!conversion_ok) paste("database_to_tiny_obs failed:", conversion_error),
    if (conversion_ok && !requireNamespace("tinyAM", quietly = TRUE)) "tinyAM is not installed; check_obs was not run",
    if (conversion_ok && requireNamespace("tinyAM", quietly = TRUE) && !check_obs_ok) "tinyAM::check_obs failed"
  )
  type_counts <- table(x$type)
  count <- function(type) if (type %in% names(type_counts)) as.integer(type_counts[[type]]) else 0L

  data.frame(
    assessment_id = id,
    stock_id = assessments$stock_id[i],
    assessment_year = assessments$assessment_year[i],
    model = assessments$model_family[i],
    framework_year = assessments$framework_year[i],
    modeled_first_year = year_start,
    modeled_terminal_year = year_end,
    minimum_age = if (length(ages)) min(ages) else NA_integer_,
    maximum_age = if (length(ages)) max(ages) else NA_integer_,
    recruitment_age = setting(a, "N", "recruitment_age", setting(a, "recruitment", "recruitment_age")),
    plus_group = setting(a, "population", "plus_group", "unknown"),
    catch_rows = count("catch"),
    catch_at_age_rows = count("catch_at_age"),
    catch_at_age_source_years = if (nrow(catch_age)) paste(range(catch_age$year), collapse = "-") else "unknown",
    catch_at_age_source_ages = if (nrow(catch_age)) paste(range(catch_age$age), collapse = "-") else "unknown",
    landings_rows = count("landings"),
    index_rows = count("index"),
    weight_rows = count("weight"),
    catch_weight_rows = count("catch_weight"),
    maturity_rows = count("maturity"),
    maturity_cohort_rows = count("maturity_cohort"),
    maturity_cohort_years = if (nrow(maturity_cohort)) paste(range(maturity_cohort$year), collapse = "-") else "unknown",
    maturity_cohort_ages = if (nrow(maturity_cohort)) paste(range(maturity_cohort$age), collapse = "-") else "unknown",
    maturity_cohort_source_grid_complete = maturity_cohort_source_grid,
    M_rows = count("M"),
    surveys = if (length(index_surveys)) paste(index_surveys, collapse = "; ") else "none",
    survey_sampling_times = if (length(index_times)) paste(index_times, collapse = "; ") else "unknown",
    catch_full_year_age_grid = catch_grid,
    numerical_catch_at_age_available = nrow(catch_age) > 0L,
    catch_at_age_source_grid_complete = catch_age_source_grid,
    catch_at_age_full_model_grid = catch_age_grid,
    index_age_structured = index_age,
    survey_timing_known = timing_exact,
    survey_timing_recorded = timing_recorded,
    weight_full_year_age_grid = weight_grid,
    catch_weight_full_year_age_grid = catch_weight_grid,
    maturity_full_year_age_grid = maturity_grid,
    M_numerical_values_available = m_grid,
    database_to_tiny_obs_succeeds = conversion_ok,
    tinyAM_check_obs_passes = check_obs_ok,
    N_at_age_output = has_age_output("N"),
    F_at_age_output = has_age_output("F"),
    SSB_output = "SSB" %in% output_types,
    recruitment_output = "recruitment" %in% output_types,
    missing_items = if (length(missing)) paste(missing, collapse = "; ") else "none",
    stringsAsFactors = FALSE
  )
})
readiness <- do.call(rbind, readiness)
dir.create(file.path(root, "results"), showWarnings = FALSE, recursive = TRUE)
write.csv(readiness, file.path(root, "results", "fit_readiness.csv"), row.names = FALSE, na = "")
print(readiness[, c("assessment_id", "modeled_first_year", "modeled_terminal_year",
                    "minimum_age", "maximum_age", "catch_rows", "catch_at_age_rows",
                    "index_rows", "weight_rows", "catch_weight_rows", "maturity_rows",
                    "maturity_cohort_rows", "M_rows",
                    "database_to_tiny_obs_succeeds", "tinyAM_check_obs_passes", "missing_items")],
      row.names = FALSE)
