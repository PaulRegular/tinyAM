root <- file.path("analysis", "comp_assessments")
read_table <- function(name) {
  read.csv(file.path(root, "database", name), stringsAsFactors = FALSE,
           na.strings = c("", "NA"), check.names = FALSE)
}

assessments <- read_table("assessments.csv")
assessments <- assessments[assessments$is_current & assessments$is_applied, , drop = FALSE]
assumptions <- read_table("assumptions.csv")
inputs <- read_table("inputs.csv")
outputs <- read_table("outputs.csv")

setting <- function(x, component, name, default = NA_character_) {
  values <- unique(x$value[x$component == component & x$setting == name])
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) values[[1]] else default
}

model_ages <- function(text) {
  if (length(text) != 1L || is.na(text) || !nzchar(text)) return(integer())
  ages <- as.integer(unlist(regmatches(text, gregexpr("[0-9]+", text))))
  if (length(ages) >= 2L) seq(min(ages), max(ages)) else integer()
}

year_bounds <- function(text) {
  if (length(text) != 1L || is.na(text) || !nzchar(text)) return(integer())
  years <- as.integer(unlist(regmatches(text, gregexpr("[0-9]{4}", text))))
  if (length(years) >= 2L) range(years) else integer()
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
  sys.source(file.path(root, "R", "database_to_tam_obs.R"), envir = env)
  env$database_to_tam_obs
}
readiness <- lapply(seq_len(nrow(assessments)), function(i) {
  id <- assessments$assessment_id[i]
  a <- assumptions[assumptions$assessment_id == id, , drop = FALSE]
  x <- inputs[inputs$assessment_id == id, , drop = FALSE]
  y <- outputs[outputs$assessment_id == id, , drop = FALSE]
  year_start <- suppressWarnings(as.integer(setting(a, "model", "time_series_start_year")))
  year_end <- suppressWarnings(as.integer(setting(a, "model", "data_terminal_year")))
  year_text <- setting(a, "population", "modeled_years")
  if (is.na(year_start) || is.na(year_end)) {
    bounds <- year_bounds(year_text)
    if (length(bounds) == 2L) {
      if (is.na(year_start)) year_start <- bounds[[1]]
      if (is.na(year_end)) year_end <- bounds[[2]]
    }
  }
  if (is.na(year_start)) year_start <- suppressWarnings(as.integer(min(x$year, na.rm = TRUE)))
  if (is.na(year_end)) year_end <- assessments$terminal_year[i]
  age_text <- setting(a, "population", "age_range",
                      setting(a, "population", "modeled_ages",
                              setting(a, "population", "ages")))
  ages <- model_ages(age_text)
  years <- if (!is.na(year_start) && !is.na(year_end)) seq(year_start, year_end) else integer()

  catch_direct <- x[x$type == "catch" & x$measure == "numbers_at_age", , drop = FALSE]
  catch_proportion <- x[x$type == "catch" & x$measure == "proportion_at_age", , drop = FALSE]
  catch_age <- rbind(catch_direct, catch_proportion)
  if (!length(ages) && nrow(catch_age) &&
      all(is.finite(catch_age$age)) && !anyNA(catch_age$age)) {
    ages <- seq.int(min(catch_age$age), max(catch_age$age))
  }
  index <- x[x$type == "index", , drop = FALSE]
  index_age_rows <- index[index$measure %in% c("numbers_at_age", "biomass_at_age",
                                                "proportion_at_age"), , drop = FALSE]
  unsupported_index <- index[!index$measure %in% c("numbers_at_age", "biomass_at_age",
                                                    "proportion_at_age", "total_numbers",
                                                    "total_biomass", "log_index_sd",
                                                    "relative_precision_weight"), , drop = FALSE]
  weight <- x[x$type == "weight" & x$measure == "weight_at_age", , drop = FALSE]
  weight_series <- as.character(weight$survey)
  stock_weight <- weight[is.na(weight_series) | !nzchar(weight_series), , drop = FALSE]
  survey_weight_series <- unique(weight_series[!is.na(weight_series) & nzchar(weight_series)])
  model_weight <- if (nrow(stock_weight)) stock_weight else if (length(survey_weight_series) == 1L) weight else weight[0, , drop = FALSE]
  catch_weight <- x[x$type == "catch_weight" & x$measure == "weight_at_age", , drop = FALSE]
  maturity_cohort <- x[!is.na(x$type) & x$type == "maturity" & !is.na(x$year_basis) & x$year_basis == "birth_cohort", , drop = FALSE]
  maturity <- x[!is.na(x$type) & x$type == "maturity" & !is.na(x$measure) &
                  x$measure == "maturity_at_age" & !is.na(x$year_basis) & x$year_basis == "calendar_year", , drop = FALSE]
  maturity_static <- x[!is.na(x$type) & x$type == "maturity" & !is.na(x$measure) &
                    x$measure == "maturity_at_age" & is.na(x$year) & is.na(x$year_basis), , drop = FALSE]
  m_input <- x[x$type == "M", , drop = FALSE]
  m_output <- y[y$measure == "natural_mortality_at_age" & !is.na(y$age), , drop = FALSE]
  documented_surveys <- a$survey[!is.na(a$survey) & nzchar(a$survey)]
  index_surveys <- unique(c(index$survey[!is.na(index$survey) & nzchar(index$survey)],
                            documented_surveys))
  input_times <- if (nrow(index)) {
    paste0(index$survey, "=", ifelse(is.na(index$sampling_time), "unknown", index$sampling_time))
  } else character()
  timing_rows <- a[a$setting %in% c("sampling_time", "rv_sampling_time"), , drop = FALSE]
  input_survey_names <- unique(index$survey[!is.na(index$survey) & nzchar(index$survey)])
  timing_rows <- timing_rows[!timing_rows$survey %in% input_survey_names, , drop = FALSE]
  assumed_times <- if (nrow(timing_rows)) {
    paste0(timing_rows$survey, "=", timing_rows$value)
  } else character()
  index_times <- unique(c(input_times, assumed_times))
  sampling_time <- suppressWarnings(as.numeric(index$sampling_time))
  timing_recorded <- nrow(index) > 0L && all(is.finite(sampling_time) & sampling_time >= 0 & sampling_time <= 1)
  timing_exact <- timing_recorded && !any(grepl("approx", index$notes, ignore.case = TRUE))
  index_age <- nrow(index_age_rows) > 0L && all(!is.na(index_age_rows$age))
  catch_grid <- grid_complete(catch_age, years, ages)
  catch_age_grid <- grid_complete(catch_age, years, ages)
  catch_age_source_grid <- if (nrow(catch_age)) {
    grid_complete(catch_age, seq(min(catch_age$year), max(catch_age$year)),
                  seq(min(catch_age$age), max(catch_age$age)))
  } else FALSE
  weight_grid <- grid_complete(model_weight, years, ages)
  catch_weight_grid <- grid_complete(catch_weight, years, ages)
  maturity_static_grid <- nrow(maturity_static) > 0L && !anyNA(maturity_static$age) && !anyDuplicated(maturity_static$age) && setequal(maturity_static$age, ages)
  maturity_grid <- grid_complete(maturity, years, ages) || maturity_static_grid
  maturity_cohort_source_grid <- if (nrow(maturity_cohort)) {
    grid_complete(maturity_cohort,
                  seq(min(maturity_cohort$year), max(maturity_cohort$year)),
                  seq(min(maturity_cohort$age), max(maturity_cohort$age)))
  } else FALSE
  m_grid <- grid_complete(m_input, years, ages)
  m_estimated <- any(grepl("estimate|random walk|time series",
                           setting(a, "M", "process", ""), ignore.case = TRUE)) ||
    any(grepl("estimate|random walk|time series",
              a$value[a$component == "M"], ignore.case = TRUE))
  m_represented <- m_grid || (nrow(m_input) > 0L && anyNA(m_input$year)) || m_estimated
  has_expected_obs <- all(c(nrow(catch_age) > 0L, nrow(index_age_rows) > 0L,
                            nrow(model_weight) > 0L, nrow(maturity) + nrow(maturity_static) > 0L))
  conversion_ok <- FALSE
  check_obs_ok <- FALSE
  conversion_error <- "Required catch, index, weight, and maturity rows are not all present."
  if (has_expected_obs) {
    conversion <- tryCatch(source_converter()(id, inputs, years = years, ages = ages,
      assumptions = assumptions), error = identity)
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
  missing <- c(
    if (!nrow(catch_age)) "numerical catch-at-age inputs",
    if (nrow(unsupported_index)) paste0(
      "index measures need a documented tinyAM mapping: ",
      paste(unique(paste0(unsupported_index$survey, " (",
                          unsupported_index$measure, ")")), collapse = "; ")
    ),
    if (!index_age) "age-structured survey indices",
    if (!length(index_surveys)) "survey identities" else if (!nrow(index)) "numerical observations from the documented survey series are not transcribed" else if (length(index_surveys) < 2L && grepl("survey|index", setting(a, "data", "assessment_inputs", ""), ignore.case = TRUE)) "other model survey series are not transcribed",
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
    if (!m_represented) "a documented M treatment: supplied numerical values or a source-estimated M structure",
    if (!conversion_ok) paste("database_to_tam_obs failed:", conversion_error),
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
    catch_at_age_rows = nrow(catch_age),
    catch_at_age_source_years = if (nrow(catch_age)) paste(range(catch_age$year), collapse = "-") else "unknown",
    catch_at_age_source_ages = if (nrow(catch_age)) paste(range(catch_age$age), collapse = "-") else "unknown",
    landings_rows = sum(x$type == "catch" &
                          x$measure %in% c("total_numbers", "total_biomass")),
    index_rows = count("index"),
    weight_rows = count("weight"),
    catch_weight_rows = count("catch_weight"),
    maturity_rows = nrow(maturity) + nrow(maturity_static),
    maturity_cohort_rows = nrow(maturity_cohort),
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
    M_estimated_in_source = m_estimated,
    database_to_tam_obs_succeeds = conversion_ok,
    tinyAM_check_obs_passes = check_obs_ok,
    N_at_age_output = any(y$measure == "numbers_at_age" & !is.na(y$age)),
    F_at_age_output = any(y$measure == "fishing_mortality_at_age" & !is.na(y$age)),
    SSB_output = any(y$measure == "SSB"),
    recruitment_output = any(y$measure == "recruitment"),
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
                    "database_to_tam_obs_succeeds", "tinyAM_check_obs_passes", "missing_items")],
      row.names = FALSE)
