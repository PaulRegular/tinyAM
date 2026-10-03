root <- file.path("analysis", "comp_assessments", "database")
read_table <- function(name) {
  read.csv(file.path(root, name), stringsAsFactors = FALSE,
           na.strings = c("", "NA"), check.names = FALSE)
}

stocks <- read_table("stocks.csv")
assessments <- read_table("assessments.csv")
assumptions <- read_table("assumptions.csv")
inputs <- read_table("inputs.csv")
outputs <- read_table("outputs.csv")
if (!"age_group" %in% names(outputs)) outputs$age_group <- ""

schemas <- list(
  stocks = c("stock_id", "charbonneau_id", "authority", "authority_stock_id",
             "scientific_name", "common_name", "area", "region", "ocean", "notes"),
  assessments = c("assessment_id", "stock_id", "assessment_year", "terminal_year",
                  "estimate_terminal_year", "assessment_type", "model_family",
                  "model_version", "is_current", "is_applied", "framework_year",
                  "assessment_url", "framework_url", "data_url", "model_url",
                  "repository_url", "assumptions_status", "inputs_status",
                  "outputs_status", "notes"),
  assumptions = c("assessment_id", "component", "fleet", "survey", "sex", "region",
                  "season", "setting", "value", "source_reference", "notes"),
  inputs = c("assessment_id", "type", "measure", "basis", "fleet", "survey", "sex",
             "region", "season", "year", "year_basis", "age", "value", "unit",
             "sampling_time", "source_type", "source_reference", "transformation", "notes"),
  outputs = c("assessment_id", "type", "measure", "fleet", "survey", "sex", "region",
              "season", "year", "age", "value", "se", "lwr", "upr", "unit",
              "source_type", "source_reference", "notes")
)
tables <- list(stocks = stocks, assessments = assessments, assumptions = assumptions,
               inputs = inputs, outputs = outputs)
for (name in names(tables)) {
  missing <- setdiff(schemas[[name]], names(tables[[name]]))
  if (length(missing)) stop(name, ".csv is missing columns: ",
                            paste(missing, collapse = ", "), call. = FALSE)
}

unique_id <- function(x, column, table) {
  id <- x[[column]]
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id)) {
    stop(table, "$", column, " must contain unique, non-empty values.", call. = FALSE)
  }
}
unique_id(stocks, "stock_id", "stocks")
unique_id(assessments, "assessment_id", "assessments")
if (any(!assessments$stock_id %in% stocks$stock_id)) {
  stop("assessments.csv has an unknown stock_id.", call. = FALSE)
}
assessment_ids <- assessments$assessment_id
for (name in c("assumptions", "inputs", "outputs")) {
  x <- tables[[name]]
  if (anyNA(x$assessment_id) || any(!x$assessment_id %in% assessment_ids)) {
    stop(name, ".csv has an unknown or empty assessment_id.", call. = FALSE)
  }
}

status_values <- c("not_started", "partial", "complete", "not_applicable")
for (field in c("assumptions_status", "inputs_status", "outputs_status")) {
  if (anyNA(assessments[[field]]) || any(!assessments[[field]] %in% status_values)) {
    stop("assessments$", field, " must use the documented status values.", call. = FALSE)
  }
}
logical_value <- function(x) tolower(as.character(x)) %in% c("true", "1")
for (field in c("is_current", "is_applied")) {
  raw <- tolower(as.character(assessments[[field]]))
  if (any(!raw %in% c("true", "false", "1", "0"))) {
    stop("assessments$", field, " must be TRUE/FALSE or 1/0.", call. = FALSE)
  }
}
current <- logical_value(assessments$is_current)
applied <- logical_value(assessments$is_applied)
if (any(current & !applied)) stop("A current assessment must be marked applied.", call. = FALSE)
for (i in seq_len(nrow(assessments))) {
  row <- assessments[i, ]
  for (table_name in c("assumptions", "inputs", "outputs")) {
    status <- row[[paste0(if (table_name == "assumptions") "assumptions" else table_name,
                          "_status")]]
    has_rows <- any(tables[[table_name]]$assessment_id == row$assessment_id)
    if (status == "complete" && !has_rows) {
      stop(row$assessment_id, " is marked complete for ", table_name,
           " but has no rows.", call. = FALSE)
    }
    if (status == "not_applicable" && has_rows) {
      stop(row$assessment_id, " is marked not_applicable for ", table_name,
           " but has rows.", call. = FALSE)
    }
  }
}

whole_year <- function(x, field, table, allow_blank = FALSE) {
  raw <- x[[field]]
  value <- suppressWarnings(as.numeric(raw))
  required <- if (allow_blank) !is.na(raw) & nzchar(as.character(raw)) else rep(TRUE, length(raw))
  if (any(required & (is.na(value) | !is.finite(value) | value != as.integer(value)))) {
    stop(table, "$", field, " must contain whole years or blanks.", call. = FALSE)
  }
  value
}
assessment_year <- whole_year(assessments, "assessment_year", "assessments")
terminal_year <- whole_year(assessments, "terminal_year", "assessments")
estimate_year <- whole_year(assessments, "estimate_terminal_year", "assessments",
                            allow_blank = TRUE)
invisible(whole_year(assessments, "framework_year", "assessments", allow_blank = TRUE))
if (any(terminal_year > assessment_year) || any(estimate_year < terminal_year, na.rm = TRUE)) {
  stop("Assessment years, data terminals, and estimate terminals are inconsistent.",
       call. = FALSE)
}
current_stock <- assessments$stock_id[current]
if (anyDuplicated(current_stock)) {
  repeated <- unique(current_stock[duplicated(current_stock)])
  for (stock in repeated) {
    rows <- assessments[assessments$stock_id == stock & current, , drop = FALSE]
    if (anyNA(rows$notes) || any(!nzchar(rows$notes))) {
      stop("Multiple current assessments for ", stock,
           " need an explanation in notes.", call. = FALSE)
    }
  }
}

duplicates <- list(
  assumptions = c("assessment_id", "component", "fleet", "survey", "sex", "region",
                  "season", "setting"),
  inputs = c("assessment_id", "type", "measure", "basis", "fleet", "survey", "sex",
             "region", "season", "year", "year_basis", "age"),
  outputs = c("assessment_id", "type", "measure", "fleet", "survey", "sex", "region",
              "season", "year", "age", "age_group")
)
composition_dimensions <- c("observation_id", "length_bin", "length_bin_lower",
                            "length_bin_upper", "sample_size", "age_error", "partition")
duplicates$inputs <- c(duplicates$inputs,
                      intersect(composition_dimensions, names(inputs)))
for (name in names(duplicates)) {
  x <- tables[[name]]
  key <- do.call(paste, c(lapply(x[duplicates[[name]]], function(z) {
    z[is.na(z)] <- ""
    as.character(z)
  }), sep = "\r"))
  if (anyDuplicated(key)) stop(name, ".csv has duplicate canonical rows.", call. = FALSE)
}

allowed_sources <- c("native_model", "official_machine_readable", "official_table",
                     "official_document",
                     "digitized", "reconstructed_source_input", "charbonneau_seed")
allowed_input_types <- c("catch", "index", "weight", "catch_weight", "maturity", "M", "covariate", "biology")
allowed_output_types <- c("population", "mortality", "biomass", "recruitment", "catch",
                          "index", "catchability")
allowed_input_measures <- c("numbers_at_age", "biomass_at_age", "total_numbers",
                            "total_biomass", "proportion_at_age", "proportion_at_length",
                            "conditional_proportion_at_age", "weight_at_age", "spawning_weight_at_age",
                            "maturity_at_age", "natural_mortality_at_age",
                            "landings_proportion", "landings_numbers_at_age", "landings_fraction_at_age",
                            "landings_weight_at_age", "discard_weight_at_age", "log_index_sd", "larval_abundance_index", "environmental_covariate", "fraction_F_before_spawning",
                            "fraction_M_before_spawning")
allowed_bases <- c("numbers", "biomass", "proportion_numbers", "proportion_biomass",
                   "kg_per_fish", "proportion", "per_year", "log_scale", "native_covariate")
for (name in c("inputs", "outputs")) {
  x <- tables[[name]]
  if (anyNA(x$source_type) || any(!x$source_type %in% allowed_sources)) {
    stop(name, "$source_type must use the controlled source labels.", call. = FALSE)
  }
  if (anyNA(x$source_reference) || any(!nzchar(x$source_reference))) {
    stop(name, "$source_reference must identify the source location.", call. = FALSE)
  }
  allowed_type <- if (name == "inputs") allowed_input_types else allowed_output_types
  if (anyNA(x$type) || any(!x$type %in% allowed_type)) {
    stop(name, "$type must use the documented broad categories.", call. = FALSE)
  }
  if (anyNA(x$measure) || any(!nzchar(x$measure))) {
    stop(name, "$measure must identify the exact quantity.", call. = FALSE)
  }
  if (name == "inputs") {
    if (any(!x$measure %in% allowed_input_measures) || anyNA(x$basis) ||
        any(!x$basis %in% allowed_bases)) {
      stop("inputs.csv has an undocumented measure or basis.", call. = FALSE)
    }
    if (any(!x$year_basis %in% c("calendar_year", "birth_cohort"))) {
      stop("inputs$year_basis must be calendar_year or birth_cohort.", call. = FALSE)
    }
  }
  year <- whole_year(x, "year", name, allow_blank = name == "outputs")
  blank_year <- is.na(x$year) | !nzchar(as.character(x$year))
  if (name == "outputs" && any(blank_year & !x$measure %in% c("q", "q_power"))) {
    stop("outputs$year may be blank only for time-invariant q or q_power estimates.", call. = FALSE)
  }
  age_required <- if (name == "inputs") {
    x$measure %in% c("numbers_at_age", "biomass_at_age", "proportion_at_age",
                     "conditional_proportion_at_age", "fraction_F_before_spawning",
                     "fraction_M_before_spawning", "weight_at_age", "spawning_weight_at_age", "maturity_at_age", "natural_mortality_at_age",
                     "landings_numbers_at_age", "landings_fraction_at_age",
                     "landings_weight_at_age", "discard_weight_at_age")
  } else {
    has_age_group <- if ("age_group" %in% names(x)) {
      !is.na(x$age_group) & nzchar(x$age_group)
    } else {
      rep(FALSE, nrow(x))
    }
    x$measure %in% c("numbers_at_age", "biomass_at_age", "fishing_mortality_at_age",
                     "natural_mortality_at_age") & !has_age_group
  }
  age <- suppressWarnings(as.numeric(x$age))
  if (any(age_required & (is.na(age) | !is.finite(age) | age != as.integer(age)))) {
    stop(name, "$age must contain whole ages for age-specific measures.", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(x$value))
  signed <- if (name == "inputs") x$type == "covariate" &
    x$measure == "environmental_covariate" & x$basis == "native_covariate" else rep(FALSE, nrow(x))
  if (anyNA(value) || any(!is.finite(value)) || any(value < 0 & !signed)) {
    stop(name, "$value must contain finite numbers; negative values are allowed only for environmental covariates.", call. = FALSE)
  }
  if (name == "inputs") {
    if (any(x$type == "maturity" & (value < 0 | value > 1))) {
      stop("Maturity values must be proportions from 0 to 1.", call. = FALSE)
    }
    sampling_time <- suppressWarnings(as.numeric(x$sampling_time))
    if (any(!is.na(x$sampling_time) & is.na(sampling_time)) ||
        any(!is.na(sampling_time) & (!is.finite(sampling_time) |
                                     sampling_time < 0 | sampling_time > 1))) {
      stop("inputs$sampling_time must be numeric in [0, 1] where present.", call. = FALSE)
    }
    proportion <- x$basis %in% c("proportion", "proportion_numbers", "proportion_biomass")
    if (any(proportion & value > 1)) {
      stop("Input proportions must lie in [0, 1].", call. = FALSE)
    }
    if (any(x$measure == "log_index_sd" & value <= 0)) {
      stop("Supplied log-index SDs must be positive.", call. = FALSE)
    }
  } else {
    if ("age_group" %in% names(x) && any(x$measure %in% c("Fbar", "Mbar") &
        (is.na(x$age_group) | !nzchar(x$age_group)))) {
      stop("Grouped Fbar and Mbar outputs must identify their age_group.", call. = FALSE)
    }
    for (field in c("se", "lwr", "upr")) {
      z <- suppressWarnings(as.numeric(x[[field]]))
      if (any(!is.na(x[[field]]) & (is.na(z) | !is.finite(z)))) {
        stop("outputs$", field, " must contain finite numbers or blanks.", call. = FALSE)
      }
      if (field == "se" && any(z < 0, na.rm = TRUE)) {
        stop("outputs$se cannot be negative.", call. = FALSE)
      }
    }
    both <- !is.na(x$lwr) & !is.na(x$upr)
    if (any(as.numeric(x$lwr[both]) > as.numeric(x$upr[both]))) {
      stop("outputs confidence limits have lwr greater than upr.", call. = FALSE)
    }
  }
}

composition <- inputs$measure %in% c("proportion_at_length", "conditional_proportion_at_age")
if (any(composition)) {
  missing <- setdiff(composition_dimensions, names(inputs))
  if (length(missing)) stop("Composition inputs lack dimensions: ", paste(missing, collapse=", "))
  z <- inputs[composition, , drop = FALSE]
  for (field in setdiff(composition_dimensions, "observation_id")) {
    values <- suppressWarnings(as.numeric(z[[field]]))
    if (any(!is.na(z[[field]]) & (is.na(values) | !is.finite(values)))) {
      stop("Composition inputs$", field, " must contain finite numbers or blanks.")
    }
    z[[field]] <- values
  }
  for (field in c("age_error", "partition")) {
    if (any(z[[field]] < 0 | z[[field]] != floor(z[[field]]), na.rm = TRUE)) {
      stop("Composition inputs$", field, " must contain non-negative integer codes.")
    }
  }
  if (anyNA(z$observation_id) || any(!nzchar(z$observation_id)) ||
      anyNA(z$sample_size) || any(z$sample_size <= 0)) {
    stop("Composition inputs need observation IDs and positive supplied sample sizes.")
  }
  length_rows <- z$measure == "proportion_at_length"
  if (any(length_rows & (is.na(z$length_bin) | z$length_bin < 0))) {
    stop("Length compositions need non-negative native length-bin labels.")
  }
  age_rows <- z$measure == "conditional_proportion_at_age"
  if (any(age_rows & (is.na(z$length_bin_lower) | is.na(z$length_bin_upper) |
                     z$length_bin_upper < z$length_bin_lower | is.na(z$age_error)))) {
    stop("Conditional age compositions need ordered conditioning bins and age-error codes.")
  }
}

cat("Database valid: ", nrow(stocks), " stocks, ", nrow(assessments),
    " assessments, ", nrow(assumptions), " assumptions, ", nrow(inputs),
    " inputs, and ", nrow(outputs), " outputs.\n", sep = "")
for (id in assessments$assessment_id) {
  x <- inputs[inputs$assessment_id == id, , drop = FALSE]
  surveys <- unique(x$survey[!is.na(x$survey) & nzchar(x$survey)])
  years <- suppressWarnings(as.numeric(x$year))
  ages <- suppressWarnings(as.numeric(x$age))
  cat(id, ": input rows by type: ",
      paste(names(table(x$type)), as.integer(table(x$type)), collapse = "; "),
      "; surveys: ", if (length(surveys)) paste(surveys, collapse = "; ") else "none",
      "; years: ", if (any(is.finite(years))) paste(range(years[is.finite(years)]), collapse = "-") else "unknown",
      "; ages: ", if (any(is.finite(ages))) paste(range(ages[is.finite(ages)]), collapse = "-") else "unknown",
      "\n", sep = "")
}

