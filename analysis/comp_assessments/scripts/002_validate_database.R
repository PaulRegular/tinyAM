root <- file.path("analysis", "comp_assessments", "database")
read_table <- function(name) read.csv(file.path(root, name), stringsAsFactors = FALSE,
                                      na.strings = c("", "NA"), check.names = FALSE)
stocks <- read_table("stocks.csv")
assessments <- read_table("assessments.csv")
assumptions <- read_table("assumptions.csv")
inputs <- read_table("inputs.csv")
outputs <- read_table("outputs.csv")

schemas <- list(
  stocks = c("stock_id", "charbonneau_id", "authority", "authority_stock_id",
             "scientific_name", "common_name", "area", "region", "ocean", "notes"),
  assessments = c("assessment_id", "stock_id", "assessment_year", "terminal_year",
                  "assessment_type", "model_family", "model_version", "is_current",
                  "is_production", "framework_year", "assessment_url", "framework_url",
                  "data_url", "model_url", "repository_url", "assumptions_status",
                  "inputs_status", "outputs_status", "notes"),
  assumptions = c("assessment_id", "component", "setting", "value", "source_reference", "notes"),
  inputs = c("assessment_id", "type", "fleet", "survey", "sex", "region", "season",
             "year", "age", "value", "unit", "samp_time", "source_type",
             "source_reference", "notes"),
  outputs = c("assessment_id", "type", "fleet", "survey", "sex", "region", "season",
              "year", "age", "age_group", "value", "se", "lwr", "upr", "unit", "source_type",
              "source_reference", "notes")
)
tables <- list(stocks = stocks, assessments = assessments, assumptions = assumptions,
               inputs = inputs, outputs = outputs)
for (name in names(tables)) {
  missing <- setdiff(schemas[[name]], names(tables[[name]]))
  if (length(missing)) stop(name, ".csv is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

unique_ids <- function(x, column, table) {
  id <- x[[column]]
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id)) {
    stop(table, "$", column, " must contain unique, non-empty values.", call. = FALSE)
  }
}
unique_ids(stocks, "stock_id", "stocks")
unique_ids(assessments, "assessment_id", "assessments")
if (any(!assessments$stock_id %in% stocks$stock_id)) stop("assessments.csv has an unknown stock_id.", call. = FALSE)
assessment_ids <- assessments$assessment_id
for (name in c("assumptions", "inputs", "outputs")) {
  x <- tables[[name]]
  if (anyNA(x$assessment_id) || any(!x$assessment_id %in% assessment_ids)) {
    stop(name, ".csv has an unknown or empty assessment_id.", call. = FALSE)
  }
}

duplicates <- list(
  assumptions = c("assessment_id", "component", "setting"),
  inputs = c("assessment_id", "type", "fleet", "survey", "sex", "region", "season", "year", "age"),
  outputs = c("assessment_id", "type", "fleet", "survey", "sex", "region", "season", "year", "age", "age_group")
)
for (name in names(duplicates)) {
  if (anyDuplicated(tables[[name]][duplicates[[name]]])) {
    stop(name, ".csv has duplicate canonical rows.", call. = FALSE)
  }
}

status_values <- c("complete", "partial", "unknown", "not_applicable")
for (field in c("assumptions_status", "inputs_status", "outputs_status")) {
  value <- assessments[[field]]
  if (any(!is.na(value) & !value %in% status_values)) {
    stop("assessments$", field, " must be one of: ", paste(status_values, collapse = ", "), call. = FALSE)
  }
}
logical_value <- function(x) tolower(as.character(x))
for (field in c("is_current", "is_production")) {
  value <- logical_value(assessments[[field]])
  if (any(!is.na(value) & !value %in% c("true", "false", "1", "0"))) {
    stop("assessments$", field, " must be TRUE/FALSE or 1/0.", call. = FALSE)
  }
}
current <- logical_value(assessments$is_current) %in% c("true", "1")
production <- logical_value(assessments$is_production) %in% c("true", "1")
if (any(current & !production, na.rm = TRUE)) stop("A current assessment must be marked production.", call. = FALSE)
for (i in seq_len(nrow(assessments))) {
  assessment_id <- assessments$assessment_id[i]
  tables_for_status <- list(assumptions = assumptions, inputs = inputs, outputs = outputs)
  status_fields <- c(assumptions = "assumptions_status", inputs = "inputs_status",
                     outputs = "outputs_status")
  for (table_name in names(status_fields)) {
    status <- assessments[[status_fields[[table_name]]]][i]
    has_rows <- any(tables_for_status[[table_name]]$assessment_id == assessment_id)
    if (!is.na(status) && status == "complete" && !has_rows) {
      stop(assessment_id, " is marked complete for ", table_name, " but has no rows.", call. = FALSE)
    }
    if (!is.na(status) && status == "not_applicable" && has_rows) {
      stop(assessment_id, " is marked not_applicable for ", table_name, " but has rows.", call. = FALSE)
    }
  }
}
year_fields <- c("assessment_year", "terminal_year", "framework_year")
for (field in year_fields) {
  value <- suppressWarnings(as.numeric(assessments[[field]]))
  raw <- assessments[[field]]
  if (any(!is.na(raw) & (is.na(value) | value != as.integer(value)))) {
    stop("assessments$", field, " must contain whole years or blanks.", call. = FALSE)
  }
}
assessment_year <- as.numeric(assessments$assessment_year)
terminal_year <- as.numeric(assessments$terminal_year)
if (any(terminal_year > assessment_year, na.rm = TRUE)) {
  stop("terminal_year cannot be later than assessment_year.", call. = FALSE)
}
current_stocks <- assessments$stock_id[current]
if (anyDuplicated(current_stocks)) {
  repeated <- unique(current_stocks[duplicated(current_stocks)])
  for (stock in repeated) {
    rows <- assessments[assessments$stock_id == stock & current, , drop = FALSE]
    if (anyNA(rows$notes) || any(!nzchar(rows$notes))) {
      stop("Multiple current assessments for ", stock, " need an explanation in notes.", call. = FALSE)
    }
  }
}

allowed_sources <- c("native_model", "official_machine_readable", "official_table",
                     "digitized", "reconstructed", "charbonneau_seed")
allowed_input_types <- c("catch", "catch_at_age", "landings", "index", "weight", "catch_weight", "maturity", "maturity_cohort", "M")
allowed_output_types <- c("N", "F", "M", "SSB", "biomass", "recruitment", "Fbar", "Mbar", "q",
                          "predicted_catch", "predicted_index")
for (name in c("inputs", "outputs")) {
  x <- tables[[name]]
  if (anyNA(x$source_type) || any(!x$source_type %in% allowed_sources)) {
    stop(name, "$source_type must use the controlled source labels.", call. = FALSE)
  }
  allowed_types <- if (name == "inputs") allowed_input_types else allowed_output_types
  if (anyNA(x$type) || any(!x$type %in% allowed_types)) {
    stop(name, "$type must use the documented assessment value labels.", call. = FALSE)
  }
  if (name == "outputs" && any(x$type %in% c("Fbar", "Mbar") &
      (is.na(x$age_group) | !nzchar(x$age_group)))) {
    stop("Grouped Fbar and Mbar outputs must identify their age_group.", call. = FALSE)
  }
  if (anyNA(x$source_reference) || any(!nzchar(x$source_reference))) {
    stop(name, "$source_reference must identify the source location.", call. = FALSE)
  }
  for (field in c("year", "age", "value")) {
    raw <- x[[field]]
    numeric <- suppressWarnings(as.numeric(raw))
    missing_index <- is.na(numeric)
    required <- if (field == "year") {
      rep(TRUE, length(numeric))
    } else if (field == "age" && name == "inputs") {
      x$type != "landings"
    } else if (field == "age" && name == "outputs") {
      x$type %in% c("N", "F", "M", "q", "predicted_catch", "predicted_index")
    } else {
      rep(FALSE, length(numeric))
    }
    if (any(!is.na(raw) & is.na(numeric)) ||
        (field %in% c("year", "age") &&
         (any(missing_index & required) || any(numeric[!missing_index] != as.integer(numeric[!missing_index]))))) {
      stop(name, "$", field, " must contain numeric", if (field == "value") " or blank" else " whole numbers", " values.", call. = FALSE)
    }
    if (field == "value" && any(numeric < 0, na.rm = TRUE)) {
      stop(name, "$value cannot be negative.", call. = FALSE)
    }
    if (any(!is.na(numeric) & !is.finite(numeric))) {
      stop(name, "$", field, " must contain finite values.", call. = FALSE)
    }
  }
  for (field in intersect(c("se", "lwr", "upr", "samp_time"), names(x))) {
    raw <- x[[field]]
    numeric <- suppressWarnings(as.numeric(raw))
    if (any(!is.na(raw) & is.na(numeric)) || any(!is.na(numeric) & !is.finite(numeric))) {
      stop(name, "$", field, " must contain numeric finite values or blanks.", call. = FALSE)
    }
    if (field == "samp_time" && any(numeric < 0 | numeric > 1, na.rm = TRUE)) {
      stop("inputs$samp_time must be between 0 and 1.", call. = FALSE)
    }
    if (field == "se" && any(numeric < 0, na.rm = TRUE)) {
      stop("outputs$se cannot be negative.", call. = FALSE)
    }
  }
  if (name == "inputs" && any(x$type %in% c("maturity", "maturity_cohort") &
      !is.na(x$value) & (as.numeric(x$value) < 0 | as.numeric(x$value) > 1))) {
    stop("maturity input values must be proportions from 0 to 1.", call. = FALSE)
  }
  both_bounds <- !is.na(x$lwr) & !is.na(x$upr)
  if (any(as.numeric(x$lwr[both_bounds]) > as.numeric(x$upr[both_bounds]))) {
    stop(name, " confidence limits have lwr greater than upr.", call. = FALSE)
  }
}

cat("Database valid: ", nrow(stocks), " stocks, ", nrow(assessments),
    " assessments, ", nrow(assumptions), " assumptions, ", nrow(inputs),
    " inputs, and ", nrow(outputs), " outputs.\n", sep = "")
