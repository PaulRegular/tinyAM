# Import the detailed 2016 4T American plaice assessment and 2024 indicator update.

root <- file.path("analysis", "comp_assessments")
db <- file.path(root, "database")
cache <- file.path(root, "source_cache", "dfo_american_plaice_4t_2016")
assessment_id <- "dfo_american_plaice_4t_2016"
indicator_id <- "dfo_american_plaice_4t_2024_indicators"
stock_id <- "ck_dfo_platessoides_4t"
assessment_url <- "https://waves-vagues.dfo-mpo.gc.ca/Library/365839.pdf"
indicator_url <- "https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41282851.pdf"
abort <- function(...) stop(paste0(...), call. = FALSE)

assessment_file <- file.path(cache, "assessment_2016_057.pdf")
indicator_file <- file.path(cache, "indicator_update_2024_2025_010.pdf")
if (!file.exists(assessment_file) || !file.exists(indicator_file)) {
  abort("A required cached PDF is missing from ", cache)
}
if (!requireNamespace("pdftools", quietly = TRUE)) {
  abort("Install pdftools to extract the published assessment tables.")
}
assessment_text <- pdftools::pdf_text(assessment_file)
indicator_text <- pdftools::pdf_text(indicator_file)
normalize_text <- function(x) gsub("[[:space:]]+", " ", paste(x, collapse = " "))
assessment_body <- normalize_text(assessment_text)
indicator_body <- normalize_text(indicator_text)

parse_number <- function(token) {
  if (tolower(token) %in% c("na", "n/a", "-", "–", "—")) return(NA_real_)
  token <- gsub(",", "", token, fixed = TRUE)
  if (!grepl("^[+-]?[0-9]+(\\.[0-9]+)?$", token)) {
    abort("Unrecognized table value: ", token)
  }
  as.numeric(token)
}

read_table <- function(number, years, n_ages, extra_columns = 0L) {
  headings <- which(vapply(assessment_text, function(page) {
    grepl(paste0("Table ", number, "."), page, fixed = TRUE)
  }, logical(1)))
  if (!length(headings)) abort("Could not locate Table ", number, ".")
  start <- tail(headings, 1L)
  lines <- trimws(unlist(strsplit(assessment_text[[start]], intToUtf8(10L))))
  rows <- lines[grepl("^(19|20)[0-9]{2}[[:space:]]+", lines)]
  row_years <- as.integer(sub("^((19|20)[0-9]{2}).*", "\\1", rows))
  rows <- rows[row_years %in% years]
  row_years <- row_years[row_years %in% years]
  if (anyDuplicated(row_years) || !identical(row_years, as.integer(years))) {
    abort("Unexpected year coverage in Table ", number, ".")
  }

  expected <- 1L + n_ages + extra_columns
  parsed <- lapply(rows, function(line) {
    fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    if (length(fields) != expected) {
      abort("Unexpected column count in Table ", number, ": ", line)
    }
    vapply(fields[-1L], parse_number, numeric(1))
  })
  values <- do.call(rbind, parsed)
  rownames(values) <- as.character(years)
  colnames(values) <- c(as.character(4L + seq_len(n_ages) - 1L),
                        if (extra_columns) "total")
  list(year = as.integer(years), value = values)
}

years <- 1976:2015
ages <- 4:16
catch_table <- read_table(2L, years, length(ages), extra_columns = 1L)
rv_table <- read_table(3L, years, length(ages))
sentinel_years <- 2003:2015
sentinel_table <- read_table(4L, sentinel_years, length(ages))
weight_table <- read_table(5L, years, length(ages))
maturity_table <- read_table(6L, years, length(ages))
n_table <- read_table(7L, years, length(ages))

if (anyNA(catch_table$value) || anyNA(sentinel_table$value) ||
    anyNA(weight_table$value) || anyNA(maturity_table$value) ||
    anyNA(n_table$value[, seq_along(ages), drop = FALSE])) {
  abort("A required numeric source table contains an unexpected missing value.")
}
if (!all(is.na(rv_table$value["2003", ])) ||
    anyNA(rv_table$value[setdiff(as.character(years), "2003"), , drop = FALSE])) {
  abort("Unexpected missing cells in the RV survey table.")
}
if (any(abs(catch_table$value[, "total"] -
            rowSums(catch_table$value[, as.character(ages), drop = FALSE])) > 3)) {
  abort("Printed catch totals differ from age sums by more than 3 thousand fish.")
}
if (any(maturity_table$value < 0 | maturity_table$value > 1)) {
  abort("Maturity values must be proportions between zero and one.")
}
if (!identical(unname(weight_table$value["2003", "16"]), 0)) {
  abort("The reported 2003 age-16+ weight cell changed; review the source before importing.")
}
if (!grepl("2015 spawning stock biomass of 55,023 t", assessment_body, fixed = TRUE)) {
  abort("The reported 2015 SSB value was not found in the assessment text.")
}
if (!grepl("currently at 1.25 kg per tow", indicator_body, fixed = TRUE)) {
  abort("The reported 2024 biomass indicator value was not found in the update text.")
}

empty_inputs <- function(n) data.frame(
  assessment_id = rep(assessment_id, n), type = "", measure = "", basis = "",
  fleet = "", survey = "", sex = "", region = "4T", season = "",
  year = NA_real_, year_basis = "calendar_year", age = NA_real_, value = NA_real_,
  unit = "", sampling_time = NA_real_, source_type = "official_table",
  source_reference = "", transformation = "", notes = "",
  observation_id = "", length_bin = "", length_bin_lower = NA_real_,
  length_bin_upper = NA_real_, sample_size = NA_real_, age_error = "",
  partition = "", stringsAsFactors = FALSE
)

long_input <- function(x, ages, type, measure, basis, unit,
                       source_reference, notes, fleet = "", survey = "",
                       season = "", sampling_time = NA_real_,
                       transformation = "") {
  grid <- expand.grid(year_row = seq_along(x$year), age_column = seq_along(ages))
  grid <- grid[order(grid$year_row, grid$age_column), , drop = FALSE]
  values <- as.vector(t(x$value[, seq_along(ages), drop = FALSE]))
  keep <- is.finite(values)
  grid <- grid[keep, , drop = FALSE]
  out <- empty_inputs(sum(keep))
  out$type <- type
  out$measure <- measure
  out$basis <- basis
  out$fleet <- fleet
  out$survey <- survey
  out$season <- season
  out$year <- x$year[grid$year_row]
  out$age <- ages[grid$age_column]
  out$value <- values[keep]
  out$unit <- unit
  out$sampling_time <- sampling_time
  out$source_reference <- source_reference
  out$notes <- notes
  out$transformation <- transformation
  out
}

input_rows <- do.call(rbind, list(
  long_input(catch_table, ages, "catch", "numbers_at_age", "numbers",
    "thousand fish", "DFO Research Document 2016/057, Table 2, p. 13",
    "Commercial catch-at-age includes landings and estimated discarded catch. Age 16 is the 16+ group." ,
    fleet = "Commercial fishery"),
  long_input(rv_table, ages, "index", "numbers_at_age", "numbers",
    "thousand fish", "DFO Research Document 2016/057, Table 3, p. 14",
    "Standardized DFO research-vessel abundance-at-age index. The 2003 cells are NA in the source and were excluded because the vessel was uncalibrated; only finite source observations are stored. Age 16 is 16+.",
    survey = "DFO Research Vessel Trawl Survey", season = "fall",
    sampling_time = 0.75,
    transformation = "September survey timing represented as 0.75 of the model year; the exact annual sampling date is not reported."),
  long_input(sentinel_table, ages, "index", "numbers_at_age", "numbers",
    "thousand fish", "DFO Research Document 2016/057, Table 4, p. 15",
    "Standardized mobile-gear sentinel survey abundance-at-age index, 2003-2015. Age 16 is 16+.",
    survey = "Mobile Gear Sentinel Survey", season = "fall",
    sampling_time = 0.75,
    transformation = "September survey timing represented as 0.75 of the model year; the exact annual sampling date is not reported."),
  long_input(weight_table, ages, "weight", "weight_at_age", "kg_per_fish",
    "kg/fish", "DFO Research Document 2016/057, Table 5, p. 16",
    "Population weight-at-age adjusted from September RV samples to January 1. The reported 2003 age-16+ value of 0.000 kg is preserved verbatim; the report does not explain whether it represents no observations or a model input.",
    transformation = "Source table's January 1 weight-at-age values retained; no interpolation applied."),
  long_input(maturity_table, ages, "maturity", "maturity_at_age", "proportion",
    "proportion", "DFO Research Document 2016/057, Table 6, p. 17",
    "Annual maturity-at-age proportions. The source applies median post-1996 logistic maturity coefficients to 1976-1996 and copies 2012 values to 2013-2015 because recent maturity coding was unreliable. Age 16 is 16+.")
))

empty_assumptions <- function(n) data.frame(
  assessment_id = rep(assessment_id, n), component = "", fleet = "",
  survey = "", sex = "", region = "4T", season = "", setting = "",
  value = "", source_reference = "", notes = "", stringsAsFactors = FALSE
)
assumption <- function(component, setting, value, source_reference, notes = "",
                       survey = "") {
  out <- empty_assumptions(1L)
  out$component <- component
  out$setting <- setting
  out$value <- as.character(value)
  out$source_reference <- source_reference
  out$notes <- notes
  out$survey <- survey
  out
}

assumptions <- do.call(rbind, list(
  assumption("model", "model_family", "Virtual population analysis (VPA)",
    "DFO Research Document 2016/057, pp. 5-7"),
  assumption("model", "software", "AD Model Builder (ADMB)",
    "DFO Research Document 2016/057, p. 6"),
  assumption("population", "model_years", "1976-2015",
    "DFO Research Document 2016/057, p. 6"),
  assumption("population", "model_ages", "4-16+ (13 age classes)",
    "DFO Research Document 2016/057, p. 6",
    "Age 16 is the plus group."),
  assumption("population", "recruitment_age", "4",
    "DFO Research Document 2016/057, pp. 7-8",
    "The separate stock-recruit analysis matches age-4 recruitment with SSB four years earlier."),
  assumption("catch", "catch_input", "Commercial catch-at-age including estimated discards",
    "DFO Research Document 2016/057, Table 2, p. 13",
    "Catch-at-age is an input to the VPA, not a catch-observation likelihood."),
  assumption("survey", "RV_index", "Age-specific abundance index; lognormal likelihood; age-specific q and observation variance",
    "DFO Research Document 2016/057, pp. 5-6",
    "There are 13 estimated catchability parameters and 13 observation-error variance parameters. 2003 is excluded because the survey vessel was uncalibrated.",
    survey = "DFO Research Vessel Trawl Survey"),
  assumption("survey", "RV_timing", "September; exact day not reported",
    "DFO Research Document 2016/057, pp. 4-5",
    "The database stores a 0.75 within-year approximation for September timing.",
    survey = "DFO Research Vessel Trawl Survey"),
  assumption("survey", "sentinel_index", "Age-specific abundance index; lognormal likelihood; age-specific q and observation variance",
    "DFO Research Document 2016/057, pp. 5-6",
    "The mobile-gear sentinel index covers 2003-2015. There are 13 estimated q and 13 observation-error variance parameters.",
    survey = "Mobile Gear Sentinel Survey"),
  assumption("survey", "sentinel_timing", "September; exact day not reported",
    "DFO Research Document 2016/057, p. 5",
    "The survey adopted the RV survey's stratified-random design in 2003; database timing 0.75 is a monthly approximation.",
    survey = "Mobile Gear Sentinel Survey"),
  assumption("mortality", "M_age_groups", "4-9 and 10+",
    "DFO Research Document 2016/057, p. 6",
    "Two independent annual M time series are estimated."),
  assumption("mortality", "M_process", "Log-scale random walk; M[y] = M[y-1] * exp(Mdev[y])",
    "DFO Research Document 2016/057, p. 6",
    "Mdev is normally distributed with mean zero and fixed SD 0.075."),
  assumption("mortality", "initial_M_prior_mean", "0.6 for ages 4-9; 0.3 for ages 10+",
    "DFO Research Document 2016/057, p. 6",
    "Prior SDs for initial M are not reported in the document."),
  assumption("mortality", "initial_M_prior_sd", "unknown",
    "DFO Research Document 2016/057, p. 6",
    "The report describes normal priors but does not state their standard deviations."),
  assumption("biology", "weight_timing", "Population weights refer to January 1",
    "DFO Research Document 2016/057, p. 4 and Table 5, p. 16",
    "Survey samples are collected in September; reported weights are adjusted to the population model's January 1 state."),
  assumption("biology", "maturity_history", "Median post-1996 logistic coefficients applied to 1976-1996; 2012 values copied to 2013-2015",
    "DFO Research Document 2016/057, p. 4 and Table 6, p. 17",
    "Earlier maturity staging and recent coding were considered unreliable."),
  assumption("model", "parameter_count", "144 estimated parameters",
    "DFO Research Document 2016/057, p. 6",
    "The source lists 12 terminal abundance parameters, 80 M parameters, 26 age-specific q parameters, and 26 age-specific survey error-variance parameters."),
  assumption("outputs", "point_estimate", "Median of 5,000 retained MCMC samples",
    "DFO Research Document 2016/057, p. 7",
    "The run used 201,000 iterations, discarded 1,000 burn-in iterations, and retained every 40th subsequent sample."),
  assumption("outputs", "annual_output_availability", "N-at-age is tabulated; annual F, M, and SSB are shown only in figures",
    "DFO Research Document 2016/057, Tables 7 and Figures 21-24",
    "No native fitted object or tabulated annual F/M/SSB output series was located in the sources reviewed."),
  assumption("model_availability", "native_run", "No public native fitted object or input files located",
    "DFO Research Document 2016/057; source review 2026-10-20",
    "The report identifies an ADMB implementation but does not link to a production run archive."),
  assumption("biology", "weight_age16_2003", "reported 0.000 kg/fish",
    "DFO Research Document 2016/057, Table 5, p. 16",
    "The value is preserved as printed; its interpretation is unresolved. The RV age-index cells for 2003 are unavailable and excluded." )
))

empty_outputs <- function(n, assessment = assessment_id) data.frame(
  assessment_id = rep(assessment, n), type = "", measure = "",
  fleet = "", survey = "", sex = "", region = "4T", season = "",
  year = NA_real_, age = NA_real_, age_group = "", value = NA_real_,
  se = NA_real_, lwr = NA_real_, upr = NA_real_, unit = "",
  source_type = "official_table", source_reference = "", notes = "",
  stringsAsFactors = FALSE
)
add_age_outputs <- function(table, years, ages, type, measure, unit,
                             source_reference, notes) {
  grid <- expand.grid(year_row = seq_along(years), age_column = seq_along(ages))
  grid <- grid[order(grid$year_row, grid$age_column), , drop = FALSE]
  out <- empty_outputs(nrow(grid))
  out$type <- type
  out$measure <- measure
  out$year <- years[grid$year_row]
  out$age <- ages[grid$age_column]
  out$value <- as.vector(t(table$value[, seq_along(ages), drop = FALSE]))
  out$unit <- unit
  out$source_reference <- source_reference
  out$notes <- notes
  out
}
add_annual_output <- function(assessment, type, measure, years, values, unit,
                              source_type, source_reference, notes) {
  out <- empty_outputs(length(years), assessment)
  out$type <- type
  out$measure <- measure
  out$year <- years
  out$value <- values
  out$unit <- unit
  out$source_type <- source_type
  out$source_reference <- source_reference
  out$notes <- notes
  out
}

n_at_age <- n_table$value[, seq_along(ages), drop = FALSE]
outputs <- list(
  add_age_outputs(n_table, years, ages, "population", "numbers_at_age",
    "thousand fish", "DFO Research Document 2016/057, Table 7, p. 18",
    "January 1 population abundance; values are MCMC marginal medians. Age 16 is the 16+ group."),
  add_annual_output(assessment_id, "population", "total_numbers", years,
    rowSums(n_at_age), "thousand fish", "derived_source_output",
    "DFO Research Document 2016/057, Table 7, p. 18",
    "Sum of the reported marginal median N-at-age values for ages 4-16+; this is not the posterior median of total abundance and has no uncertainty."),
  add_annual_output(assessment_id, "recruitment", "recruitment", years,
    n_at_age[, "4"], "thousand fish", "derived_source_output",
    "DFO Research Document 2016/057, Table 7, p. 18",
    "Recruitment represented by January 1 abundance at age 4 in the same calendar year; no uncertainty is available."),
  add_annual_output(assessment_id, "biomass", "SSB", 2015L, 55023,
    "t", "official_document", "DFO Research Document 2016/057, p. 7",
    "2015 beginning-of-year spawning stock biomass, reported as the MCMC median. The source shows uncertainty only as a figure, not a numeric interval.")
)
outputs <- do.call(rbind, outputs)

indicator_output <- add_annual_output(
  indicator_id, "index", "relative_biomass_index", 2024L, 1.25,
  "kg/tow", "official_document", "DFO Science Response 2025/010, p. 3",
  "Reported three-year moving-average survey biomass indicator for commercial-size fish (at least 30 cm), shown for 2024. This survey indicator is not SSB and is kept separate from the 2016 VPA outputs.")

assessments <- data.frame(
  assessment_id = c(assessment_id, indicator_id),
  stock_id = stock_id,
  assessment_year = c(2016L, 2024L),
  terminal_year = c(2015L, 2024L),
  estimate_terminal_year = c(2015L, 2024L),
  assessment_type = c("full", "indicator_update"),
  model_family = c("VPA", ""),
  model_version = c("ADMB", ""),
  is_current = c(TRUE, FALSE),
  is_applied = c(TRUE, TRUE),
  framework_year = c(NA_integer_, NA_integer_),
  assessment_url = c(assessment_url, indicator_url),
  framework_url = c("", ""), data_url = c("", ""),
  model_url = c("", ""), repository_url = c("", ""),
  assumptions_status = c("partial", "partial"),
  inputs_status = c("complete", "partial"),
  outputs_status = c("partial", "partial"),
  notes = c(
    "Most recent detailed production assessment recovered. Official model inputs are represented, but annual F, M and SSB values and the ADMB production files are unavailable as machine-readable outputs.",
    "The 2024 interim indicator update was published in 2025. The year fields identify the indicator endpoint, not a fitted population-model terminal year. It reports a survey-derived three-year biomass-index average and is recorded separately from the 2016 VPA."),
  stringsAsFactors = FALSE
)

indicator_assumptions <- data.frame(
  assessment_id = indicator_id, component = "index", fleet = "", survey = "DFO Research Vessel Trawl Survey",
  sex = "", region = "4T", season = "fall", setting = "indicator_definition",
  value = "Three-year moving average of survey-derived biomass for fish >=30 cm; 2024 value 1.25 kg/tow",
  source_reference = "DFO Science Response 2025/010, pp. 3-4",
  notes = "Indicator and scaled reference-point update only; not an SSB estimate or a fitted population model.",
  stringsAsFactors = FALSE
)

replace_rows <- function(name, ids, rows) {
  path <- file.path(db, name)
  existing_names <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  missing <- setdiff(existing_names, names(rows))
  if (length(missing)) abort(name, " rows missing columns: ", paste(missing, collapse = ", "))
  rows <- rows[, existing_names, drop = FALSE]
  rows_file <- tempfile(fileext = ".csv")
  on.exit(unlink(rows_file), add = TRUE)
  write.csv(rows, rows_file, row.names = FALSE, na = "")
  appended <- readLines(rows_file, warn = FALSE, encoding = "bytes")[-1L]
  write_bytes_lines <- function(lines, connection) {
    if (length(lines)) {
      writeBin(charToRaw(paste0(paste(lines, collapse = "\r\n"), "\r\n")),
               connection)
    }
  }

  input <- file(path, open = "rb")
  on.exit(try(close(input), silent = TRUE), add = TRUE)
  temp_path <- tempfile(tmpdir = dirname(path))
  output <- file(temp_path, open = "wb")
  on.exit({
    try(close(output), silent = TRUE)
    unlink(temp_path)
  }, add = TRUE)
  header <- readLines(input, n = 1L, warn = FALSE, encoding = "bytes")
  write_bytes_lines(header, output)
  repeat {
    chunk <- readLines(input, n = 10000L, warn = FALSE, encoding = "bytes")
    if (!length(chunk)) break
    first_field <- sub('^"?([^",]*).*$','\\1', chunk)
    keep <- !first_field %in% ids
    if (any(keep)) write_bytes_lines(chunk[keep], output)
  }
  write_bytes_lines(appended, output)
  close(output)
  close(input)
  if (!file.copy(temp_path, path, overwrite = TRUE)) abort("Could not replace ", path)
  unlink(temp_path)
}

stocks <- read.csv(file.path(db, "stocks.csv"), stringsAsFactors = FALSE,
                   check.names = FALSE)
stock_row <- stocks[stocks$stock_id == stock_id, , drop = FALSE]
if (nrow(stock_row) != 1L) abort("Expected exactly one existing 4T American plaice stock record.")
stock_row$authority_stock_id <- "4T"
stock_row$scientific_name <- "Hippoglossoides platessoides"
stock_row$common_name <- "American plaice"
stock_row$notes <- paste(
  "Southern Gulf of St. Lawrence NAFO Division 4T American plaice.",
  "The most recent detailed production assessment recovered is the 2016 VPA; a later 2024 interim indicator update is recorded separately.",
  "Charbonneau-Keith Master.csv identifier DFO_platessoides_4T retained."
)

replace_rows("stocks.csv", stock_id, stock_row)
replace_rows("assessments.csv", c(assessment_id, indicator_id), assessments)
replace_rows("assumptions.csv", c(assessment_id, indicator_id),
             rbind(assumptions, indicator_assumptions))
replace_rows("inputs.csv", assessment_id, input_rows)
replace_rows("outputs.csv", c(assessment_id, indicator_id),
             rbind(outputs, indicator_output))

cat("Imported ", nrow(input_rows), " detailed 2016 inputs and ",
    nrow(outputs), " 2016 outputs, plus one 2024 indicator value.\n", sep = "")
