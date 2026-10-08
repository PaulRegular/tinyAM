assessment_id <- "dfo_cod_3pn4rs_2025"
database_dir <- file.path("analysis", "comp_assessments", "database")
source_dir <- file.path("analysis", "comp_assessments", "source_cache",
                        assessment_id)
source_file <- file.path(source_dir, "assessment.txt")
assessment_url <- paste0(
  "https://publications.gc.ca/collections/collection_2026/mpo-dfo/",
  "fs70-5/Fs70-5-2026-010-eng.pdf"
)
framework_url <- paste0(
  "https://publications.gc.ca/collections/collection_2025/mpo-dfo/",
  "fs70-5/Fs70-5-2025-074-eng.pdf"
)

if (!file.exists(source_file)) stop("Cached assessment text is missing.")
lines <- readLines(source_file, warn = FALSE)

read_table <- function(number, next_number, value_count) {
  start <- grep(paste0("^Table ", number, "\\."), trimws(lines))
  end <- grep(paste0("^Table ", next_number, "\\."), trimws(lines))
  start <- start[[1]]
  end <- end[end > start][[1]]
  rows <- lines[seq.int(start + 1L, end - 1L)]
  rows <- rows[grepl("^[[:space:]]*(19|20)[0-9]{2}[[:space:]]+", rows)]
  values <- strsplit(gsub(",", "", trimws(rows)), "[[:space:]]+")
  if (!length(values) || any(lengths(values) != value_count + 1L)) {
    stop("Unexpected row structure in Table ", number, ".")
  }
  values <- unlist(values)
  values[values == "-"] <- NA_character_
  matrix(as.numeric(values), ncol = value_count + 1L, byrow = TRUE)
}

catch_weight <- read_table(24, 25, 10)
sentinel <- read_table(28, 29, 11)
catch <- read_table(30, 31, 10)
biomass <- read_table(32, 33, 12)
abundance <- read_table(33, 34, 12)
fishing_mortality <- read_table(35, 36, 12)
natural_mortality <- read_table(36, 37, 10)

stopifnot(
  all(catch_weight[, 1] == 1974:2024),
  all(sentinel[, 1] == 1995:2024),
  all(catch[, 1] == 1974:2024),
  all(biomass[, 1] == 1973:2024),
  all(abundance[, 1] == 1973:2024),
  all(fishing_mortality[, 1] == 1973:2024),
  all(natural_mortality[, 1] == 1973:2024)
)

age_number <- 2:11
age_group <- ifelse(age_number == 11, "11+", "")
assessment_ref <- "DFO Research Document 2026/010"

long_at_age <- function(table, columns, table_number, measure_type, measure,
                        unit, note) {
  rows <- expand.grid(
    year = table[, 1],
    index = seq_along(age_number),
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- rows[order(rows$year, rows$index), ]
  data.frame(
    assessment_id = assessment_id,
    type = measure_type,
    measure = measure,
    fleet = "",
    survey = "",
    sex = "",
    region = "",
    season = "",
    year = rows$year,
    age = age_number[rows$index],
    age_group = age_group[rows$index],
    value = as.vector(t(table[, columns, drop = FALSE])),
    se = NA_real_,
    lwr = NA_real_,
    upr = NA_real_,
    unit = unit,
    source_type = "official_table",
    source_reference = paste(assessment_ref, "Table", table_number),
    notes = note,
    stringsAsFactors = FALSE
  )
}

catch_rows <- expand.grid(
  year = catch[, 1],
  index = seq_along(age_number),
  KEEP.OUT.ATTRS = FALSE
)
catch_rows <- catch_rows[order(catch_rows$year, catch_rows$index), ]
inputs <- data.frame(
  assessment_id = assessment_id,
  type = "catch",
  measure = "numbers_at_age",
  basis = "numbers",
  fleet = "",
  survey = "",
  sex = "",
  region = "",
  season = "",
  year = catch_rows$year,
  year_basis = "calendar_year",
  age = age_number[catch_rows$index],
  value = as.vector(t(catch[, 2:11, drop = FALSE])),
  unit = "thousand fish",
  sampling_time = NA_real_,
  source_type = "official_table",
  source_reference = paste(assessment_ref, "Table 30"),
  transformation = "",
  notes = ifelse(
    age_number[catch_rows$index] == 11,
    "Published catch-at-age input; age 11 is the 11+ group. The table does not separate source fleets.",
    "Published catch-at-age input. The table does not separate source fleets."
  ),
  observation_id = "",
  length_bin = NA_real_,
  length_bin_lower = NA_real_,
  length_bin_upper = NA_real_,
  sample_size = NA_real_,
  age_error = NA_real_,
  partition = NA_real_,
  stringsAsFactors = FALSE
)

table_inputs <- function(table, ages, type, measure, basis, unit, number, note) {
  rows <- inputs[rep(1L, nrow(table) * length(ages)), , drop = FALSE]
  rows$type <- type
  rows$measure <- measure
  rows$basis <- basis
  rows$year <- rep(table[, 1], each = length(ages))
  rows$age <- rep(ages, times = nrow(table))
  rows$value <- as.vector(t(table[, -1, drop = FALSE]))
  rows$unit <- unit
  rows$source_reference <- paste(assessment_ref, "Table", number)
  rows$notes <- paste(note, "Age 11 represents 11+.")
  rows[!is.na(rows$value), , drop = FALSE]
}

weights <- table_inputs(
  catch_weight, 2:11, "catch_weight", "weight_at_age", "kg_per_fish", "kg", 24,
  paste("Commercial catch weights, not beginning-of-year stock weights.",
        "Printed zeros are retained; seven '-' cells are missing and omitted.")
)
index <- table_inputs(
  sentinel, 1:11, "index", "numbers_at_age", "numbers", "mean numbers per tow", 28,
  paste("Published July Sentinel mobile index; tow geometry is standardized",
        "to 54 ft horizontal opening and 1.25 nautical miles.",
        "Age 1 is reported but the accepted model fits ages 2-11+.",
        "Sampling time 0.54 is a mid-July approximation, not the recovered model setting.")
)
index$survey <- "Sentinel mobile"
index$region <- "3Pn4RS"
index$season <- "July"
index$sampling_time <- 0.54
inputs <- rbind(inputs, weights, index)

fixed_M <- do.call(rbind, lapply(1973:2024, function(year) {
  fixed_ages <- c(2, 3)
  fixed_values <- c(1, 0.65)
  if (year <= 1983) {
    fixed_ages <- c(fixed_ages, 4:11)
    fixed_values <- c(fixed_values, 0.45, rep(0.15, 7))
  }
  data.frame(year = year, age = fixed_ages, value = fixed_values)
}))
inputs <- rbind(inputs, data.frame(
  assessment_id = assessment_id,
  type = "M",
  measure = "natural_mortality_at_age",
  basis = "per_year",
  fleet = "",
  survey = "",
  sex = "",
  region = "",
  season = "",
  year = fixed_M$year,
  year_basis = "calendar_year",
  age = fixed_M$age,
  value = fixed_M$value,
  unit = "per_year",
  sampling_time = NA_real_,
  source_type = "official_table",
  source_reference = paste(assessment_ref, "Table 36"),
  transformation = "",
  notes = ifelse(
    fixed_M$age %in% 2:3,
    "Fixed model assumption for all years; age 11 represents 11+.",
    "Fixed initial M assumption through 1983; age 11 represents 11+."
  ),
  observation_id = "",
  length_bin = NA_real_,
  length_bin_lower = NA_real_,
  length_bin_upper = NA_real_,
  sample_size = NA_real_,
  age_error = NA_real_,
  partition = NA_real_,
  stringsAsFactors = FALSE
))

make_output <- function(type, measure, year, value, unit, source_table,
                        age = NA_real_, age_group = "", notes = "") {
  data.frame(
    assessment_id = assessment_id,
    type = type,
    measure = measure,
    fleet = "",
    survey = "",
    sex = "",
    region = "",
    season = "",
    year = year,
    age = age,
    age_group = age_group,
    value = value,
    se = NA_real_,
    lwr = NA_real_,
    upr = NA_real_,
    unit = unit,
    source_type = "official_table",
    source_reference = paste(assessment_ref, "Table", source_table),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

outputs <- rbind(
  long_at_age(
    abundance, 2:11, 33, "population", "numbers_at_age", "thousand fish",
    "Beginning-of-year abundance; age 11 is the 11+ group."
  ),
  long_at_age(
    biomass, 2:11, 32, "biomass", "biomass_at_age", "tonnes",
    "Beginning-of-year biomass; age 11 is the 11+ group."
  ),
  long_at_age(
    fishing_mortality, 2:11, 35, "mortality", "fishing_mortality_at_age",
    "per_year", "Age 11 is the 11+ group."
  ),
  make_output(
    "mortality", "natural_mortality_at_age",
    rep(natural_mortality[natural_mortality[, 1] >= 1984, 1], each = 8),
    as.vector(t(natural_mortality[natural_mortality[, 1] >= 1984,
                                  4:11, drop = FALSE])),
    "per_year", 36, rep(4:11, times = 2024 - 1984 + 1),
    ifelse(rep(4:11, times = 2024 - 1984 + 1) == 11, "11+", ""),
    "Estimated M output; fixed ages 2-3 and pre-1984 values are stored as inputs."
  ),
  make_output("biomass", "SSB", biomass[, 1], biomass[, 13], "tonnes", 32,
              notes = "Spawning stock biomass."),
  make_output("biomass", "total_biomass", biomass[, 1], biomass[, 12],
              "tonnes", 32, age_group = "2+",
              notes = "Beginning-of-year biomass for ages 2+."),
  make_output("population", "total_numbers", abundance[, 1], abundance[, 12],
              "thousand fish", 33, age_group = "2+",
              notes = "Beginning-of-year abundance for ages 2+."),
  make_output("population", "total_numbers", abundance[, 1], abundance[, 13],
              "thousand fish", 33, age_group = "5+",
              notes = "Beginning-of-year abundance for ages 5+."),
  make_output("recruitment", "recruitment", abundance[, 1], abundance[, 2],
              "thousand fish", 33, age = 2,
              notes = "Recruitment is reported at model age 2; year is the year fish reach age 2."),
  make_output("mortality", "Fbar", fishing_mortality[, 1],
              fishing_mortality[, 12], "per_year", 35, age_group = "4-6",
              notes = "Source-reported mean F for ages 4-6."),
  make_output("mortality", "Fbar", fishing_mortality[, 1],
              fishing_mortality[, 13], "per_year", 35, age_group = "6-9",
              notes = "Source-reported mean F for ages 6-9.")
)

stock <- data.frame(
  stock_id = "dfo_cod_3pn4rs",
  charbonneau_id = "DFO-GSL_morhua_3Pn, 4RS",
  authority = "DFO",
  authority_stock_id = "NAFO 3Pn4RS",
  scientific_name = "Gadus morhua",
  common_name = "Northern Gulf of St. Lawrence cod",
  area = "NAFO Subdivision 3Pn and Divisions 4RS",
  region = "Quebec Region",
  ocean = "Northwest Atlantic",
  notes = "Charbonneau-Keith identifier is DFO-GSL_morhua_3Pn, 4RS."
)

assessment <- data.frame(
  assessment_id = assessment_id,
  stock_id = stock$stock_id,
  assessment_year = 2025,
  terminal_year = 2024,
  estimate_terminal_year = 2024,
  assessment_type = "full_assessment",
  model_family = "state-space statistical catch-at-age",
  model_version = "2025 accepted assessment; RTMB model with 16 mean-F blocks",
  is_current = TRUE,
  is_applied = TRUE,
  framework_year = 2022,
  assessment_url = assessment_url,
  framework_url = framework_url,
  data_url = "",
  model_url = "",
  repository_url = "",
  assumptions_status = "partial",
  inputs_status = "partial",
  outputs_status = "partial",
  notes = paste(
    "Peer-reviewed February 18-19, 2025; detailed report published February 2026.",
    "English title says stock in 2025, while the model data terminal is 2024.",
    "The French title and body also identify 2024. Table 28 Sentinel mobile indices",
    "and Table 24 commercial catch weights are recovered. Other processed survey",
    "indices, annual stock weights and revised maturity ogives remain unavailable."
  ),
  stringsAsFactors = FALSE
)

assumptions <- data.frame(
  assessment_id = assessment_id,
  component = c(
    "population", "population", "population", "population", "recruitment",
    "recruitment", "F", "F", "M", "M", "catch", "catch",
    rep("index", 6), "biology", "biology", "biology"
  ),
  fleet = "",
  survey = c(
    rep("", 12), "DFO August", "Sentinel mobile", "Sentinel gillnet",
    "Sentinel summer longline", "Sentinel fall longline", "Minet bottom trawl",
    rep("", 3)
  ),
  sex = "",
  region = "",
  season = "",
  setting = c(
    "modeled_years", "modeled_ages", "recruitment_age", "plus_group",
    "mean_recruitment_periods", "process", "mean_structure", "process",
    "fixed_ages", "process", "catch_age_input", "catch_uncertainty",
    rep("series", 6), "stock_weights", "maturity", "maturity_age_11_plus"
  ),
  value = c(
    "1973-2024", "2-11+", "2", "11+", "1973-1990; 1991-2024",
    "Log recruitment deviations have a temporal correlation; Table 34 reports sigma and phi.",
    "Sixteen mean-F values defined by age and year blocks.",
    "Stochastic deviations correlated across age and year; variance differs for ages 2-3, 4, 5 and 6+.",
    "M at ages 2 and 3 is fixed at 1.0 and 0.65; ages 4+ are estimated from 1984.",
    "M deviations are correlated across age and year; Table 34 reports shared process parameters.",
    "One published catch-at-age matrix covers ages 2-11+ from 1974 to 2024.",
    "The assessment uses censored catch observations for some years from 2006 onward.",
    "DFO August survey; 1985-2024, ages 2-11+.",
    "Sentinel mobile survey; 1995-2024, ages 2-11+.",
    "Sentinel gillnet survey; 1995-2024, ages 4-11+.",
    "Sentinel summer longline (LLS1); 1995-2024, ages 3-11+.",
    "Sentinel fall longline (LLS2); 1995-2020, ages 3-11+.",
    "Minet bottom-trawl survey; 1973-1976, ages 3-11+.",
    "Beginning-of-year stock weights are estimated annually; numerical age-year values were not found in the cached report tables.",
    "Female maturity ogives are used for SSB; numerical age-year values were not found in the cached report tables.",
    "The 11+ mature proportion is an equilibrium-abundance-weighted average of ages 11-15."
  ),
  source_reference = c(
    rep(paste(assessment_ref, "Table 33"), 4),
    paste(assessment_ref, "Figure 95 and Table 34"),
    paste(assessment_ref, "Figure 95 and Table 34"),
    paste(assessment_ref, "Section 2.4 and Figure 18"),
    paste(assessment_ref, "Table 34"),
    paste(assessment_ref, "Table 36"),
    paste(assessment_ref, "Table 34"),
    paste(assessment_ref, "Table 30"),
    paste(assessment_ref, "Table 34"),
    rep(paste(assessment_ref, "Section 2.4.1"), 6),
    paste(assessment_ref, "Section 2.3"),
    paste(assessment_ref, "Section 2.2.3.1"),
    paste(assessment_ref, "Section 2.2.3.1")
  ),
  notes = rep("", 21)
)

assumptions$source_reference[assumptions$survey == "Sentinel mobile"] <-
  paste(assessment_ref, "Sections 2.1.3.3 and 2.4.1, Table 28")
assumptions$notes[assumptions$survey == "Sentinel mobile"] <-
  "Table 28 also reports age 1; the accepted model uses ages 2-11+. The standardized series is retained as one survey."
assumptions$notes[assumptions$setting == "maturity"] <- paste(
  "The revised annual ogives are fitted using a cohort-effect beta-binomial model.",
  "Technical Report 3671 (2025), Figures 5 and 7, provides no numeric annual matrix;",
  "raw maturity-stage samples are not the accepted fitted ogives."
)
add_assumption <- function(component, survey, setting, value, reference, note = "") {
  row <- assumptions[1, , drop = FALSE]
  row$component <- component
  row$survey <- survey
  row$setting <- setting
  row$value <- value
  row$source_reference <- reference
  row$notes <- note
  row
}
assumptions <- rbind(
  assumptions,
  add_assumption("index", "Sentinel mobile", "standardization",
    "July survey; 54 ft horizontal trawl opening and 1.25 nautical mile tow distance.",
    paste(assessment_ref, "Section 2.1.3.3"),
    "Reported indices already correct changing tow geometry; separate vessel q groups are not imposed on Table 28."),
  add_assumption("index", "DFO August", "vessel_and_gear",
    paste("Lady Hammond / Western IIA (1984-1990); Alfred Needler / URI (1990-2004);",
          "Teleost / Campelen (2004-2022); John Cabot / modified Campelen",
          "(comparative fishing 2021-2022; regular survey from 2023)."),
    paste(assessment_ref, "Section 2.1.3.1"),
    "Published accepted indices are calibrated to John Cabot equivalents. Uncalibrated alternatives would require separate vessel/gear q groups, survey design weighting and age-length expansion."),
  add_assumption("biology", "", "catch_weights",
    "Commercial fishery weights at age 2-11+, 1974-2024, in kg per fish (Table 24).",
    paste(assessment_ref, "Table 24"),
    "Not substituted for stock weights; seven missing cells remain missing, and printed zero weights are retained.")
)

write_updated <- function(name, new_rows, key) {
  path <- file.path(database_dir, name)
  columns <- names(read.csv(path, nrows = 0, check.names = FALSE))
  new_rows <- new_rows[columns]
  id <- as.character(new_rows[[key]][1])
  lines <- readLines(path, warn = FALSE)
  starts_with_id <- function(line) {
    startsWith(line, paste0(id, ",")) || startsWith(line, paste0('"', id, '",'))
  }
  existing <- lines[vapply(lines, starts_with_id, logical(1))]

  csv_field <- function(value) {
    if (is.na(value)) return("")
    value <- enc2utf8(as.character(value))
    if (grepl('[",\r\n]', value) || grepl("^\\s|\\s$", value)) {
      paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
    } else {
      value
    }
  }
  new_lines <- vapply(seq_len(nrow(new_rows)), function(i) {
    paste(vapply(new_rows, function(column) csv_field(column[i]), character(1)),
          collapse = ",")
  }, character(1))

  if (length(existing)) {
    if (identical(existing, new_lines)) return(invisible(NULL))
    selected <- vapply(lines, starts_with_id, logical(1))
    insert_at <- which(selected)[1]
    replacement <- c(lines[seq_len(insert_at - 1L)], new_lines,
                     lines[seq_along(lines) > insert_at & !selected])
    writeLines(replacement, path, useBytes = TRUE)
    return(invisible(NULL))
  }

  original <- readBin(path, "raw", n = file.info(path)[["size"]])
  line_ending <- if (any(original[-1] == as.raw(10) &
                         original[-length(original)] == as.raw(13))) {
    "\r\n"
  } else {
    "\n"
  }
  prefix <- if (length(original) && tail(original, 1) != as.raw(10)) {
    line_ending
  } else {
    ""
  }
  appended <- paste0(prefix, paste(new_lines, collapse = line_ending), line_ending)
  connection <- file(path, open = "ab")
  on.exit(close(connection))
  writeBin(charToRaw(enc2utf8(appended)), connection)
}

write_updated("stocks.csv", stock, "stock_id")
write_updated("assessments.csv", assessment, "assessment_id")
write_updated("assumptions.csv", assumptions, "assessment_id")
write_updated("inputs.csv", inputs, "assessment_id")
write_updated("outputs.csv", outputs, "assessment_id")

cat("Imported 3Pn4RS catch-at-age, catch weights, Sentinel indices, fixed M, and available model outputs.\n")
