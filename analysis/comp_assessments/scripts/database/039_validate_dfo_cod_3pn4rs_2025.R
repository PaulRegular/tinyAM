assessment_id <- "dfo_cod_3pn4rs_2025"
root <- file.path("analysis", "comp_assessments")
source_file <- file.path(root, "source_cache", assessment_id, "assessment.txt")
lines <- readLines(source_file, warn = FALSE)

read_table <- function(number, next_number, value_count) {
  start <- grep(paste0("Table ", number, "."), lines, fixed = TRUE)[[1]]
  end <- grep(paste0("Table ", next_number, "."), lines, fixed = TRUE)
  end <- end[end > start][[1]]
  rows <- lines[seq.int(start + 1L, end - 1L)]
  rows <- rows[grepl("^[[:space:]]*(19|20)[0-9]{2}[[:space:]]+", rows)]
  values <- strsplit(gsub(",", "", trimws(rows)), "[[:space:]]+")
  if (!length(values) || any(lengths(values) != value_count + 1L)) {
    stop("Unexpected source row structure in Table ", number, ".")
  }
  matrix(as.numeric(unlist(values)), ncol = value_count + 1L, byrow = TRUE)
}

tables <- list(
  catch = read_table(30, 31, 10),
  biomass = read_table(32, 33, 12),
  abundance = read_table(33, 34, 12),
  fishing_mortality = read_table(35, 36, 12),
  natural_mortality = read_table(36, 37, 10)
)
database <- file.path(root, "database")
inputs <- read.csv(file.path(database, "inputs.csv"), stringsAsFactors = FALSE,
                   check.names = FALSE)
outputs <- read.csv(file.path(database, "outputs.csv"), stringsAsFactors = FALSE,
                    check.names = FALSE)
inputs <- inputs[inputs$assessment_id == assessment_id, , drop = FALSE]
outputs <- outputs[outputs$assessment_id == assessment_id, , drop = FALSE]

same_values <- function(actual, expected) {
  isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = 0))
}

check_age_surface <- function(measure, source, columns) {
  actual <- outputs[outputs$measure == measure & !is.na(outputs$age), , drop = FALSE]
  actual <- actual[order(actual$year, actual$age), , drop = FALSE]
  expected <- as.vector(t(source[, columns, drop = FALSE]))
  nrow(actual) == length(expected) && same_values(actual$value, expected)
}

check_annual <- function(measure, source, column, age_group = "") {
  actual <- outputs[outputs$measure == measure & outputs$age_group == age_group,
                    , drop = FALSE]
  actual <- actual[order(actual$year), , drop = FALSE]
  nrow(actual) == nrow(source) &&
    same_values(actual$value, source[, column])
}

catch <- inputs[inputs$type == "catch", , drop = FALSE]
catch <- catch[order(catch$year, catch$age), , drop = FALSE]
fixed_m <- inputs[inputs$type == "M", , drop = FALSE]
m_estimated <- outputs[outputs$measure == "natural_mortality_at_age", , drop = FALSE]
m_source <- tables$natural_mortality[tables$natural_mortality[, 1] >= 1984,
                                      c(1, 4:11), drop = FALSE]
m_estimated <- m_estimated[order(m_estimated$year, m_estimated$age), , drop = FALSE]

stopifnot(
  all(tables$catch[, 1] == 1974:2024),
  all(tables$biomass[, 1] == 1973:2024),
  all(tables$abundance[, 1] == 1973:2024),
  all(tables$fishing_mortality[, 1] == 1973:2024),
  all(tables$natural_mortality[, 1] == 1973:2024),
  nrow(catch) == 510L,
  same_values(catch$value, as.vector(t(tables$catch[, 2:11, drop = FALSE]))),
  nrow(fixed_m) == 192L,
  all(fixed_m$value[fixed_m$age == 2] == 1),
  all(fixed_m$value[fixed_m$age == 3] == 0.65),
  nrow(m_estimated) == 328L,
  same_values(m_estimated$value, as.vector(t(m_source[, -1, drop = FALSE]))),
  check_age_surface("numbers_at_age", tables$abundance, 2:11),
  check_age_surface("biomass_at_age", tables$biomass, 2:11),
  check_age_surface("fishing_mortality_at_age", tables$fishing_mortality, 2:11),
  check_annual("SSB", tables$biomass, 13),
  check_annual("total_biomass", tables$biomass, 12, "2+"),
  check_annual("total_numbers", tables$abundance, 12, "2+"),
  check_annual("total_numbers", tables$abundance, 13, "5+"),
  check_annual("recruitment", tables$abundance, 2),
  check_annual("Fbar", tables$fishing_mortality, 12, "4-6"),
  check_annual("Fbar", tables$fishing_mortality, 13, "6-9")
)

cat("3Pn4RS source-table validation passed.\n")
