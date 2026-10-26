.libPaths(c("C:/Users/paulr/AppData/Local/R/win-library/4.4", .libPaths()))
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
record <- read_assessment("asmfc_striped_bass_atlantic_2018", db)
years <- 1982:2017
inputs <- record$inputs
outputs <- record$outputs
stopifnot(
  record$stocks$authority == "ASMFC",
  record$stocks$ocean == "Northwest Atlantic",
  !record$assessment$is_current,
  record$assessment$is_applied,
  record$assessment$terminal_year == 2017,
  record$assessment$assessment_type == "benchmark",
  all(unlist(record$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(inputs) == 6657L,
  nrow(outputs) == 1764L,
  all(inputs$source_type %in% c("official_table")),
  all(outputs$source_type == "official_table")
)
cell <- function(x, year, age, measure, type = NULL, fleet = NULL, survey = NULL) {
  keep <- x$measure == measure & x$year == year
  keep <- keep & if (is.na(age)) is.na(x$age) else x$age == age
  if (!is.null(type)) keep <- keep & x$type == type
  if (!is.null(fleet)) keep <- keep & x$fleet == fleet
  if (!is.null(survey)) keep <- keep & x$survey == survey
  z <- x[keep, , drop = FALSE]
  stopifnot(nrow(z) == 1L)
  z$value[[1L]]
}
catch <- inputs[inputs$type == "catch" & inputs$measure == "proportion_at_age", ]
catch_sums <- aggregate(value ~ year + fleet, catch, sum)
stopifnot(
  nrow(catch) == 1080L,
  max(abs(catch_sums$value - 1)) < 0.002,
  abs(cell(inputs, 1982, 2, "proportion_at_age", "catch", "Chesapeake Bay") - 0.19419) < 1e-8,
  cell(inputs, 1982, NA, "total_numbers", "catch", "Chesapeake Bay") == 228642
)
stopifnot(
  nrow(inputs[inputs$measure == "weight_at_age", ]) == 540L,
  nrow(inputs[inputs$measure == "spawning_weight_at_age", ]) == 540L,
  nrow(inputs[inputs$measure == "catch_weight_at_age", ]) == 540L,
  cell(inputs, 1982, 1, "weight_at_age", "weight") == 0.06,
  cell(inputs, 1982, 4, "spawning_weight_at_age", "weight") == 1.38,
  cell(inputs, 1982, 1, "catch_weight_at_age", "catch_weight") == 0.13
)
linked <- inputs[inputs$type == "index" & inputs$measure == "total_numbers" & !is.na(inputs$age), ]
stopifnot(
  all(linked$year >= 1982 & linked$year <= 2017),
  all(linked$age == ifelse(linked$survey %in% c("NY Age-1", "MD Age-1"), 2, 1)),
  all(inputs$year[inputs$type == "index" & !is.na(inputs$year)] >= 1982),
  abs(cell(inputs, 1985, 2, "total_numbers", "index", survey = "NY Age-1") -
    inputs$value[inputs$type == "index" & inputs$measure == "total_numbers" &
      inputs$survey == "NY Age-1" & inputs$year == 1985]) < 1e-8
)
de30 <- inputs[inputs$type == "index" & inputs$measure == "proportion_at_age" &
  inputs$survey == "DE 30 Trawl" & inputs$year == 2017, ]
stopifnot(abs(sum(de30$value) - 1.105) < 1e-8)
N <- outputs[outputs$type == "population" & outputs$measure == "numbers_at_age", ]
F <- outputs[outputs$type == "mortality" & outputs$measure == "fishing_mortality_at_age", ]
SSB <- outputs[outputs$type == "biomass" & outputs$measure == "SSB", ]
SSB_age <- outputs[outputs$type == "biomass" & outputs$measure == "mature_biomass_at_age", ]
stopifnot(
  nrow(N) == 540L, nrow(F) == 540L, nrow(SSB_age) == 540L,
  cell(N, 1982, 1, "numbers_at_age", "population") == 37879200,
  cell(F, 1982, 15, "fishing_mortality_at_age", "mortality") == 0.171,
  SSB$value[SSB$year == 2017] == 68476,
  SSB$se[SSB$year == 2017] > 0,
  all(SSB$unit == "metric tons"),
  all(abs(aggregate(value ~ year, SSB_age, sum)$value - SSB$value) <= 8)
)
translated <- run_assessment(record$assessment$assessment_id, database = db, fit = FALSE, silent = TRUE)
comparison_outputs <- translated$comparison$comparison_outputs
total_biomass <- comparison_outputs[
  comparison_outputs$type == "biomass" & comparison_outputs$measure == "total_biomass", ]
stopifnot(
  translated$diagnostics$status == "not_fitted",
  nrow(translated$obs$catch) == 540L,
  nrow(translated$obs$weight) == 540L,
  nrow(translated$obs$maturity) == 540L,
  nrow(total_biomass) == length(years),
  !anyDuplicated(total_biomass$year),
  setequal(total_biomass$year, years),
  tinyAM::check_obs(translated$obs)
)
cat("2018 striped-bass source extraction and tinyAM observation preparation passed.\n")


