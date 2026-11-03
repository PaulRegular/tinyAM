pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

x <- run_assessment("nefsc_scup_2015", fit = FALSE, silent = TRUE)
comparison_outputs <- x$comparison$comparison_outputs
derived <- comparison_outputs[
  !is.na(comparison_outputs$source_type) &
    comparison_outputs$source_type == "derived_source_output", , drop = FALSE
]
biomass <- derived[derived$measure == "biomass_at_age", , drop = FALSE]
ssb <- derived[derived$measure == "mature_biomass_at_age", , drop = FALSE]

source_n <- x$source$outputs[
  !is.na(x$source$outputs$measure) & x$source$outputs$measure == "numbers_at_age" &
    !is.na(x$source$outputs$year) & as.integer(x$source$outputs$year) == 2014 &
    !is.na(x$source$outputs$age) &
    as.character(x$source$outputs$age) == "3", , drop = FALSE
]
source_f <- x$source$outputs[
  !is.na(x$source$outputs$measure) & x$source$outputs$measure == "fishing_mortality_at_age" &
    !is.na(x$source$outputs$year) & as.integer(x$source$outputs$year) == 2014 &
    !is.na(x$source$outputs$age) &
    as.character(x$source$outputs$age) == "3", , drop = FALSE
]
source_m <- x$source$inputs[
  !is.na(x$source$inputs$type) & x$source$inputs$type == "M" &
    !is.na(x$source$inputs$age) & as.character(x$source$inputs$age) == "3", , drop = FALSE
]
weight <- x$obs$weight[
  !is.na(x$obs$weight$year) & as.integer(x$obs$weight$year) == 2014 &
    !is.na(x$obs$weight$age) & as.character(x$obs$weight$age) == "3", , drop = FALSE
]
maturity <- x$obs$maturity[
  !is.na(x$obs$maturity$year) & as.integer(x$obs$maturity$year) == 2014 &
    !is.na(x$obs$maturity$age) & as.character(x$obs$maturity$age) == "3", , drop = FALSE
]
expected_biomass <- as.numeric(source_n$value) * weight$obs
expected_ssb <- expected_biomass * maturity$obs *
  exp(-(as.numeric(source_f$value) + as.numeric(source_m$value)) * 0.5)
biomass_age <- as.integer(as.character(biomass$age))
ssb_age <- as.integer(as.character(ssb$age))
biomass_check <- biomass$value[!is.na(biomass$year) & biomass$year == 2014 & biomass_age == 3]
ssb_check <- ssb$value[!is.na(ssb$year) & ssb$year == 2014 & ssb_age == 3]

stopifnot(
  identical(x$diagnostics$status, "not_fitted"),
  is.null(x$fit),
  isTRUE(tinyAM::check_obs(x$obs)),
  nrow(biomass) == 31L * 8L,
  nrow(ssb) == 31L * 8L,
  identical(sort(unique(as.integer(biomass$year))), 1984:2014),
  identical(sort(unique(biomass_age)), 0:7),
  all(biomass$unit == "metric tons"),
  all(ssb$unit == "metric tons"),
  abs(as.numeric(biomass_check) - expected_biomass) < 1e-6 * expected_biomass,
  abs(as.numeric(ssb_check) - expected_ssb) < 1e-6 * expected_ssb,
  x$comparison$comparison_definitions$ssb$status == "approximate",
  any(grepl("directly reported aggregate SSB", x$background, ignore.case = TRUE))
)

message("Scup fit-disabled translation and common-definition biomass/SSB derivations verified.")
