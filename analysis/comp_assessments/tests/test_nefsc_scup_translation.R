pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

x <- run_assessment("nefsc_scup_2015", fit = FALSE, silent = TRUE)
ctdeep_name <- "CTDEEP fall trawl survey"
ctdeep <- x$obs$index[x$obs$index$survey == ctdeep_name, ]
raw_ctdeep <- x$source$inputs[x$source$inputs$survey %in% ctdeep_name, ]
native_plus <- raw_ctdeep[raw_ctdeep$measure %in% "index_by_age_group", ]
raw_plus <- raw_ctdeep[raw_ctdeep$measure %in% c("numbers_at_age", "index_by_age_group") &
  as.integer(raw_ctdeep$age) %in% 7:10, ]
expected_plus <- aggregate(as.numeric(raw_plus$value), list(year = raw_plus$year), sum)
expected_plus <- expected_plus[expected_plus$x > 0, ]
observed_plus <- ctdeep[ctdeep$age == 7L, ]

discard_name <- "Recreational discards"
raw_discard <- x$source$inputs[x$source$inputs$type %in% "catch" &
  x$source$inputs$measure %in% "numbers_at_age" &
  x$source$inputs$fleet %in% discard_name, ]
expected_discard <- aggregate(as.numeric(raw_discard$value) * 0.15,
  list(year = raw_discard$year, age = pmin(as.integer(raw_discard$age), 7L)), sum)
observed_discard <- x$obs$catch[x$obs$catch$fleet == discard_name, ]
key <- function(d) paste(d$year, d$age)

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
  identical(x$diagnostics$translation_tier, "curated"),
  is.null(x$fit),
  isTRUE(tinyAM::check_tam_data(x$obs, quiet = TRUE)),
  setequal(x$obs$catch$fleet, c("Commercial landings", "Commercial discards",
    "Recreational landings", discard_name)),
  nrow(x$obs$catch) == 31L * 8L * 4L,
  identical(observed_discard$unit, rep("thousand fish", nrow(observed_discard))),
  isTRUE(all.equal(observed_discard$obs,
    expected_discard$x[match(key(observed_discard), key(expected_discard))])),
  setequal(x$obs$index$survey, c("NEFSC fall trawl survey", ctdeep_name)),
  identical(sort(unique(as.integer(ctdeep$year))), setdiff(1984:2014, 2010L)),
  all(ctdeep$samp_time == 0.75), all(ctdeep$unit == "thousand fish per tow"),
  all(x$obs$index$samp_time[x$obs$index$survey == "NEFSC fall trawl survey"] == 0.8),
  nrow(native_plus) == 30L, all(as.integer(native_plus$age) == 10L),
  nrow(observed_plus) == 13L, 2014L %in% observed_plus$year,
  isTRUE(all.equal(observed_plus$obs, expected_plus$x[match(observed_plus$year,
    expected_plus$year)] / 1000)),
  !any(as.character(ctdeep$q_key) %in%
    as.character(x$obs$index$q_key[x$obs$index$survey != ctdeep_name])),
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

message("Curated Scup fleet, survey and common-definition biomass/SSB preparation verified.")
