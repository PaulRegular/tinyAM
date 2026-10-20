source("analysis/comp_assessments/R/read_database.R")

x <- read_assessment("dfo_american_plaice_4t_2016")
stock <- x$stocks[x$stocks$stock_id == "ck_dfo_platessoides_4t", , drop = FALSE]
events <- x$assessments[x$assessments$stock_id == stock$stock_id, , drop = FALSE]
i <- x$inputs
o <- x$outputs

stopifnot(
  nrow(stock) == 1L,
  stock$authority_stock_id == "4T",
  stock$scientific_name == "Hippoglossoides platessoides",
  stock$common_name == "American plaice",
  nrow(events) == 2L,
  events$is_current[events$assessment_id == "dfo_american_plaice_4t_2016"],
  !events$is_current[events$assessment_id == "dfo_american_plaice_4t_2024_indicators"],
  events$inputs_status[events$assessment_id == "dfo_american_plaice_4t_2016"] == "complete",
  events$outputs_status[events$assessment_id == "dfo_american_plaice_4t_2016"] == "partial"
)

catch <- i[i$type == "catch" & i$measure == "numbers_at_age", , drop = FALSE]
rv <- i[!is.na(i$survey) & i$survey == "DFO Research Vessel Trawl Survey", , drop = FALSE]
sentinel <- i[!is.na(i$survey) & i$survey == "Mobile Gear Sentinel Survey", , drop = FALSE]
weight <- i[i$type == "weight" & i$measure == "weight_at_age", , drop = FALSE]
maturity <- i[i$type == "maturity" & i$measure == "maturity_at_age", , drop = FALSE]
stopifnot(
  nrow(catch) == 40L * 13L,
  identical(sort(unique(catch$year)), 1976:2015),
  identical(sort(unique(catch$age)), 4:16),
  all(catch$unit == "thousand fish"),
  all(catch$fleet == "Commercial fishery"),
  nrow(rv) == 39L * 13L,
  !any(rv$year == 2003),
  all(rv$sampling_time == 0.75),
  all(rv$unit == "thousand fish"),
  nrow(sentinel) == 13L * 13L,
  identical(sort(unique(sentinel$year)), 2003:2015),
  all(sentinel$sampling_time == 0.75),
  nrow(weight) == 40L * 13L,
  weight$value[weight$year == 2003 & weight$age == 16] == 0,
  weight$value[weight$year == 1976 & weight$age == 4] == 0.047,
  nrow(maturity) == 40L * 13L,
  all(maturity$value >= 0 & maturity$value <= 1),
  maturity$value[maturity$year == 1976 & maturity$age == 4] == 0.350,
  identical(
    maturity$value[maturity$year == 2012][order(maturity$age[maturity$year == 2012])],
    maturity$value[maturity$year == 2015][order(maturity$age[maturity$year == 2015])]
  )
)

n <- o[o$measure == "numbers_at_age", , drop = FALSE]
total <- o[o$measure == "total_numbers", , drop = FALSE]
recruitment <- o[o$measure == "recruitment", , drop = FALSE]
ssb <- o[o$measure == "SSB", , drop = FALSE]
stopifnot(
  nrow(n) == 40L * 13L,
  n$value[n$year == 2015 & n$age == 16] == 804,
  nrow(total) == 40L,
  nrow(recruitment) == 40L,
  all(recruitment$type == "recruitment"),
  recruitment$value[recruitment$year == 2015] == 528392,
  total$value[total$year == 2015] == sum(n$value[n$year == 2015]),
  all(is.na(total$se)), all(is.na(recruitment$se)),
  nrow(ssb) == 1L,
  ssb$year == 2015,
  ssb$value == 55023,
  all(!o$measure %in% c("fishing_mortality_at_age", "natural_mortality_at_age"))
)

indicator <- read_assessment("dfo_american_plaice_4t_2024_indicators")
stopifnot(
  nrow(indicator$outputs) == 1L,
  indicator$outputs$measure == "relative_biomass_index",
  indicator$outputs$year == 2024,
  indicator$outputs$value == 1.25,
  indicator$outputs$unit == "kg/tow"
)

cat("4T American plaice source coverage, provenance, and output checks passed.\n")
