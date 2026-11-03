pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("dfo_herring_4vwx_2025_indicators", read_database())
inputs <- record$inputs; outputs <- record$outputs
stopifnot(
  record$stocks$charbonneau_id == "DFO-MAR_harengus_4VWX",
  record$stocks$scientific_name == "Clupea harengus",
  record$stocks$ocean == "Northwest Atlantic",
  record$assessment$assessment_year == 2025,
  record$assessment$terminal_year == 2024,
  !record$assessment$is_current, record$assessment$is_applied,
  all(unlist(record$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(inputs) == 68, all(is.na(inputs$age)), all(is.na(inputs$sampling_time)),
  !any(inputs$measure %in% c("numbers_at_age", "proportion_at_age")),
  all(grepl("4130293x.pdf", inputs$source_reference, fixed = TRUE)),
  !any(grepl("2018", inputs$source_reference, fixed = TRUE))
)
catch <- inputs[inputs$type == "catch", ]
sw <- catch[catch$region == "Southwest Nova Scotia/Bay of Fundy", ]
stopifnot(
  identical(sw$year, 2020:2024), identical(sw$value, c(34.4, 32.6, 23.9, 20.7, 15.7)),
  all(sw$unit == "kt"), all(sw$year_basis == "calendar_year"),
  !any(grepl("New Brunswick", catch$region, fixed = TRUE)),
  all(catch$year_basis[catch$year < 2020] == "model_year"),
  all(catch$value[catch$region == "Bras d'Or Lakes"] == 0),
  catch$value[catch$region == "Offshore Scotian Shelf"] == 10,
  all(catch$source_type[catch$region == "Offshore Scotian Shelf"] == "official_document")
)
index <- inputs[inputs$type == "index", ]
little_hope <- index[index$region == "Little Hope/Port Mouton", ]
halifax <- index[index$region == "Halifax/Eastern Shore", ]
stopifnot(
  all(index$measure == "relative_biomass_index"), all(index$basis == "index_scale"),
  all(index$unit == "tonnes"), all(index$year_basis == "calendar_year"),
  identical(little_hope$value, c(145395, 61408, 66815, 168164, 92019, 35739, 82297, 41383, 68573, 64189)),
  identical(halifax$value, c(68562, 54312, 58681, 42416, 141198, 26205, 18341, 20313, 28057, 15430)),
  round(mean(little_hope$value[little_hope$year >= 2020])) == 58436,
  round(mean(halifax$value[halifax$year >= 2020])) == 21669,
  index$value[index$region == "German Bank and Scots Bay"] == 301290,
  sum(index$region == "Glace Bay") == 1,
  index$value[index$region == "Glace Bay"] == 50,
  nrow(outputs) == 1, outputs$year == 2024, outputs$value == 296430,
  outputs$measure == "acoustic_SSB_index_three_year_mean", outputs$unit == "tonnes",
  all(is.na(outputs$se)), all(is.na(outputs$lwr)), all(is.na(outputs$upr))
)
reference <- database_to_tam_ref(record$assessment$assessment_id, outputs,
  years = 2020:2024, ages = 1:11, terminal_year = 2024, assumptions = record$assumptions)
stopifnot(all(is.na(reference$pop$ssb$est)), all(is.na(reference$pop$N$est)),
  all(is.na(reference$pop$F$est)))
blocked <- run_assessment(record$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("4VWX herring indicator values, units, component separation and absence of fitted states passed.\n")
