pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("nefsc_yellowtail_georges_bank_2025", read_database())
inputs <- record$inputs
outputs <- record$outputs
stopifnot(
  record$stocks$charbonneau_id == "DFO-NEFSC_ferruginea_GB",
  record$stocks$ocean == "Northwest Atlantic",
  record$stocks$scientific_name == "Limanda ferruginea",
  record$assessment$assessment_year == 2025,
  record$assessment$terminal_year == 2024,
  record$assessment$estimate_terminal_year == 2028,
  !record$assessment$is_current, record$assessment$is_applied,
  all(unlist(record$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(inputs) == 10, identical(inputs$year, 2015:2024),
  all(inputs$measure == "total_biomass"), all(inputs$unit == "tonnes"),
  all(inputs$source_type == "official_table"),
  all(outputs$source_type[!grepl("rho_adjusted", outputs$measure)] == "official_table"),
  all(outputs$source_type[grepl("rho_adjusted", outputs$measure)] == "official_document"),
  all(is.na(inputs$age)), inputs$value[inputs$year == 2024] == 12,
  nrow(outputs) == 40,
  !any(outputs$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "Fbar", "natural_mortality_at_age")),
  all(is.na(outputs$se)),
  !any(grepl("Constant_M|constantM|fit_input_2024", c(inputs$source_reference, outputs$source_reference))),
  any(grepl("updated", record$assumptions$value[record$assumptions$setting == "maturity_revision"]))
)
select_output <- function(measure) outputs[outputs$measure == measure, ]
ssb <- select_output("SSB")
f <- select_output("Ffull")
recruitment <- select_output("recruitment")
historical <- ssb$year <= record$assessment$terminal_year
stopifnot(
  identical(ssb$value[historical], c(3324, 2330, 1401, 910, 747, 445, 841, 614, 633, 713)),
  identical(f$value[f$year <= 2024], c(0.088, 0.03, 0.082, 0.063, 0.011, 0.216, 0.068, 0.023, 0.039, 0.017)),
  identical(recruitment$value, c(5143, 2514, 2524, 4495, 8455, 1742, 519, 1387, 508, 534)),
  all(recruitment$age == 1), all(recruitment$unit == "thousand fish"),
  all(f$age_group == "4-6"), all(is.na(f$age)),
  ssb$lwr[ssb$year == 2024] == 491, ssb$upr[ssb$year == 2024] == 1034,
  f$lwr[f$year == 2024] == 0.012, f$upr[f$year == 2024] == 0.025,
  select_output("SSB_rho_adjusted")$value == 684,
  select_output("SSB_rho_adjusted")$year == 2024,
  identical(ssb$year[!historical], 2025:2028),
  all(grepl("Projection, not a fitted", ssb$notes[!historical], fixed = TRUE)),
  all(grepl("Interval level and construction are not stated", ssb$notes[!historical], fixed = TRUE)),
  all(is.na(ssb$lwr[ssb$year < 2024])),
  all(grepl("90% confidence region", ssb$notes[historical], fixed = TRUE))
)
reference <- database_to_tam_ref(record$assessment$assessment_id, outputs,
  years = 2015:2028, ages = 1:6, terminal_year = 2024,
  assumptions = record$assumptions)
stopifnot(
  reference$pop$ssb$est[reference$pop$ssb$year == 2024] == 713,
  all(reference$pop$ssb$is_proj == (reference$pop$ssb$year > 2024)),
  all(is.na(reference$pop$F$est)), all(is.na(reference$pop$F_bar$est)),
  all(is.na(reference$pop$N$est)),
  !any(reference$pop$ssb$est == 684),
  all(is.na(reference$pop$ssb$se))
)
blocked <- run_assessment(record$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("Georges Bank yellowtail 2025 summary values, run separation and uncertainty passed.\n")
