pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
s <- read_assessment("ices_golden_redfish_northeast_arctic_2026")
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$stock_id == "ck_ices_afwg_norvegicus_deep1_2",
  s$stocks$charbonneau_id == "ICES-AFWG_norvegicus_DEEP1,2",
  s$stocks$authority_stock_id == "reg.27.1-2", s$stocks$scientific_name == "Sebastes norvegicus",
  s$stocks$ocean == "Northeast Atlantic", s$assessment$model_family == "SPiCT",
  s$assessment$assessment_year == 2026, s$assessment$framework_year == 2026,
  s$assessment$terminal_year == 2025, s$assessment$estimate_terminal_year == 2026,
  s$assessment$is_current, s$assessment$is_applied,
  all(unlist(s$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(i) == 118L, nrow(o) == 237L, nrow(s$assumptions) == 14L)
stopifnot(all(i$type == "catch"), all(i$measure == "total_biomass"), all(i$basis == "biomass"),
  identical(i$year, 1908:2025), all(i$year_basis == "calendar_year"),
  all(i$unit == "tonnes"), all(is.na(i$age)), all(is.na(i$sampling_time)),
  length(unique(i$fleet)) == 1L, all(i$region == "ICES subareas 1 and 2"),
  all(i$source_type == "official_machine_readable"),
  identical(i$value[i$year %in% c(1908, 1945, 1957, 2015, 2025)], c(650, 560, 51610, 3917, 11435)),
  grepl("2069 t equals 2024", i$notes[i$year == 2025], fixed = TRUE),
  !any(i$type %in% c("index", "weight", "M", "maturity")))
b <- o[o$measure == "B_over_Bmsy", ]
f <- o[o$measure == "F_over_Fmsy", ]
stopifnot(identical(b$year, 1908:2026), identical(f$year, 1908:2025),
  all(o$unit == "ratio"), all(is.na(o$age)), all(is.na(o$se)),
  all(o$lwr <= o$value), all(o$value <= o$upr),
  identical(c(tail(b$value, 1), tail(b$lwr, 1), tail(b$upr, 1)), c(.44, .26, .75)),
  identical(c(tail(f$value, 1), tail(f$lwr, 1), tail(f$upr, 1)), c(1.59, .91, 2.76)),
  all(grepl("does not state their confidence level", o$notes, fixed = TRUE)),
  !any(o$measure %in% c("SSB", "total_biomass", "Fbar", "fishing_mortality_at_age", "recruitment", "numbers_at_age")),
  any(grepl("1.59", s$assumptions$value, fixed = TRUE)),
  any(grepl("1.71", s$assumptions$value, fixed = TRUE)),
  any(grepl("alternative high-biomass", s$assumptions$value, fixed = TRUE)))
reference <- database_to_tam_ref(s$assessment$assessment_id, o,
  years = 1908:2026, ages = 1:2, terminal_year = 2025, assumptions = s$assumptions)
stopifnot(all(is.na(reference$pop$ssb$est)), all(is.na(reference$pop$F$est)),
  all(is.na(reference$pop$N$est)), all(is.na(reference$pop$recruitment$est)),
  identical(reference$dat$is_proj, c(rep(FALSE, 118), TRUE)))
blocked <- run_assessment(s$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("Golden redfish identity, catch history, relative-state units/bounds, reported horizon and absent age-based comparison passed.\n")
