pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
s <- read_assessment("ices_plaice_western_channel_7e_2025")
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$stock_id == "ck_ices_wgcse_platessa_is7e",
  s$stocks$charbonneau_id == "ICES-WGCSE_platessa_IS7e", s$stocks$authority_stock_id == "ple.27.7e",
  s$stocks$scientific_name == "Pleuronectes platessa", s$stocks$ocean == "Northeast Atlantic",
  s$assessment$assessment_year == 2025, s$assessment$terminal_year == 2024,
  s$assessment$estimate_terminal_year == 2024, s$assessment$framework_year == 2024,
  s$assessment$is_current, s$assessment$is_applied,
  all(unlist(s$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(i) == 134L, nrow(o) == 66L, nrow(s$assumptions) == 13L,
  all(is.na(i$age)), all(is.na(i$sampling_time)), max(i$year) == 2024,
  !any(i$type %in% c("M", "weight", "maturity")))
l <- i[i$measure == "landings_biomass", ]
d <- i[i$measure == "discard_biomass", ]
c <- i[i$type == "catch" & i$measure == "total_biomass", ]
index <- i[i$type == "index", ]
stopifnot(identical(l$year, 1980:2024), identical(d$year, 1980:2024),
  sum(is.na(d$value)) == 22L, all(is.na(d$value[d$year < 2002])),
  all(is.finite(d$value[d$year >= 2002])),
  identical(c$year, 2003:2024), tail(c$value, 1) == 1336,
  all(i$unit[i$type == "catch"] == "tonnes"),
  length(unique(i$fleet[i$type == "catch"])) == 1L,
  all(grepl("migration", i$notes[i$type == "catch"], fixed = TRUE)),
  identical(index$year, 2003:2024), all(index$measure == "relative_biomass_index"),
  all(index$unit == "kg per hour per metre beam"), all(index$season == "Q3"),
  all(index$survey == "UK-FSP-Q3 (B4381)"),
  all(grepl("ages 2-8", index$notes, fixed = TRUE)), abs(tail(index$value, 1) - .80) <= .005)
hr <- o[o$measure == "relative_harvest_rate", ]
dead <- o[o$measure == "assumed_dead_catch", ]
biomass <- o[o$measure == "survey_biomass_index", ]
stopifnot(identical(hr$year, 2003:2024), tail(hr$value, 1) == 1291,
  all(hr$unit == "tonnes per (kg per hour per metre beam)"),
  identical(dead$year, 2003:2024), tail(dead$value, 1) == 1030,
  tail(dead$value, 1) < tail(c$value, 1),
  identical(biomass$value, index$value),
  all(is.na(o$se)), all(is.na(o$lwr)), all(is.na(o$upr)),
  !any(o$measure %in% c("SSB", "total_biomass", "Fbar", "fishing_mortality_at_age", "numbers_at_age", "recruitment")),
  any(grepl("50%", s$assumptions$value, fixed = TRUE)),
  any(grepl("no new advice", s$assumptions$value, fixed = TRUE)),
  any(grepl("fca97de", s$assessment$model_version, fixed = TRUE)))
ref <- database_to_tam_ref(s$assessment$assessment_id, o,
  years = 1980:2025, ages = 2:8, terminal_year = 2024, assumptions = s$assumptions)
stopifnot(all(is.na(ref$pop$ssb$est)), all(is.na(ref$pop$F$est)), all(is.na(ref$pop$N$est)),
  all(is.na(ref$pop$recruitment$est)), identical(ref$dat$is_proj, c(rep(FALSE, 45), TRUE)))
blocked <- run_assessment(s$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("Western Channel plaice identity, catch components/survival, indicator units, production vintage and unavailable age-based reference passed.\n")
