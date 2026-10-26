pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
s <- read_assessment("ices_plaice_celtic_sea_7fg_2026")
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$stock_id == "ck_ices_wgsce_platessa_cs7f_g",
  s$stocks$charbonneau_id == "ICES-WGSCE_platessa_CS7f-g",
  s$stocks$authority_stock_id == "ple.27.7fg",
  s$stocks$scientific_name == "Pleuronectes platessa", s$stocks$ocean == "Northeast Atlantic",
  s$assessment$assessment_year == 2026, s$assessment$terminal_year == 2025,
  s$assessment$is_current, s$assessment$is_applied, s$assessment$framework_year == 2022,
  all(unlist(s$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(i) == 685L, nrow(o) == 46L, nrow(s$assumptions) == 11L,
  all(i$year_basis == "calendar_year"), all(is.na(i$sampling_time)),
  !any(i$type %in% c("M", "weight", "maturity")))
catch <- i[i$type == "catch", ]
landings <- catch[catch$fleet == "Combined landings", ]
discards <- catch[catch$fleet == "Combined discards", ]
stopifnot(nrow(catch) == 76L, all(is.na(catch$age)),
  all(catch$measure == "total_biomass"), all(catch$unit == "tonnes"),
  identical(landings$year, 1988:2025), tail(landings$value, 1) == 107,
  sum(is.na(discards$value)) == 16L, all(is.na(discards$value[discards$year < 2004])),
  discards$value[discards$year == 2012] == 1123,
  tail(discards$value, 1) == 180,
  tail(landings$value + discards$value, 1) == 287,
  all(grepl("Table 8", catch$source_reference, fixed = TRUE)))
biomass <- i[i$measure == "relative_biomass_index", ]
uk <- i[i$measure == "numbers_at_age" & i$survey == "UK(E&W)-BTS-Q3", ]
irish <- i[i$measure == "numbers_at_age" & i$survey == "IGFS-WIBTS-Q4", ]
effort <- i[i$measure == "effort", ]
stopifnot(identical(biomass$year, 1988:2025), all(biomass$unit == "kg per 100 km towed"),
  all(biomass$basis == "index_scale"), tail(biomass$value, 1) == 44,
  nrow(uk) == 33L * 11L, setequal(uk$year, 1993:2025), setequal(uk$age, 0:10),
  all(uk$unit == "fish per 100 km towed"), all(uk$basis == "numbers"),
  all(uk$season == "Q3"), all(uk$region == "Division 7.f"),
  uk$value[uk$year == 2025 & uk$age == 2] == 127,
  uk$value[uk$year == 2025 & uk$age == 10] == 0,
  nrow(irish) == 23L * 7L, setequal(irish$year, 2003:2025), setequal(irish$age, 1:7),
  all(irish$basis == "index_scale"), all(irish$season == "Q4"),
  all(grepl("normalization unresolved", irish$unit, fixed = TRUE)),
  irish$value[irish$year == 2025 & irish$age == 2] == 108,
  !any(irish$age == 8), !any(grepl("age 7 the plus", irish$notes, fixed = TRUE)),
  all(effort$survey == "IGFS-WIBTS-Q4"), all(effort$unit == "fished minutes"),
  identical(effort$year, 2003:2025), tail(effort$value, 1) == 1530,
  all(grepl("Ancillary", c(uk$notes, irish$notes), fixed = TRUE)))
pressure <- o[o$measure == "length_based_fishing_pressure_proxy", ]
lengths <- i[i$type == "biology", ]
mean_length <- lengths[lengths$measure == "mean_catch_length", ]
stopifnot(nrow(lengths) == 24L, all(lengths$unit == "cm"), all(is.na(lengths$age)),
  identical(mean_length$value, c(29.2, 29.8, 28.4, 27.7, 29.0, 27.3, 28.3, 27.7)),
  all(lengths$value[lengths$measure == "MSY_proxy_length"] == 30),
  all(abs(30 / mean_length$value - pressure$value) < .006))
stopifnot(identical(pressure$year, 2018:2025),
  identical(pressure$value, c(1.03, 1.01, 1.06, 1.08, 1.03, 1.10, 1.06, 1.08)),
  all(pressure$unit == "ratio"),
  identical(o$value[o$measure == "survey_biomass_index"], biomass$value),
  !any(o$measure %in% c("SSB", "numbers_at_age", "Fbar", "fishing_mortality_at_age", "recruitment")),
  all(is.na(o$se)), all(is.na(o$lwr)), all(is.na(o$upr)),
  any(grepl("ass_yr=2024", s$assumptions$notes, fixed = TRUE)),
  any(grepl("r=0.55", s$assumptions$notes, fixed = TRUE)))
reference <- database_to_tam_ref(s$assessment$assessment_id, o,
  years = 1988:2025, ages = 0:10, terminal_year = 2025, assumptions = s$assumptions)
stopifnot(all(is.na(reference$pop$ssb$est)), all(is.na(reference$pop$N$est)),
  all(is.na(reference$pop$F$est)), all(is.na(reference$pop$recruitment$est)))
blocked <- run_assessment(s$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("Celtic Sea plaice source identity, published indicators, units, age-table ambiguities and missing fitted states passed.\n")
