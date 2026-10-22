source("analysis/comp_assessments/R/read_database.R")
source_data <- read_assessment("ices_herring_gulf_bothnia_2026")

inputs <- source_data$inputs
outputs <- source_data$outputs
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", ]
catch_prop <- inputs[inputs$type == "catch" & inputs$measure == "proportion_at_age", ]
catch_total <- inputs[inputs$type == "catch" & inputs$measure == "total_biomass", ]
catch_weight <- inputs[inputs$type == "catch_weight" & inputs$measure == "weight_at_age", ]
maturity <- inputs[inputs$type == "maturity" & inputs$measure == "maturity_at_age", ]
m <- inputs[inputs$type == "M" & inputs$measure == "natural_mortality_at_age", ]
acoustic <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age", ]
recruitment <- outputs[outputs$measure == "recruitment", ]
ssb <- outputs[outputs$measure == "SSB", ]
fbar <- outputs[outputs$measure == "Fbar", ]

proportion_sums <- aggregate(value ~ year, catch_prop, sum)

stopifnot(
  source_data$stocks$stock_id == "ices_herring_gulf_bothnia",
  source_data$stocks$authority_stock_id == "her.27.3031",
  source_data$stocks$area == "ICES subdivisions 30 and 31",
  source_data$stocks$ocean == "Northeast Atlantic",
  source_data$assessment$assessment_year == 2026,
  source_data$assessment$terminal_year == 2025,
  source_data$assessment$estimate_terminal_year == 2025,
  source_data$assessment$model_version == "3.30.22.1",
  source_data$assessment$is_current,
  source_data$assessment$is_applied,
  source_data$assessment$assumptions_status == "partial",
  source_data$assessment$inputs_status == "partial",
  source_data$assessment$outputs_status == "partial",
  nrow(catch_total) == 63L,
  identical(sort(unique(catch_total$year)), 1963:2025),
  catch_total$value[catch_total$year == 2025] == 67344,
  nrow(catch) == 46L * 15L,
  identical(sort(unique(catch$year)), 1980:2025),
  identical(sort(unique(catch$age)), 1:15),
  catch$value[catch$year == 2025 & catch$age == 1] == 426132,
  nrow(catch_prop) == 46L * 15L,
  all(abs(proportion_sums$value - 1) < 1e-12),
  all(catch_prop$basis == "proportion_numbers"),
  nrow(catch_weight) == 63L * 18L,
  identical(sort(unique(catch_weight$age)), 0:17),
  catch_weight$value[catch_weight$year == 2025 & catch_weight$age == 1] == 0.00885,
  nrow(maturity) == 46L * 17L,
  identical(sort(unique(maturity$age)), 1:17),
  maturity$value[maturity$year == 2025 & maturity$age == 2] == 0.73,
  nrow(m) == 21L,
  all(is.na(m$year)),
  m$value[m$age == 0] == 0.563,
  m$value[m$age == 15] == 0.257,
  nrow(acoustic) == 19L * 16L,
  identical(sort(unique(acoustic$year)), 2007:2025),
  identical(sort(unique(acoustic$age)), 0:15),
  acoustic$value[acoustic$year == 2007 & acoustic$age == 5] == 5449.22,
  acoustic$value[acoustic$year == 2025 & acoustic$age == 15] == 24.63,
  all(acoustic$unit == "million fish"),
  nrow(recruitment) == 63L,
  nrow(ssb) == 63L,
  nrow(fbar) == 63L,
  recruitment$value[recruitment$year == 2025] == 24418100,
  recruitment$lwr[recruitment$year == 2025] == 9692671.83925113,
  ssb$value[ssb$year == 2025] == 603383,
  ssb$upr[ssb$year == 2025] == 763375.042,
  fbar$value[fbar$year == 2025] == 0.138778,
  all(is.na(c(recruitment$se, ssb$se, fbar$se))),
  !any(outputs$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age")),
  all(inputs$year[!is.na(inputs$year)] <= 2025),
  all(outputs$year[!is.na(outputs$year)] <= 2025)
)

cat("Gulf of Bothnia herring source identity, coverage, units, transformations, and reported summaries passed.\n")
