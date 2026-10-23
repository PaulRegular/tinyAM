source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))

assessment_id <- "ices_saithe_iceland_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
inputs <- source_data$inputs
outputs <- source_data$outputs
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", ]
index <- inputs[inputs$type == "index" & inputs$survey == "IS-SMB", ]
weight <- inputs[inputs$type == "weight", ]
maturity <- inputs[inputs$type == "maturity", ]
mortality <- inputs[inputs$type == "M", ]
n_at_age <- outputs[outputs$measure == "numbers_at_age", ]
f_at_age <- outputs[outputs$measure == "fishing_mortality_at_age", ]
ssb <- outputs[outputs$measure == "SSB", ]
reference_biomass <- outputs[outputs$measure == "biomass_by_age_group", ]
recruitment <- outputs[outputs$measure == "recruitment", ]
fbar <- outputs[outputs$measure == "Fbar", ]

stopifnot(
  nrow(source_data$stock) == 1L,
  source_data$stock$authority_stock_id == "pok.27.5a",
  source_data$stock$common_name == "Icelandic saithe",
  source_data$stock$ocean == "Northeast Atlantic",
  grepl("source_reviews/ices_saithe_iceland_2026.md", source_data$stock$notes, fixed = TRUE),
  nrow(source_data$assessment) == 1L,
  isTRUE(source_data$assessment$is_current),
  isTRUE(source_data$assessment$is_applied),
  source_data$assessment$assessment_year == 2026L,
  source_data$assessment$terminal_year == 2025L,
  source_data$assessment$estimate_terminal_year == 2026L,
  source_data$assessment$inputs_status == "complete",
  source_data$assessment$assumptions_status == "partial",
  source_data$assessment$outputs_status == "partial",
  nrow(catch) == 46L * 12L,
  setequal(unique(catch$year), 1980:2025),
  setequal(unique(catch$age), 3:14),
  all(catch$unit == "thousand fish"),
  nrow(index) == 42L * 12L,
  setequal(unique(index$year), 1985:2026),
  all(is.na(index$sampling_time)),
  sum(is.na(index$value)) == 28L,
  all(index$basis == "index_scale"),
  all(index$unit == "survey index (unit unresolved)"),
  nrow(weight) == 47L * 12L,
  all(weight$basis == "g_per_fish"),
  all(weight$unit == "g per fish"),
  nrow(maturity) == 47L * 12L,
  all(maturity$basis == "proportion"),
  nrow(mortality) == 12L,
  all(is.na(mortality$year)),
  all(is.na(mortality$year_basis)),
  all(mortality$value == 0.2),
  nrow(n_at_age) == 47L * 12L,
  setequal(unique(n_at_age$year), 1980:2026),
  nrow(f_at_age) == 46L * 12L,
  setequal(unique(f_at_age$year), 1980:2025),
  all(n_at_age$age_group[n_at_age$age == 14L] == "14+"),
  all(f_at_age$age_group[f_at_age$age == 14L] == "14+"),
  nrow(ssb) == 47L,
  nrow(reference_biomass) == 47L,
  nrow(recruitment) == 47L,
  nrow(fbar) == 46L,
  max(fbar$year) == 2025L,
  !any(fbar$year == 2026L),
  abs(ssb$value[ssb$year == 2026] - 150803) < 1,
  abs(reference_biomass$value[reference_biomass$year == 2026] - 249434) < 1,
  abs(recruitment$value[recruitment$year == 2026] - 29920) < 1,
  all(ssb$unit == "tonnes"),
  all(recruitment$unit == "thousand fish"),
  all(ssb$lwr < ssb$value & ssb$value < ssb$upr),
  abs(fbar$value[fbar$year == 2025] - 0.167390866) < 1e-8,
  !any(outputs$measure == "Fbar" & outputs$year == 2026L)
)

review <- paste(readLines(file.path(root, "source_reviews", "ices_saithe_iceland_2026.md"), warn = FALSE), collapse = " ")
stopifnot(
  grepl("28 blank age-year cells", review, fixed = TRUE),
  grepl("Table 2 `F` column's age range and calculation are not defined", review, fixed = TRUE),
  grepl("Through 2025, `HR` agrees to rounding", review, fixed = TRUE),
  grepl("exact 2026 production settings are not fully published", review, fixed = TRUE),
  grepl("0.20 as a month-midpoint approximation", review, fixed = TRUE)
)

cat("Icelandic saithe 2026 source records passed.\n")