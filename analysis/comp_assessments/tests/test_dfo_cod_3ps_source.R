root <- file.path("analysis", "comp_assessments")
read_table <- function(file) {
  read.csv(file.path(root, "database", file), na.strings = c("", "NA"),
    stringsAsFactors = FALSE, check.names = FALSE)
}

stocks <- read_table("stocks.csv")
assessments <- read_table("assessments.csv")
inputs <- read_table("inputs.csv")
outputs <- read_table("outputs.csv")
stock <- stocks[stocks$stock_id == "ck_dfo_morhua_3ps", ]
a <- assessments[assessments$assessment_id == "dfo_cod_3ps_2024", ]
i <- inputs[inputs$assessment_id == "dfo_cod_3ps_2024", ]
o <- outputs[outputs$assessment_id == "dfo_cod_3ps_2024", ]

stopifnot(
  nrow(stock) == 1L,
  stock$authority_stock_id == "3Ps",
  stock$common_name == "Atlantic cod",
  nrow(a) == 1L,
  a$assessment_year == 2024,
  a$terminal_year == 2024,
  a$estimate_terminal_year == 2025,
  a$model_family == "HYBRID",
  a$is_current,
  a$is_applied,
  a$assumptions_status == "partial",
  a$inputs_status == "partial",
  a$outputs_status == "partial",
  nrow(i) == 4911L,
  sum(i$type == "catch" & i$measure == "numbers_at_age") == 845L,
  all(i$unit[i$type == "catch" & i$measure == "numbers_at_age"] == "thousand fish"),
  all(grepl("unit convention", i$transformation[i$type == "catch" & i$measure == "numbers_at_age"], fixed = TRUE)),
  all(i$unit[i$type == "index"] == "survey index (unit unresolved)"),
  sum(i$type == "catch" & i$measure == "total_biomass") == 66L,
  sum(i$type == "catch" & i$measure == "landings_lower_bound") == 66L,
  sum(i$type == "catch" & i$measure == "landings_upper_bound") == 66L,
  sum(i$type == "M" & i$measure == "condition_index_for_mortality") == 66L,
  sum(i$type == "index") == 1228L,
  sum(i$type == "maturity") == 858L,
  sum(i$type == "weight") == 858L,
  sum(i$type == "catch_weight") == 858L,
  all(is.finite(i$value)),
  all(i$source_type == "official_table"),
  all(grepl("DFO Research Document 2026/051", i$source_reference, fixed = TRUE)),
  max(i$year) == 2024,
  !any(i$year == 2025),
  !any(i$type == "index" & i$survey == "DFO RV Spring" & i$year %in% c(2006, 2020, 2023)),
  i$value[i$type == "catch" & i$measure == "numbers_at_age" & i$year == 1959 & i$age == 3] == 1001,
  i$value[i$type == "catch" & i$measure == "landings_lower_bound" & i$year == 1959] == 54153,
  i$value[i$type == "maturity" & i$year == 1959 & i$age == 2] == 0.001,
  i$value[i$type == "weight" & i$year == 1959 & i$age == 14] == 15.02,
  i$value[i$type == "catch_weight" & i$year == 1959 & i$age == 14] == 9.16,
  i$value[i$type == "index" & i$survey == "DFO RV Spring" & i$year == 2024 & i$age == 14] == 0,
  i$value[i$type == "M" & i$measure == "condition_index_for_mortality" & i$year == 1959] == -0.927,
  all(is.na(i$sampling_time[i$type == "index"])),
  nrow(o) == 5L,
  o$value[o$measure == "SSB" & o$year == 2024] == 29.9,
  o$lwr[o$measure == "SSB" & o$year == 2024] == 24,
  o$upr[o$measure == "SSB" & o$year == 2024] == 37.3,
  o$value[o$measure == "SSB" & o$year == 2025] == 33.2,
  o$value[o$measure == "recruitment" & o$year == 2024] == 36.6,
  o$age[o$measure == "recruitment" & o$year == 2024] == 2,
  o$value[o$measure == "Fbar" & o$year == 2024] == 0.04,
  o$age_group[o$measure == "Fbar"] == "5-8",
  o$value[o$measure == "Mbar" & o$year == 2024] == 0.32,
  !any(o$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age"))
)

cat("3Ps 2024 HYBRID source coverage, units, gaps, and reported summaries passed.\n")

