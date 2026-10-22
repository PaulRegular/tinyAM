source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")

assessment_id <- "ices_haddock_irish_sea_7a_2025"
advice_id <- "ices_haddock_irish_sea_7a_2026_summary"
database <- read_database()
assessment <- read_assessment(assessment_id, database)
i <- assessment$inputs
o <- assessment$outputs

stopifnot(
  identical(assessment$stocks$common_name, "Irish Sea haddock"),
  identical(assessment$stocks$ocean, "Northeast Atlantic"),
  as.integer(assessment$assessment$terminal_year) == 2024L,
  identical(assessment$assessment$inputs_status, "partial"),
  identical(assessment$assessment$outputs_status, "partial"),
  nrow(i[i$type == "catch" & i$measure == "numbers_at_age", ]) == 32L * 6L,
  nrow(i[i$type == "weight" & i$measure == "weight_at_age", ]) == 32L * 6L,
  nrow(i[i$type == "maturity" & i$measure == "maturity_at_age", ]) == 32L * 6L,
  identical(i$value[i$type == "catch" & i$year == 1993 & i$age == 0], 790),
  identical(i$value[i$type == "weight" & i$year == 1993 & i$age == 0], 0.02),
  identical(i$value[i$type == "maturity" & i$year == 1993 & i$age == 2], 0.762),
  identical(i$value[i$type == "index" & i$survey == "NIGFS-WIBTS-Q1" & i$year == 1995 & i$age == 3], 0.1),
  all(is.na(i$value[i$type == "index" & i$survey == "NIGFS-WIBTS-Q4" & i$year == 2024])),
  all(is.na(i$value[i$type == "index" & i$survey == "NI-MIK" & i$year %in% c(2020, 2022)])),
  all(is.na(i$value[i$type == "index" & i$survey == "UK-FSPW" & i$year == 2014])),
  identical(i$value[i$type == "M" & i$age == 0], 2.44),
  sum(i$type == "index" & i$survey == "NIGFS-WIBTS-Q1") == 32L * 5L,
  sum(i$type == "index" & i$survey == "NIGFS-WIBTS-Q4") == 31L * 5L,
  sum(i$type == "index" & i$survey == "NI-MIK") == 31L,
  sum(i$type == "index" & i$survey == "UK-FSPW") == 20L * 5L,
  identical(unique(i$unit[i$type == "index"]), "survey index (unit unresolved)"),
  nrow(o[o$measure == "SSB" & o$year >= 1993 & o$year <= 2024, ]) == 32L,
  identical(o$value[o$measure == "SSB" & o$year == 2024], 6308),
  identical(o$lwr[o$measure == "SSB" & o$year == 2024], 4217),
  identical(o$upr[o$measure == "SSB" & o$year == 2024], 8399),
  identical(o$value[o$measure == "Fbar" & o$year == 2024], 0.08),
  identical(o$age_group[o$measure == "Fbar"], rep("2-4", 32L)),
  !any(o$measure %in% c("numbers_at_age", "fishing_mortality_at_age")),
  isTRUE(database$assessments$is_current[database$assessments$assessment_id == advice_id]),
  identical(database$outputs$value[database$outputs$assessment_id == advice_id], c(586, 586))
)

analysis_inputs <- i[!(i$type == "index" & is.na(i$value)), , drop = FALSE]
obs <- database_to_tam_obs(
  assessment_id, analysis_inputs, years = 1993:2024, ages = 0:5,
  surveys = c("NIGFS-WIBTS-Q1", "NIGFS-WIBTS-Q4", "NI-MIK", "UK-FSPW"),
  assumptions = assessment$assumptions
)
stopifnot(
  identical(sort(unique(obs$index$survey)), c("NI-MIK", "NIGFS-WIBTS-Q1", "NIGFS-WIBTS-Q4", "UK-FSPW")),
  identical(as.integer(sort(unique(obs$index$year))), 1993:2024),
  all(is.finite(obs$pop_covariates$M_assumption)),
  nrow(obs$weight) == 32L * 6L,
  setequal(unique(obs$weight$age), 0:5),
  any(obs$index$obs == 0)
)

source("analysis/comp_assessments/R/run_assessment.R")
prepared <- run_assessment(assessment_id, fit = FALSE, dashboard = FALSE, cache = FALSE)
stopifnot(
  identical(prepared$diagnostics$status, "not_fitted"),
  !is.null(prepared$obs),
  nrow(prepared$obs$catch) == 32L * 6L,
  identical(prepared$comparison$comparison_definitions$F_bar$status, "non_equivalent"),
  identical(prepared$comparison$comparison_definitions$M$status, "approximate"),
  setequal(unique(prepared$obs$index$survey), c(
    "NI-MIK", "NIGFS-WIBTS-Q1", "NIGFS-WIBTS-Q4", "UK-FSPW"
  ))
)

cat("Irish Sea haddock source coverage, units, missingness, outputs and observation translation passed.\n")
