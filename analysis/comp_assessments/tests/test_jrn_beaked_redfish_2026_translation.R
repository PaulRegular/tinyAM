pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

assessment_id <- "jrn_redfish_beaked_2026"
database <- read_database()
record <- read_assessment(assessment_id, database)
original_inputs <- record$inputs
original_outputs <- record$outputs
stock <- .assessment_stock(record, do_fit = FALSE, silent = TRUE)
prepared <- stock$env
obs <- prepared$obs
years <- 1992:2025
ages <- 6:19

stopifnot(
  is.null(stock$error),
  is.null(prepared$fit),
  identical(prepared$fit_stage, "translation"),
  identical(prepared$years, years),
  identical(prepared$ages, ages),
  identical(deparse(prepared$index_settings$q_form),
    deparse(~ 0 + survey + survey:age)),
  identical(deparse(prepared$index_settings$sd_form),
    deparse(~ 0 + survey)),
  isTRUE(tinyAM::check_obs(obs)),
  nrow(obs$catch) == length(years) * length(ages),
  identical(as.integer(sort(unique(obs$catch$year))), years),
  identical(as.integer(sort(unique(obs$catch$age))), ages),
  all(obs$catch$unit == "thousand fish"),
  setequal(unique(obs$index$survey), c("Ecosystem survey", "Winter bottom trawl")),
  identical(as.integer(sort(unique(obs$index$age))), 6:15),
  all(obs$index$unit == "thousand fish"),
  all(obs$index$samp_time[obs$index$survey == "Ecosystem survey"] == 0.75),
  all(obs$index$samp_time[obs$index$survey == "Winter bottom trawl"] == 0.12),
  all(obs$pop_covariates$M_assumption == 0.05),
  nrow(obs$weight) == length(years) * length(ages),
  nrow(obs$maturity) == length(years) * length(ages),
  identical(record$inputs, original_inputs),
  identical(record$outputs, original_outputs)
)

comparison <- prepared$comparison_outputs
derived_ssb <- comparison[
  comparison$measure == "mature_biomass_at_age", , drop = FALSE
]
total_f <- comparison[
  comparison$measure == "fishing_mortality_at_age", , drop = FALSE
]
source_n <- original_outputs[
  original_outputs$measure == "numbers_at_age" &
    original_outputs$year %in% years & original_outputs$age %in% ages,
  , drop = FALSE
]
stock_weight <- original_inputs[
  original_inputs$type == "weight" & original_inputs$measure == "weight_at_age",
  , drop = FALSE
]
maturity <- original_inputs[
  original_inputs$type == "maturity" & original_inputs$measure == "maturity_at_age",
  , drop = FALSE
]
expected_ssb <- source_n$value *
  stock_weight$value[match(source_n$age, stock_weight$age)] *
  maturity$value[match(source_n$age, maturity$age)]

stopifnot(
  nrow(derived_ssb) == length(years) * length(ages),
  isTRUE(all.equal(derived_ssb$value, expected_ssb, tolerance = 1e-12)),
  all(derived_ssb$unit == "tonnes"),
  all(derived_ssb$source_type == "derived_source_output"),
  all(is.na(derived_ssb$se)),
  all(is.na(total_f$fleet)),
  prepared$comparison_definitions$N$status == "approximate",
  prepared$comparison_definitions$F$status == "approximate",
  prepared$comparison_definitions$ssb$status == "approximate",
  prepared$comparison_definitions$recruitment$status == "non_equivalent",
  prepared$comparison_definitions$abundance$status == "non_equivalent",
  prepared$comparison_definitions$M$status == "matched"
)

cat("Beaked redfish age-6+ translation and comparison outputs passed.\n")
