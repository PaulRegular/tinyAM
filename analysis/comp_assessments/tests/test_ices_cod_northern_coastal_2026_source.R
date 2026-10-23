root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

id <- "ices_cod_northern_coastal_2026"
database <- read_database()
assessment <- read_assessment(id, database)
inputs <- assessment$inputs
outputs <- assessment$outputs
legacy <- database$stocks[database$stocks$stock_id == "ck_ices_afwg_morhua_coastnor_1_2", , drop = FALSE]

input_rows <- function(type, measure = NULL) {
  x <- inputs[inputs$type == type, , drop = FALSE]
  if (!is.null(measure)) x <- x[x$measure == measure, , drop = FALSE]
  x
}
output_rows <- function(measure) outputs[outputs$measure == measure, , drop = FALSE]

stopifnot(
  nrow(assessment$stocks) == 1L,
  assessment$stocks$stock_id == "ices_cod_northern_coastal",
  assessment$stocks$authority_stock_id == "cod.27.1-2.coastN",
  is.na(assessment$stocks$charbonneau_id),
  nrow(legacy) == 1L,
  !is.na(legacy$charbonneau_id),
  assessment$assessment$assessment_year == 2026,
  assessment$assessment$terminal_year == 2025,
  assessment$assessment$inputs_status == "complete",
  assessment$assessment$outputs_status == "complete",
  assessment$assessment$assumptions_status == "partial",
  nrow(inputs) == 1700L,
  nrow(input_rows("catch")) == 32L * 9L,
  all(input_rows("catch")$unit == "thousand fish"),
  nrow(input_rows("catch_weight")) == 32L * 9L,
  all(input_rows("catch_weight")$unit == "kg per fish"),
  nrow(input_rows("weight")) == 32L * 9L,
  nrow(input_rows("M")) == 32L * 9L,
  all(input_rows("M")$unit == "per year"),
  nrow(input_rows("maturity")) == 31L * 10L,
  all(input_rows("maturity")$sampling_time == 0.8),
  nrow(input_rows("index", "index_by_age_group")) == 31L,
  input_rows("index", "index_by_age_group")$value[
    input_rows("index", "index_by_age_group")$year == 2000
  ] == 35255,
  is.na(input_rows("index", "index_by_age_group")$value[
    input_rows("index", "index_by_age_group")$year == 2021
  ]),
  nrow(input_rows("index", "numbers_at_age")) == 23L * 9L,
  all(input_rows("index", "numbers_at_age")$unit == "thousand fish"),
  all(input_rows("index", "numbers_at_age")$sampling_time == 0.8),
  setequal(unique(as.integer(input_rows("index", "numbers_at_age")$age)), 2:10),
  nrow(output_rows("numbers_at_age")) == 32L * 9L,
  nrow(output_rows("fishing_mortality_at_age")) == 32L * 9L,
  all(is.na(output_rows("numbers_at_age")$se)),
  all(is.na(output_rows("fishing_mortality_at_age")$se)),
  nrow(output_rows("recruitment")) == 32L,
  all(output_rows("recruitment")$age == 3),
  nrow(output_rows("SSB")) == 32L,
  nrow(output_rows("Fbar")) == 32L,
  nrow(output_rows("total_biomass")) == 32L,
  nrow(output_rows("q")) == 7L,
  isTRUE(all.equal(
    output_rows("q")$value,
    c(0.219, 0.141, 0.065, 0.109, 0.150, 0.166, 0.201)
  ))
)

prepared <- run_assessment(id, database = database, fit = FALSE, silent = TRUE)
obs <- prepared$obs
comparison_biomass <- prepared$comparison$comparison_outputs[
  prepared$comparison$comparison_outputs$measure == "biomass_at_age",
  , drop = FALSE
]
source_n <- outputs[
  outputs$measure == "numbers_at_age" & outputs$year %in% 1995:2025 &
    outputs$age %in% 3:10,
  , drop = FALSE
]
source_weight <- inputs[
  inputs$type == "weight" & inputs$measure == "weight_at_age" &
    inputs$year %in% 1995:2025 & inputs$age %in% 3:10,
  , drop = FALSE
]
expected_biomass <- source_n$value[match(
  paste(comparison_biomass$year, comparison_biomass$age),
  paste(source_n$year, source_n$age)
)] * source_weight$value[match(
  paste(comparison_biomass$year, comparison_biomass$age),
  paste(source_weight$year, source_weight$age)
)]
tam_data <- tinyAM::prepare_tam(
  data = obs,
  years = 1995:2025,
  ages = 3:10,
  N_settings = list(process = "rw", init = "exp"),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 4:8),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~1)
)

stopifnot(
  prepared$diagnostics$status == "not_fitted",
  nrow(comparison_biomass) == 31L * 8L,
  all(comparison_biomass$unit == "tonnes"),
  all(comparison_biomass$source_type == "derived_source_output"),
  all(is.na(comparison_biomass$se)),
  isTRUE(all.equal(comparison_biomass$value, expected_biomass)),
  nrow(obs$catch) == 31L * 8L,
  nrow(obs$index) == 23L * 8L,
  all(obs$index$survey == "Coastal swept-area index"),
  all(obs$index$samp_time == 0.8),
  identical(levels(obs$index$q_key), c("age3", "age4", "ages5_9", "age10_plus")),
  all(obs$maturity$spawn_time == 0.8),
  is.null(tam_data$plus_ages),
  tam_data$N_settings$process == "rw",
  tam_data$rec_settings$process == "rw",
  tam_data$F_settings$process == "cor_rw",
  all(is.finite(obs$pop_covariates$M_assumption))
)

cat("Northern Norwegian coastal cod source translation checks passed.
")
