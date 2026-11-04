source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")

assessment_id <- "ices_sole_celtic_7fg_2026"
database <- read_database()
assessment <- read_assessment(assessment_id, database)
inputs <- assessment$inputs
outputs <- assessment$outputs
belgian_recent <- inputs[
  inputs$type == "index" & inputs$survey == "BE-CBT3-2006-2025", , drop = FALSE
]

stopifnot(
  identical(assessment$stocks$stock_id, "ck_ices_wgsce_solea_cs7f_7g"),
  identical(assessment$stocks$authority_stock_id, "sol.27.7fg"),
  identical(assessment$stocks$ocean, "Northeast Atlantic"),
  as.integer(assessment$assessment$terminal_year) == 2025L,
  identical(assessment$assessment$model_family, "SAM"),
  identical(assessment$assessment$inputs_status, "complete"),
  identical(assessment$assessment$outputs_status, "complete"),
  sum(inputs$type == "catch" & inputs$measure == "numbers_at_age") == 55L * 10L,
  sum(inputs$type == "catch_weight" & inputs$measure == "weight_at_age") == 55L * 10L,
  sum(inputs$type == "weight" & inputs$measure == "weight_at_age") == 55L * 10L,
  sum(inputs$type == "index" & inputs$measure == "numbers_at_age") == 38L * 10L,
  sum(inputs$type == "index" & inputs$measure == "relative_biomass_index") == 88L,
  nrow(belgian_recent) == 20L,
  all(belgian_recent$value >= 0),
  isTRUE(all.equal(belgian_recent$value[belgian_recent$year == 2013], 1.9697236004292)),
  isTRUE(all.equal(belgian_recent$value[belgian_recent$year == 2019], 3.5890134575644)),
  isTRUE(all.equal(belgian_recent$value[belgian_recent$year == 2025], 3.248300029998741)),
  sum(inputs$type == "maturity" & inputs$measure == "maturity_at_age") == 15L,
  sum(inputs$type == "M" & inputs$measure == "natural_mortality_at_age") == 10L,
  all(inputs$value[inputs$type == "M"] == 0.1),
  all(inputs$source_type %in% c("official_table", "official_document")),
  setequal(as.integer(unique(inputs$age[inputs$type == "index" & inputs$measure == "numbers_at_age"])), 0:9),
  setequal(as.integer(unique(inputs$year[inputs$type == "index" & inputs$measure == "numbers_at_age"])), 1988:2025),
  !any(inputs$value[inputs$type == "index" & inputs$measure == "numbers_at_age" & inputs$age %in% 1:5] == 0),
  all(grepl("Table 22\\.", inputs$source_reference[inputs$assessment_id == assessment_id])),
  sum(outputs$type == "population" & outputs$measure == "numbers_at_age") == 55L * 10L,
  sum(outputs$type == "mortality" & outputs$measure == "fishing_mortality_at_age") == 55L * 10L,
  sum(outputs$measure == "SSB") == 55L,
  sum(outputs$measure == "recruitment") == 55L,
  identical(outputs$value[outputs$measure == "SSB" & outputs$year == "2025"], 4838),
  identical(outputs$value[outputs$measure == "recruitment" & outputs$year == "2025"], 3780),
  identical(outputs$value[outputs$measure == "Fbar" & outputs$year == "2025"], 0.291),
  all(grepl("interval level is not stated", outputs$notes[!is.na(outputs$lwr)]))
)

translated <- .test_stock(assessment)
stopifnot(
  identical(translated$years, 1971:2025),
  identical(translated$ages, 1:10),
  identical(translated$age_plus_group, 10L),
  nrow(translated$obs$catch) == 55L * 10L,
  sum(is.na(translated$obs$catch$obs)) == 66L,
  all(translated$obs$catch$unit == "thousand fish"),
  nrow(translated$obs$index) == 38L * 5L,
  setequal(unique(translated$obs$index$age), 1:5),
  all(translated$obs$index$samp_time == 0.625),
  all(translated$obs$pop_covariates$M_assumption == 0.1),
  identical(translated$settings$N_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$rec_settings$process, "rw"),
  identical(translated$settings$F_settings$process, "ar1"),
  identical(translated$settings$M_settings$process, "off"),
  is.null(translated$warm_call),
  !is.null(translated$fit_call)
)

cat("Celtic Sea sole source coverage, units, outputs, and translation preparation passed.\n")
