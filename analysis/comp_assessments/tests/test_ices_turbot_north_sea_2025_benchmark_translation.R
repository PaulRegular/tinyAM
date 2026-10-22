source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_turbot_north_sea_2025_benchmark"
source_data <- read_assessment(assessment_id, read_database())
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat
inputs <- source_data$inputs

tinyAM::check_obs(obs)

catch <- inputs[
  inputs$type == "catch" & inputs$measure == "numbers_at_age",
  , drop = FALSE
]
index <- inputs[
  inputs$type == "index" & inputs$measure == "numbers_at_age",
  , drop = FALSE
]
catch_expected <- match(paste(obs$catch$year, obs$catch$age),
  paste(catch$year, catch$age))
index_expected <- match(
  paste(obs$index$survey, obs$index$year, obs$index$age),
  paste(index$survey, index$year, index$age)
)
m_input <- inputs[inputs$type == "M" & inputs$measure == "natural_mortality_at_age", ]
m_expected <- m_input$value[match(obs$pop_covariates$age, m_input$age)]
weight_input <- inputs[inputs$type == "weight" & inputs$measure == "weight_at_age", ]
weight_output <- source_data$outputs[
  source_data$outputs$type == "biology" & source_data$outputs$measure == "weight_at_age", ]
key <- function(x) paste(x$year, x$age, sep = "\034")
weight_input_index <- match(key(weight_input), key(obs$weight))
weight_supplement_index <- match(key(translated$weight_supplement), key(obs$weight))
weight_output_index <- match(key(translated$weight_supplement), key(weight_output))
observed_catch <- is.finite(obs$catch$obs)

stopifnot(
  identical(translated$years, 1981:2023),
  identical(translated$ages, 1:8),
  identical(translated$age_plus_group, 8L),
  nrow(obs$catch) == 43L * 8L,
  nrow(obs$index) == nrow(index),
  nrow(obs$weight) == 43L * 8L,
  nrow(obs$maturity) == 43L * 8L,
  nrow(translated$weight_supplement) == 102L,
  all(is.finite(obs$weight$obs) & obs$weight$obs > 0),
  all(!is.na(weight_input_index)), all(!is.na(weight_supplement_index)),
  all(!is.na(weight_output_index)),
  isTRUE(all.equal(obs$weight$obs[weight_input_index], weight_input$value)),
  isTRUE(all.equal(obs$weight$obs[weight_supplement_index], weight_output$value[weight_output_index])),
  all(!is.na(catch_expected[observed_catch])),
  all(is.na(catch_expected[!observed_catch])),
  all(is.na(obs$catch$obs[!observed_catch])),
  all(!is.na(index_expected)),
  isTRUE(all.equal(obs$catch$obs[observed_catch], catch$value[catch_expected[observed_catch]])),
  isTRUE(all.equal(obs$index$obs, index$value[index_expected])),
  setequal(unique(obs$index$survey), c("BTS+COAST", "BSAS")),
  all(obs$index$samp_time == 0.705),
  all(obs$maturity$spawn_time == 0),
  length(unique(obs$index$q_key)) == 13L,
  length(unique(obs$catch$sd_group)) == 3L,
  length(unique(obs$index$sd_group)) == 5L,
  length(unique(obs$pop_covariates$F_sd_group)) == 4L,
  all(obs$index$q_key[obs$index$survey == "BTS+COAST" & obs$index$age %in% 6:7] ==
    obs$index$q_key[obs$index$survey == "BTS+COAST" & obs$index$age == 6][1]),
  all(obs$index$q_key[obs$index$survey == "BSAS" & obs$index$age %in% 7:8] ==
    obs$index$q_key[obs$index$survey == "BSAS" & obs$index$age == 7][1]),
  isTRUE(all.equal(obs$pop_covariates$M_assumption, m_expected)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$rec$type, "rw"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$F_settings$mean_ages, 2:6),
  identical(dat$M_settings$process, "off"),
  is.null(translated$start_par),
  identical(source_data$inputs, inputs),
  identical(translated$comparison_aggregates, c("abundance", "biomass", "ssb")),
  identical(translated$comparison_definitions$F_bar$status, "matched"),
  identical(translated$comparison_definitions$biomass$status, "matched"),
  grepl("official N", translated$comparison_definitions$biomass$reason, fixed = TRUE),
  identical(translated$comparison_definitions$ssb$status, "matched"),
  grepl("not native SAM SSB", translated$comparison_definitions$ssb$reason, fixed = TRUE),
  any(grepl("tinyAM::fit_tam", deparse(translated$fit_call), fixed = TRUE))
)

cat("North Sea turbot translation structure passed.\n")

