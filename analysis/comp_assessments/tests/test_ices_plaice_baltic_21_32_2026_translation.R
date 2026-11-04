source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_plaice_baltic_21_32_2026"
database <- read_database()
source_record <- read_assessment(assessment_id, database)
translated <- .test_stock(source_record)
obs <- translated$obs
dat <- translated$dat
q1 <- "Q1IBTS+BITS1-7+"
q4 <- "Q4IBTS+BITS1-7+"

tinyAM::check_tam_data(obs)
m_input <- source_record$inputs[source_record$inputs$measure == "natural_mortality_at_age", ]
m_index <- match(paste(obs$pop_covariates$year, obs$pop_covariates$age), paste(m_input$year, m_input$age))
stopifnot(
  identical(translated$years, 2002:2025),
  identical(translated$ages, 1:7),
  identical(translated$age_plus_group, 7L),
  nrow(obs$catch) == 168L,
  sum(is.finite(obs$catch$obs)) == 161L,
  sum(is.na(obs$catch$obs)) == 7L,
  setequal(unique(obs$index$survey), c(q1, q4)),
  sum(obs$index$survey == q1) == 168L,
  sum(obs$index$survey == q4) == 161L,
  all(obs$index$samp_time[obs$index$survey == q1] == 0.1507111),
  all(obs$index$samp_time[obs$index$survey == q4] == 0.8173179),
  length(unique(obs$index$q_key[obs$index$survey == q1])) == 5L,
  setequal(unique(obs$index$q_block[obs$index$survey == q1]), c("1", "2", "3-5", "6", "7+")),
  setequal(unique(obs$index$q_block[obs$index$survey == q4]), c("1", "2", "3", "4+")),
  length(unique(obs$index$sd_block[obs$index$survey == q1])) == 1L,
  length(unique(obs$index$sd_block[obs$index$survey == q4])) == 3L,
  identical(translated$settings$N_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$rec_settings$process, "rw"),
  identical(translated$settings$F_settings$process, "ar1"),
  identical(translated$settings$F_settings$mean_ages, 3:5),
  identical(translated$settings$M_settings$process, "off"),
  is.null(translated$warm_call),
  !"start_par" %in% names(translated$fit_call),
  all(is.finite(obs$pop_covariates$M_assumption)),
  nrow(obs$pop_covariates) == length(translated$years) * length(translated$ages),
  !anyNA(m_index),
  max(abs(obs$pop_covariates$M_assumption - as.numeric(m_input$value[m_index]))) < 1e-12,
  identical(source_record$inputs,
    database$inputs[database$inputs$assessment_id == assessment_id, ])
)

cat("ICES Baltic plaice translation structure passed.\n")
