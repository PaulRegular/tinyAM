source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_haddock_iceland_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
original <- source_data
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat

tinyAM::check_tam_data(obs)
stopifnot(
  identical(translated$years, 1979:2026),
  identical(translated$ages, 1:12),
  identical(translated$age_plus_group, 12L),
  nrow(obs$catch) == 48L * 12L,
  sum(is.na(obs$catch$obs)) == 12L,
  max(obs$catch$year[is.finite(obs$catch$obs)]) == 2025L,
  setequal(unique(obs$index$survey), c("IS-SMB", "IS-SMH")),
  max(obs$index$year[obs$index$survey == "IS-SMB"]) == 2026L,
  max(obs$index$year[obs$index$survey == "IS-SMH"]) == 2025L,
  all(obs$index$samp_time[obs$index$survey == "IS-SMB"] == 0.20),
  all(obs$index$samp_time[obs$index$survey == "IS-SMH"] == 0.80),
  nrow(obs$weight) == 48L * 12L,
  nrow(obs$maturity) == 48L * 12L,
  nrow(obs$pop_covariates) == 48L * 12L,
  all(obs$pop_covariates$M_assumption == 0.2),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "free"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$M_settings$process, "off"),
  identical(paste(deparse(dat$index_settings$q_form), collapse = ""), "~0 + q_key"),
  is.null(translated$fit),
  !"start_par" %in% names(translated$fit_call),
  identical(source_data$inputs,
    database$inputs[database$inputs$assessment_id == assessment_id, ])
)

cat("Icelandic haddock 2026 translation preparation passed.\n")
