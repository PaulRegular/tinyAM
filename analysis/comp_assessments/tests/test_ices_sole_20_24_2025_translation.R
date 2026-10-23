source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

database <- read_database()
source_data <- read_assessment("ices_sole_20_24_2025", database)
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat
survey_names <- c("Q4Solesurvey1-9+", "Gillnet", "Trawl")

tinyAM::check_obs(obs)

stopifnot(
  identical(translated$years, 1984:2024),
  identical(translated$ages, 1:9),
  identical(translated$age_plus_group, 9L),
  nrow(obs$catch) == 369L,
  all(obs$catch$year %in% 1984:2024),
  setequal(unique(obs$catch$age), 1:9),
  sum(is.finite(obs$catch$obs)) == 328L,
  sum(is.na(obs$catch$obs)) == 41L,
  nrow(obs$index) == 379L,
  !anyNA(obs$index$obs),
  sum(obs$index$survey == survey_names[[1L]]) == 171L,
  setequal(unique(obs$index$survey), survey_names),
  all(obs$index$samp_time[obs$index$survey == survey_names[[1L]]] == 0.8754645),
  all(obs$index$samp_time[obs$index$survey == survey_names[[2L]]] == 0.56),
  all(obs$index$samp_time[obs$index$survey == survey_names[[3L]]] == 0.875),
  length(unique(obs$index$q_key)) == 12L,
  setequal(unique(obs$index$q_block[obs$index$survey == survey_names[[1L]]]), c("1", "2", "3", "4+")),
  all(as.character(obs$index$q_block[obs$index$survey != survey_names[[1L]] & obs$index$age >= 5]) == "5+"),
  nrow(obs$weight) == 369L,
  nrow(obs$maturity) == 369L,
  all(obs$maturity$spawn_time == 0),
  all(obs$pop_covariates$M_assumption == 0.1),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$rec$type, "rw"),
  identical(dat$F_settings$process, "ar1"),
  identical(deparse(dat$F_settings$mu_form), "~1"),
  identical(dat$F_settings$mean_ages, 4:8),
  identical(dat$M_settings$process, "off"),
  is.null(translated$start_par),
  is.null(translated$warm_call),
  !"start_par" %in% names(translated$fit_call),
  identical(source_data$inputs, database$inputs[database$inputs$assessment_id == "ices_sole_20_24_2025", ])
)

cat("ICES sole 20–24 translation structure passed.\n")
