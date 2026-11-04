source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_saithe_iceland_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
original_inputs <- source_data$inputs
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat

stopifnot(
  tinyAM::check_tam_data(obs),
  identical(translated$years, 1980:2026),
  identical(translated$ages, 3:14),
  identical(translated$age_plus_group, 14L),
  nrow(obs$catch) == 47L * 12L,
  sum(is.na(obs$catch$obs)) == 12L,
  max(obs$catch$year[is.finite(obs$catch$obs)]) == 2025L,
  nrow(obs$index) == 476L,
  sum(is.na(obs$index$obs)) == 0L,
  setequal(unique(obs$index$survey), "IS-SMB"),
  identical(range(as.integer(obs$index$year)), c(1985L, 2026L)),
  all(obs$index$samp_time == 0.20),
  nrow(obs$weight) == 47L * 12L,
  nrow(obs$maturity) == 47L * 12L,
  nrow(obs$pop_covariates) == 47L * 12L,
  all(obs$pop_covariates$M_assumption == 0.2),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "free"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$F_settings$mean_ages, 4:9),
  identical(dat$M_settings$process, "off"),
  identical(paste(deparse(dat$index_settings$q_form), collapse = ""), "~0 + q_block"),
  is.null(translated$fit),
  is.function(translated$fit_model),
  is.function(translated$fit_start),
  "start_par" %in% names(translated$fit_call),
  identical(as.character(translated$fit_call$start_par), "start_par"),
  identical(original_inputs, source_data$inputs)
)

source_start <- translated$fit_start()
source_N <- .translation_start_surface(
  source_data$outputs, "population", "numbers_at_age", translated$years, translated$ages
)
source_F <- .translation_start_surface(
  source_data$outputs, "mortality", "fishing_mortality_at_age", 1980:2025, translated$ages
)
positive_F <- is.finite(source_F) & source_F > 0
stopifnot(
  all(is.finite(source_start$log_n0)),
  all(is.finite(source_start$log_r)),
  all(is.finite(source_start$log_n)),
  all(is.finite(source_start$log_f)),
  isTRUE(all.equal(unname(source_start$log_r0), log(source_N[1L, 1L]))),
  isTRUE(all.equal(unname(source_start$log_n0), unname(log(source_N[1L, -1L])))),
  isTRUE(all.equal(unname(source_start$log_r), unname(log(source_N[-1L, 1L])))),
  isTRUE(all.equal(unname(source_start$log_n), unname(log(source_N[-1L, -1L, drop = FALSE])))),
  isTRUE(all.equal(
    unname(source_start$log_f[as.character(1980:2025), , drop = FALSE][positive_F]),
    unname(log(source_F)[positive_F])
  ))
)

cat("Icelandic saithe 2026 translation preparation and source starts passed.\n")