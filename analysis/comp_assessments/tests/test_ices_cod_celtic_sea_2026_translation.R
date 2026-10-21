source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/read_database.R")

source_data <- read_assessment("ices_cod_celtic_sea_2026")
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat
source_start <- translated$fit_start()
source_N <- .translation_start_surface(
  source_data$outputs, "population", "numbers_at_age", 1980:2025, 1:7
)
source_F <- .translation_start_surface(
  source_data$outputs, "mortality", "fishing_mortality_at_age", 1980:2025, 1:7
)

stopifnot(
  tinyAM::check_obs(obs),
  identical(translated$years, 1980:2025),
  identical(translated$ages, 1:7),
  translated$age_plus_group == 7L,
  nrow(obs$catch) == 46L * 7L,
  sum(is.na(obs$catch$obs)) == 48L,
  nrow(obs$weight) == 46L * 7L,
  nrow(obs$maturity) == 46L * 7L,
  nrow(obs$index) == 68L,
  setequal(unique(obs$index$survey), "VAST_IGFS_EVHOE"),
  all(obs$index$samp_time == 0.855),
  all(obs$index$age %in% 1:3),
  all(is.finite(obs$index$obs) & obs$index$obs > 0),
  all(obs$index$q_key[obs$index$age == 1] == "age1"),
  all(obs$index$q_key[obs$index$age == 2] == "age2"),
  all(obs$index$q_key[obs$index$age == 3] == "age3"),
  all(obs$index$obs_sd_group == "VAST"),
  all(obs$catch$obs_sd_group[obs$catch$age == 1] == "age1"),
  all(obs$catch$obs_sd_group[obs$catch$age == 2] == "age2"),
  all(obs$catch$obs_sd_group[obs$catch$age >= 3] == "age3plus"),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 1] == 0.501),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 2] == 0.330),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 3] == 0.264),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age >= 5] == 0.211),
  identical(sort(unique(obs$maturity$spawn_time)), 0),
  dat$N_settings$init == "exp",
  dat$N_settings$process == "iid",
  dat$rec_settings$process == "rw",
  dat$F_settings$process == "cor_rw",
  identical(dat$F_settings$mean_ages, 2:5),
  dat$M_settings$process == "off",
  setequal(unique(obs$index$q_key), c("age1", "age2", "age3")),
  identical(dimnames(source_start$log_f), list(as.character(1980:2025), as.character(1:7))),
  identical(dimnames(source_start$log_n), list(as.character(1981:2025), as.character(2:7))),
  identical(names(source_start$log_r), as.character(1981:2025)),
  isTRUE(all.equal(unname(source_start$log_r0),
    log(source_data$outputs$value[source_data$outputs$measure == "numbers_at_age" &
      source_data$outputs$year == 1980 & source_data$outputs$age == 1]))),
  isTRUE(all.equal(source_start$log_r, log(source_N[-1L, 1L]))),
  isTRUE(all.equal(source_start$log_n, log(source_N[-1L, -1L, drop = FALSE]))),
  isTRUE(all.equal(source_start$log_f, log(source_F))),
  all(is.finite(source_start$log_f)),
  all(is.finite(source_start$log_n)),
  all(is.finite(source_start$log_r)),
  is.null(translated$fit_call$start_par),
  grepl("aggregate commercial LPUE likelihood", paste(translated$background, collapse = " ")),
  grepl("stock-annex PDF was attempted but returned HTTP 401", paste(translated$background, collapse = " "))
)

cat("Celtic Sea cod tinyAM observation translation and settings passed.\n")
