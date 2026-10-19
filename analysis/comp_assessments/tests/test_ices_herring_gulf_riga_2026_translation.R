pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_record <- read_assessment("ices_herring_gulf_riga_2026", database)
stock <- .assessment_stock(source_record, do_fit = FALSE, silent = TRUE)
stopifnot(is.null(stock$error))
obs <- stock$env$obs

stopifnot(
  nrow(obs$catch) == 49L * 9L,
  nrow(obs$index) == 27L * 8L,
  nrow(obs$maturity) == 49L * 9L,
  all(obs$index$samp_time == 0.675),
  all(obs$catch$obs_sd_group[obs$catch$age == 0] == "age0"),
  all(obs$catch$obs_sd_group[obs$catch$age > 0] == "age1plus"),
  all(obs$index$q_key[obs$index$age == 1] == "age1"),
  all(obs$index$q_key[obs$index$age == 2] == "age2"),
  all(obs$index$q_key[obs$index$age >= 3] == "age3plus"),
  all(obs$index$obs_sd_group[obs$index$age == 1] == "age1"),
  all(obs$index$obs_sd_group[obs$index$age %in% 2:7] == "age2to7"),
  all(obs$index$obs_sd_group[obs$index$age == 8] == "age8plus"),
  identical(sort(unique(obs$maturity$spawn_time)), c(0.2, 0.25)),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$year %in% 1979:1983] == 0.25),
  all(obs$pop_covariates$M_assumption[!obs$pop_covariates$year %in% 1979:1983] == 0.20)
)

cat("Gulf of Riga tinyAM observation translation passed.\n")
