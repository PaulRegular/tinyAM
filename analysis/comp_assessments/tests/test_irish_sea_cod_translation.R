source("analysis/comp_assessments/tests/helper_stock.R")
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_data <- read_assessment("ices_cod_irish_sea_2026", database)
translated <- .test_stock(source_data)
obs <- translated$obs
catch <- obs$catch
index <- obs$index

excluded <- catch$year %in% c(2003:2005, 2020)

stopifnot(
  tinyAM::check_obs(obs),
  identical(translated$years, 1968:2025),
  identical(translated$ages, 1:6),
  nrow(catch) == 58L * 6L,
  all(is.na(catch$obs[excluded])),
  all(is.finite(catch$obs[!excluded])),
  nrow(obs$weight) == 58L * 6L,
  nrow(obs$maturity) == 58L * 6L,
  all(is.finite(obs$weight$obs)),
  all(is.finite(obs$maturity$obs)),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 1] == 0.714),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 2] == 0.68),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age >= 3] == 0.65),
  nrow(index) == 224L,
  setequal(unique(index$survey), c("NIGFS-WIBTS-Q1 (G7144)", "UK FSP (B7897)")),
  all(index$samp_time == 0.125),
  nlevels(index$q_key) == 9L,
  all(obs$maturity$spawn_time == 0),
  identical(translated$dat$N_settings$process, "iid"),
  identical(translated$dat$N_settings$init, "exp"),
  identical(translated$dat$rec_settings$process, "rw"),
  identical(translated$dat$F_settings$process, "ar1"),
  identical(translated$dat$F_settings$mean_ages, 2:4),
  identical(translated$dat$M_settings$process, "off"),
  identical(translated$comparison_aggregates, "ssb"),
  identical(translated$comparison_definitions$ssb$status, "approximate"),
  !"start_par" %in% names(translated$fit_call)
)

cat("Irish Sea cod age-1+ translation structure passed.\n")
