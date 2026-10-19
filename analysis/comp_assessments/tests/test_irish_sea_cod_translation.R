source("analysis/comp_assessments/tests/helper_stock.R")
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_data <- read_assessment("ices_cod_irish_sea_2026", database)
translated <- .test_stock(source_data)
obs <- translated$obs
catch <- obs$catch
index <- obs$index

excluded <- catch$year %in% 2003:2005

stopifnot(
  tinyAM::check_obs(obs),
  identical(translated$years, 1968:2018),
  identical(translated$ages, 0:6),
  nrow(catch) == 51L * 7L,
  all(is.na(catch$obs[excluded])),
  all(is.finite(catch$obs[!excluded])),
  nrow(obs$weight) == 51L * 7L,
  nrow(obs$maturity) == 51L * 7L,
  all(is.finite(obs$weight$obs)),
  all(is.finite(obs$maturity$obs)),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 0] == 1.7),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 1] == 0.714),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 2] == 0.68),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age >= 3] == 0.65),
  nrow(index) == 185L,
  setequal(unique(index$survey), c(
    "NIGFS-WIBTS-Q1 (G7144)",
    "NIGFS-WIBTS-Q4 (G7655)",
    "UK FSP (B7897)"
  )),
  all(index$samp_time[index$survey != "NIGFS-WIBTS-Q4 (G7655)"] == 0.125),
  all(index$samp_time[index$survey == "NIGFS-WIBTS-Q4 (G7655)"] == 0.875),
  nlevels(index$q_key) == 10L,
  all(obs$maturity$spawn_time == 0),
  setequal(unique(obs$pop_covariates$recruitment_block), c("pre_2002", "2002_plus")),
  identical(translated$dat$N_settings$process, "iid"),
  identical(translated$dat$N_settings$init, "exp"),
  identical(translated$dat$rec_settings$process, "iid"),
  identical(translated$dat$rec$curve$type, "bh"),
  identical(translated$dat$rec$curve$lag, 1L),
  identical(translated$dat$F_settings$process, "ar1"),
  identical(translated$dat$F_settings$mean_ages, 2:4),
  identical(translated$dat$M_settings$process, "off"),
  identical(translated$comparison_aggregates, "ssb"),
  identical(translated$comparison_definitions$ssb$status, "approximate"),
  identical(translated$comparison_definitions$recruitment$status, "approximate"),
  !"start_par" %in% names(translated$fit_call)
)

cat("Irish Sea cod 1968-2018 source-aligned translation structure passed.\n")
