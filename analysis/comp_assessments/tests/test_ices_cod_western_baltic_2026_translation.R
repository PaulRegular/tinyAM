pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_data <- read_assessment("ices_cod_western_baltic_2026", database)
translation <- new.env(parent = environment())
translation$source <- source_data
translation$do_fit <- FALSE
translation$silent <- TRUE
sys.source(
  "analysis/comp_assessments/scripts/translation/stocks/ices_cod_western_baltic_2026.R",
  envir = translation
)

obs <- translation$obs
years <- translation$years
ages <- translation$ages

stopifnot(
  identical(years, 1985:2023),
  identical(ages, 1:7),
  nrow(obs$catch) == 273L,
  nrow(obs$index) == 138L,
  setequal(unique(obs$index$survey), c("BITS Q1", "BITS Q4")),
  identical(as.integer(sort(unique(obs$index$year))), 2001:2023),
  identical(as.integer(sort(unique(obs$index$age))), 1:3),
  all(obs$index$samp_time[obs$index$survey == "BITS Q1"] == 0.1653494),
  all(obs$index$samp_time[obs$index$survey == "BITS Q4"] == 0.8477816),
  length(unique(obs$index$q_key)) == 6L,
  all(obs$maturity$spawn_time == 0),
  nrow(obs$pop_covariates) == length(years) * length(ages)
)

prepared <- tinyAM::prepare_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "iid", mu_form = ~0 + f_age_group, mean_ages = 3:5),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + sd_block),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
)

stopifnot(
  identical(prepared$years, years),
  identical(prepared$ages, ages),
  prepared$N_settings$init == "exp",
  prepared$F_settings$process == "iid",
  prepared$M_settings$process == "off"
)

cat("Western Baltic cod 2026 restricted observations and model preparation passed.\n")
