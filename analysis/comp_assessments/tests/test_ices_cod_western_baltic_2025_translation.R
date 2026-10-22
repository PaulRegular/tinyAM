pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_data <- read_assessment("ices_cod_western_baltic_2025", database)
translation <- new.env(parent = environment())
translation$source <- source_data
translation$do_fit <- FALSE
translation$silent <- TRUE
sys.source(
  "analysis/comp_assessments/scripts/translation/stocks/ices_cod_western_baltic_2025.R",
  envir = translation
)

obs <- translation$obs
years <- translation$years
ages <- translation$ages
catch <- obs$catch
index <- obs$index

stopifnot(
  identical(years, 1985:2025),
  identical(ages, 0:7),
  nrow(index) == 263L,
  setequal(unique(index$survey), c("BITS Q4", "BITS Q1", "FEJUCS")),
  all(index$samp_time[index$survey == "BITS Q4"] == 0.8477816),
  all(index$samp_time[index$survey == "BITS Q1"] == 0.1653494),
  all(index$samp_time[index$survey == "FEJUCS"] == 0.85),
  sum(catch$year == 2023L & catch$age == 7L & is.finite(catch$obs)) == 1L,
  all(catch$sd_block[catch$age == 1L] == "age1"),
  all(catch$sd_block[catch$age >= 2L] == "age2plus"),
  length(unique(index$q_key)) == 10L,
  all(!is.na(index$sd_block)),
  nrow(obs$pop_covariates) == length(years) * length(ages),
  setequal(levels(obs$pop_covariates$f_age_group), c("age0", "age1", "age2", "age3_4", "age5plus")),
  all(obs$maturity$spawn_time == 0)
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
  prepared$F_settings$process == "iid",
  prepared$M_settings$process == "off"
)

cat("Western Baltic cod tinyAM observations and model preparation passed.\n")
