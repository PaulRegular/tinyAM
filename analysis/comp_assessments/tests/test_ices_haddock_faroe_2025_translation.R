pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
assessment <- read_assessment("ices_haddock_faroe_2025", database)
env <- new.env(parent = environment())
env$source <- assessment
env$do_fit <- FALSE
env$silent <- TRUE
sys.source(
  "analysis/comp_assessments/scripts/translation/stocks/ices_haddock_faroe_2025.R",
  envir = env
)
model_data <- tinyAM::prepare_tam(
  data = env$obs,
  years = env$years,
  ages = env$ages,
  N_settings = list(process = "iid", init = "exp", sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:7),
  M_settings = list(
    process = "off", mu_form = NULL,
    mu_supplied = ~M_assumption
  ),
  catch_settings = list(sd_form = ~0 + sd_scale, sd_supplied = ~relative_sd),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_group)
)

stopifnot(
  is.null(env$fit),
  identical(env$years, 1957:2025),
  identical(env$ages, 1:10),
  nrow(env$obs$catch) == 690L,
  sum(env$obs$catch$obs == 0) == 41L,
  all(env$obs$catch$unit == "thousand fish"),
  all(env$obs$catch$relative_sd[env$obs$catch$year < 2025] == 1),
  isTRUE(all.equal(
    env$obs$catch$relative_sd[env$obs$catch$year == 2025],
    rep(sqrt(2), sum(env$obs$catch$year == 2025))
  )),
  ncol(model_data$sd_catch_modmat) == 1L,
  ncol(model_data$q_modmat) == 13L,
  ncol(model_data$sd_index_modmat) == 4L,
  length(unique(env$obs$index$q_key)) == 13L,
  length(unique(env$obs$index$sd_group)) == 4L,
  all(env$obs$index$samp_time[env$obs$index$survey == "FO-GFS-Q3"] == 0.65),
  all(env$obs$index$samp_time[env$obs$index$survey == "FO-GFS-Q1"] == 0.225),
  length(unique(env$comparison_aggregates)) == 4L
)

cat("Faroe haddock 2025 tinyAM translation preparation checks passed.\n")
