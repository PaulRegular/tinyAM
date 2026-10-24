pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_data <- read_assessment("ices_saithe_faroe_2025", database)
env <- new.env(parent = environment())
env$source <- source_data
env$do_fit <- FALSE
env$silent <- TRUE
sys.source(
  "analysis/comp_assessments/scripts/translation/stocks/ices_saithe_faroe_2025.R",
  envir = env
)

prepared <- tinyAM::prepare_tam(
  data = env$obs,
  years = env$years,
  ages = env$ages,
  N_settings = list(process = "iid", init = "exp", sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 4:8),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + sd_scale, sd_supplied = ~relative_sd),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
)

stopifnot(
  identical(env$years, 1961:2025),
  identical(env$ages, 3:15),
  is.null(env$fit),
  nrow(env$obs$catch) == 845L,
  nrow(env$obs$index) == 496L,
  env$obs$catch$obs[env$obs$catch$year == 1962 & env$obs$catch$age == 3] == 562,
  env$obs$catch$obs[env$obs$catch$year == 1961 & env$obs$catch$age == 4] == 379,
  all(env$obs$catch$unit == "thousand fish"),
  sum(env$obs$catch$obs == 0) == 28L,
  all(env$obs$catch$relative_sd[env$obs$catch$year < 2025] == 1),
  isTRUE(all.equal(
    env$obs$catch$relative_sd[env$obs$catch$year == 2025],
    rep(sqrt(2), sum(env$obs$catch$year == 2025))
  )),
  ncol(prepared$sd_catch_modmat) == 1L,
  ncol(prepared$q_modmat) == 14L,
  ncol(prepared$sd_index_modmat) == 16L,
  length(unique(env$obs$index$q_key)) == 14L,
  length(unique(env$obs$index$sd_block)) == 16L,
  abs(env$obs$index$obs[env$obs$index$survey == "SummerSurvey" &
    env$obs$index$year == 1996 & env$obs$index$age == 3] - 0.001425) < 1e-12,
  abs(env$obs$index$obs[env$obs$index$survey == "SpringSurvey" &
    env$obs$index$year == 1994 & env$obs$index$age == 3] - 0.00118) < 1e-12,
  abs(env$obs$index$obs[env$obs$index$survey == "SummerSurvey" &
    env$obs$index$year == 1997 & env$obs$index$age == 3] - 0.006105) < 1e-12,
  abs(env$obs$index$obs[env$obs$index$survey == "SpringSurvey" &
    env$obs$index$year == 1994 & env$obs$index$age == 4] - 0.0082) < 1e-12,
  all(env$obs$index$samp_time[env$obs$index$survey == "SummerSurvey"] == 0.65),
  all(env$obs$index$samp_time[env$obs$index$survey == "SpringSurvey"] == 0.225),
  all(env$obs$pop_covariates$M_assumption == 0.2),
  all(is.finite(env$obs$pop_covariates$weight)),
  all(is.finite(env$obs$catch$weight)),
  env$obs$pop_covariates$weight[env$obs$pop_covariates$year == 2025,
    env$obs$pop_covariates$age == 3] ==
    source_data$outputs$value[source_data$outputs$type == "biology" &
      source_data$outputs$measure == "weight_at_age" &
      source_data$outputs$year == 2025 & source_data$outputs$age == 3],
  env$age_plus_group == 15L,
  identical(env$comparison_aggregates,
    c("abundance", "biomass", "ssb", "recruitment", "F_bar")),
  env$comparison_definitions$F$status == "approximate",
  env$comparison_definitions$F_bar$status == "approximate"
)

native_index <- source_data$outputs[
  source_data$outputs$type == "index" &
    source_data$outputs$measure == "predicted_index",
  , drop = FALSE
]
dashboard_index <- env$comparison_outputs[
  env$comparison_outputs$type == "index" &
    env$comparison_outputs$measure == "predicted_index",
  , drop = FALSE
]
stopifnot(
  nrow(dashboard_index) == nrow(native_index),
  all(dashboard_index$unit == "thousand fish per effort-hour"),
  isTRUE(all.equal(dashboard_index$value, native_index$value * 0.001)),
  isTRUE(all.equal(dashboard_index$se, native_index$se * 0.001)),
  isTRUE(all.equal(dashboard_index$lwr, native_index$lwr * 0.001)),
  isTRUE(all.equal(dashboard_index$upr, native_index$upr * 0.001)),
  isTRUE(all.equal(source_data$outputs$value[
    source_data$outputs$type == "index" &
      source_data$outputs$measure == "predicted_index"
  ], native_index$value))
)

native_q <- source_data$outputs[
  source_data$outputs$type == "catchability" &
    source_data$outputs$measure == "q",
  , drop = FALSE
]
dashboard_q <- env$comparison_outputs[
  env$comparison_outputs$type == "catchability" &
    env$comparison_outputs$measure == "q",
  , drop = FALSE
]
stopifnot(
  nrow(dashboard_q) == nrow(native_q),
  all(dashboard_q$unit == "per effort-hour"),
  isTRUE(all.equal(dashboard_q$value, native_q$value * 0.001)),
  isTRUE(all.equal(dashboard_q$se, native_q$se * 0.001)),
  isTRUE(all.equal(dashboard_q$lwr, native_q$lwr * 0.001)),
  isTRUE(all.equal(dashboard_q$upr, native_q$upr * 0.001)),
  isTRUE(all.equal(source_data$outputs$value[
    source_data$outputs$type == "catchability" &
      source_data$outputs$measure == "q"
  ], native_q$value))
)

cat("Faroe saithe 2025 tinyAM translation preparation checks passed.\n")
