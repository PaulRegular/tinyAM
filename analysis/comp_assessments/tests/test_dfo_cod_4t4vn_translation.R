pkgload::load_all(".", quiet = TRUE)

root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "database_to_tam_ref.R"))

stock_env <- new.env(parent = environment())
sys.source(
  file.path(root, "scripts", "translation", "stocks", "dfo_cod_4t4vn_2019.R"),
  envir = stock_env
)
source_data <- read_assessment("dfo_cod_4t4vn_2019", read_database())
landings_at_age <- source_data$inputs[
  source_data$inputs$type == "catch" &
    source_data$inputs$measure == "landings_numbers_at_age",
  , drop = FALSE
]
stopifnot(nrow(landings_at_age) == 480L)
stopifnot(!any(source_data$inputs$type == "catch" &
                 source_data$inputs$measure == "numbers_at_age"))
translated <- stock_env$translate_stock(source_data)
obs <- translated$obs

stopifnot(identical(translated$years, 1971:2018))
stopifnot(identical(translated$ages, 2:12))
stopifnot(check_obs(obs))
expected_key <- paste(landings_at_age$year, landings_at_age$age, sep = ":")
observed_catch <- obs$catch[!is.na(obs$catch$obs), , drop = FALSE]
catch_key <- paste(observed_catch$year, observed_catch$age, sep = ":")
expected_rows <- match(catch_key, expected_key)
stopifnot(!anyNA(expected_rows))
stopifnot(nrow(observed_catch) == 480L)
stopifnot(all(observed_catch$obs == landings_at_age$value[expected_rows] * 1000))
stopifnot(setequal(unique(obs$index$survey), c(
  "DFO September RV survey", "Mobile Sentinel August survey"
)))
stopifnot(!any(obs$index$survey == "DFO September RV survey" &
                 obs$index$year %in% c(1980, 1985, 2003)))
stopifnot(all(obs$index$samp_time[obs$index$survey == "DFO September RV survey"] == 0.75))
stopifnot(all(obs$index$samp_time[obs$index$survey == "Mobile Sentinel August survey"] == 0.625))
stopifnot(identical(
  obs$weight$obs[obs$weight$age == 12],
  obs$weight$obs[obs$weight$age == 11]
))
stopifnot(all(obs$weight$M_prior_mean[obs$weight$age <= 4] == 0.65))
stopifnot(all(obs$weight$M_prior_mean[obs$weight$age %in% 5:8] == 0.15))
stopifnot(all(obs$weight$M_prior_mean[obs$weight$age >= 9] == 0.15))
stopifnot(identical(translated$settings$M_settings$process, "rw"))
stopifnot(identical(translated$settings$M_settings$age_breaks, c(2, 5, 9, 12)))
stopifnot(identical(translated$settings$M_settings$first_dev_year, 1971L))
dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)
stopifnot(identical(rownames(par$log_m), as.character(1971:2018)))
stopifnot(identical(colnames(par$log_m), c("2-4", "5-8", "9-12")))
stopifnot(isTRUE(all.equal(
  unname(par$log_m[1, ]), log(c(0.65, 0.15, 0.15))
)))
stopifnot(identical(
  translated$warm_start_settings$M_settings$process, "iid"
))
stopifnot(all(is.finite(obs$weight$obs[obs$weight$year %in% c(1980, 1985)])))
stopifnot(any(grepl("translation_assumption",
                    attr(obs, "translation")$source_provenance$source_type)))
catch_provenance <- attr(obs, "translation")$source_provenance
stopifnot(any(catch_provenance$component == "catch" &
                catch_provenance$source_type == "translation_assumption"))
stopifnot(setequal(unique(translated$comparison_outputs$measure), c(
  "numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age",
  "SSB", "recruitment"
)))
ref <- database_to_tam_ref(
  "dfo_cod_4t4vn_2019", translated$comparison_outputs,
  obs = obs, years = translated$years, ages = translated$ages,
  terminal_year = 2018, age_plus_group = 12
)
stopifnot(setequal(ref$pop$M$age, 5:12))
stopifnot(all(ref$pop$M$year == 2018L))
stopifnot(all(ref$pop$M$est[ref$pop$M$age %in% 5:8] == 0.81))
stopifnot(all(ref$pop$M$est[ref$pop$M$age %in% 9:12] == 0.85))

readiness_env <- new.env(parent = globalenv())
sys.source(
  file.path(root, "scripts", "database", "003_observation_readiness.R"),
  envir = readiness_env
)
readiness <- utils::read.csv(
  file.path(root, "results", "observation_readiness.csv"),
  stringsAsFactors = FALSE
)
southern_gulf <- readiness[
  readiness$assessment_id == "dfo_cod_4t4vn_2019", , drop = FALSE
]
stopifnot(nrow(southern_gulf) == 1L)
stopifnot(southern_gulf$catch_at_age_rows == 0L)
stopifnot(southern_gulf$landings_at_age_rows == 480L)
stopifnot(isTRUE(southern_gulf$database_to_tam_obs_succeeds))
stopifnot(isTRUE(southern_gulf$tinyAM_check_obs_passes))
