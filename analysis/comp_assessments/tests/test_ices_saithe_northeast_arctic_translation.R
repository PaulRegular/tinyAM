source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "database_to_tam_ref.R"))

assessment_id <- "ices_saithe_northeast_arctic_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
translated <- .test_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

inputs <- source_data$inputs
outputs <- source_data$outputs
precision <- inputs[
  inputs$type == "index" & inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
key <- function(x) paste(x$survey, x$year, x$age, sep = "\034")
precision_key <- key(precision)
weighted <- key(translated$obs$index) %in% precision_key
weighted_key <- key(translated$obs$index[weighted, ])
source_weight <- as.numeric(precision$value[match(weighted_key, precision_key)])

ref <- database_to_tam_ref(
  assessment_id,
  outputs,
  obs = translated$obs,
  years = translated$years,
  ages = translated$ages,
  age_plus_group = 12L,
  assumptions = source_data$assumptions,
  comparison_aggregates = translated$comparison_aggregates,
  comparison_definitions = translated$comparison_definitions
)

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1960:2025),
  identical(translated$ages, 3:12),
  identical(translated$age_plus_group, 12L),
  nrow(translated$obs$catch) == 660L,
  nrow(translated$obs$index) == 160L,
  nrow(translated$obs$weight) == 660L,
  nrow(translated$obs$maturity) == 660L,
  identical(as.integer(table(translated$obs$index$survey)), c(40L, 120L)),
  all(translated$obs$index$samp_time == 0.8),
  all(is.finite(translated$obs$index$relative_sd)),
  all(translated$obs$index$relative_sd > 0),
  isTRUE(all.equal(source_weight, translated$obs$index$relative_precision_weight[weighted])),
  isTRUE(all.equal(
    translated$obs$index$relative_sd[weighted], 1 / sqrt(source_weight)
  )),
  all(translated$obs$index$relative_precision_weight[!weighted] == 1),
  all(translated$obs$pop_covariates$M_assumption == 0.2),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "free"),
  identical(dat$rec_settings$process, "iid"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$F_settings$mean_ages, 4:7),
  identical(dat$M_settings$process, "off"),
  is.null(translated$start_par),
  !"start_par" %in% names(translated$fit_call),
  identical(dim(par$log_f), c(66L, 10L)),
  all(is.finite(par$log_f)),
  nrow(translated$obs$index[translated$obs$index$q_age == "age6_7", ]) == 64L,
  length(unique(translated$obs$index$q_key)) == 8L,
  sum(outputs$type == "population" & outputs$measure == "numbers_at_age") == 660L,
  sum(outputs$type == "mortality" & outputs$measure == "fishing_mortality_at_age") == 660L,
  sum(outputs$type == "catch" & outputs$measure == "predicted_catch") == 660L,
  sum(outputs$type == "index" & outputs$measure == "predicted_index") == 160L,
  sum(outputs$type == "catchability" & outputs$measure == "q") == 8L,
  all(is.finite(outputs$se[outputs$measure %in% c("numbers_at_age", "fishing_mortality_at_age")])),
  all(is.finite(ref$pop$N$est)),
  all(is.finite(ref$pop$F$est)),
  any(ref$pop$N$year == 2025L),
  any(outputs$year == 2026L & outputs$measure == "SSB")
)

cat("Northeast Arctic saithe translation structure passed.\n")
