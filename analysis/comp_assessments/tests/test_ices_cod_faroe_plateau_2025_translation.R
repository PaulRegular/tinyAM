root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

id <- "ices_cod_faroe_plateau_2025"
source_data <- read_assessment(id, read_database())
env <- new.env(parent = globalenv())
env$source <- source_data
env$do_fit <- FALSE
env$silent <- TRUE
env$print_sources <- function(x) character()
sys.source(file.path(root, "scripts", "translation", "stocks", paste0(id, ".R")), env)

stopifnot(
  identical(env$years, 1959:2025),
  identical(env$ages, 1:10),
  is.null(env$fit),
  sum(is.finite(env$obs$catch$obs)) == 590L,
  sum(is.finite(env$obs$index$obs)) == 550L,
  setequal(unique(env$obs$index$survey), c("FO-GFS-Q3", "FO-GFS-Q1")),
  all(env$obs$index$samp_time[env$obs$index$survey == "FO-GFS-Q3"] == 0.65),
  all(env$obs$index$samp_time[env$obs$index$survey == "FO-GFS-Q1"] == 0.225),
  length(unique(env$obs$index$q_key[env$obs$index$survey == "FO-GFS-Q3"])) == 8L,
  length(unique(env$obs$index$q_key[env$obs$index$survey == "FO-GFS-Q1"])) == 8L,
  length(unique(env$obs$index$sd_block[env$obs$index$survey == "FO-GFS-Q3"])) == 2L,
  length(unique(env$obs$index$sd_block[env$obs$index$survey == "FO-GFS-Q1"])) == 2L,
  setequal(unique(env$obs$catch$sd_block), c("ages2_9", "age10_plus")),
  all(env$obs$pop_covariates$M_assumption == 0.2),
  env$age_plus_group == 10L,
  identical(env$comparison_aggregates, c("abundance", "biomass", "ssb")),
  env$comparison_definitions$F$status == "approximate"
)

prepared <- tinyAM::prepare_tam(
  data = env$obs,
  years = env$years,
  ages = env$ages,
  N_settings = list(process = "iid", init = "exp", sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:7),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + sd_block),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
)

stopifnot(
  ncol(prepared$process_sd$N$matrix) == 1L,
  all(is.finite(prepared$obs$catch$obs[is.finite(prepared$obs$catch$obs)])),
  all(is.finite(prepared$obs$index$obs[is.finite(prepared$obs$index$obs)]))
)

cat("Faroe Plateau cod tinyAM observation translation and preparation checks passed.\n")
