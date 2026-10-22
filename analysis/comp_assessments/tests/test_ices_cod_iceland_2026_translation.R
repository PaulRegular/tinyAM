root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source_data <- read_assessment("ices_cod_iceland_2026", read_database())
env <- new.env(parent = globalenv())
env$source <- source_data
env$do_fit <- FALSE
env$silent <- TRUE
env$print_sources <- function(x) character()
sys.source(file.path(root, "scripts", "translation", "stocks", "ices_cod_iceland_2026.R"), env)

stopifnot(
  identical(env$years, 1955:2026),
  identical(env$ages, 3:14),
  is.null(env$fit),
  sum(is.finite(env$obs$catch$obs)) == 852L,
  all(is.na(env$obs$catch$obs[env$obs$catch$year == 2026])),
  nrow(env$obs$index) == 823L,
  !any(env$obs$index$survey == "IS-SMH" & env$obs$index$age == 14),
  !any(env$obs$index$survey == "IS-SMH" & env$obs$index$year == 2011),
  all(env$obs$index$samp_time[env$obs$index$survey == "IS-SMB"] == 0.2),
  all(env$obs$index$samp_time[env$obs$index$survey == "IS-SMH"] == 0.8),
  length(unique(env$obs$index$q_key)) == 23L,
  setequal(unique(env$obs$index$sd_block), c("IS-SMB", "IS-SMH")),
  all(env$obs$pop_covariates$M_assumption == 0.2),
  all(env$obs$weight$obs > 0),
  is.null(env$age_plus_group),
  env$comparison_definitions$ssb$status == "approximate"
)

prepared <- tinyAM::prepare_tam(
  data = env$obs,
  years = env$years,
  ages = env$ages,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "rw", mu_form = NULL, mean_ages = 5:10),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
)
stopifnot(
  identical(prepared$years, env$years),
  identical(prepared$ages, env$ages),
  all(is.finite(prepared$obs$catch$obs) | is.na(prepared$obs$catch$obs)),
  all(is.finite(prepared$obs$index$obs) | is.na(prepared$obs$index$obs)),
  all(prepared$obs$index$obs >= 0)
)

cat("Icelandic cod tinyAM observation conversion and model preparation checks passed.\n")
