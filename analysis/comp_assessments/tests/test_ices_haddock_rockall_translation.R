root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
.translation_abort <- function(message) stop(message, call. = FALSE)
id <- "ices_haddock_rockall_2026"
source_data <- read_assessment(id, read_database())
env <- new.env(parent = globalenv())
env$source <- source_data
env$do_fit <- FALSE
env$silent <- TRUE
env$print_sources <- function(x) character()
sys.source(file.path(root, "scripts", "translation", "stocks", paste0(id, ".R")), env)

model_data <- tinyAM::prepare_tam(
  data = env$obs,
  years = env$years,
  ages = env$ages,
  N_settings = list(process = "iid", init = "exp", sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(
    q_form = ~0 + q_key,
    sd_form = ~0 + survey,
    sd_supplied = ~relative_sd
  )
)

stopifnot(
  identical(env$years, 1991:2025), identical(env$ages, 1:8),
  is.null(env$fit),
  nrow(env$obs$catch) == 280L,
  nrow(env$obs$weight) == 280L,
  nrow(env$obs$maturity) == 280L,
  nrow(env$obs$index) == 247L,
  sum(env$obs$index$survey == "Rockall pt1") == 55L,
  sum(env$obs$index$survey == "Modelled Q3") == 192L,
  all(is.finite(env$obs$index$relative_sd)),
  all(env$obs$index$relative_sd[env$obs$index$survey == "Rockall pt1"] == 1),
  length(unique(env$obs$index$q_key[env$obs$index$survey == "Rockall pt1"])) == 5L,
  length(unique(env$obs$index$q_key[env$obs$index$survey == "Modelled Q3"])) == 5L,
  !"N_sd_key" %in% names(env$obs$pop_covariates),
  ncol(model_data$process_sd$N$matrix) == 1L,
  all(env$obs$pop_covariates$M_assumption > 0),
  env$age_plus_group == 8L,
  identical(env$comparison_aggregates, "ssb"),
  any(attr(env$obs, "translation")$source_provenance$source_type == "analysis_approximation")
)
cat("Rockall haddock translation observations, q groups, precision and fit guard passed.\n")
