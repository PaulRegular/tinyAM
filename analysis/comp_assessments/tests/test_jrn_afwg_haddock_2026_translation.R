pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "jrn_afwg_haddock_northeast_arctic_2026"
source_data <- read_assessment(assessment_id, read_database())
env <- new.env(parent = globalenv())
env$source <- source_data
env$do_fit <- FALSE
env$silent <- TRUE
env$print_sources <- function(x) character()
sys.source(file.path(root, "scripts", "translation", "stocks", paste0(assessment_id, ".R")), env)

settings <- list(
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(
    process = "cor_rw", mu_form = NULL,
    sd_form = ~0 + F_sd_group, mean_ages = 4:7
  ),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + catch_sd_group),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey)
)
prepared <- do.call(tinyAM::prepare_tam, c(
  list(data = env$obs, years = env$years, ages = env$ages), settings
))
index_counts <- table(env$obs$index$survey)
q_counts <- tapply(env$obs$index$q_key, env$obs$index$survey,
  function(x) length(unique(x)))
expected_counts <- c(
  "RU-BTr-Q4" = 130L, "BS-NoRU-Q1(Aco)" = 198L,
  "BS-NoRu-Q1(BTr)" = 230L, "FLT007:Eco-NoRu-Q3(Btr)" = 120L
)
expected_times <- c(
  "RU-BTr-Q4" = 0.95, "BS-NoRU-Q1(Aco)" = 0.133,
  "BS-NoRu-Q1(BTr)" = 0.133, "FLT007:Eco-NoRu-Q3(Btr)" = 0.7
)

stopifnot(
  is.null(env$fit),
  identical(env$years, 1950:2026),
  identical(env$ages, 3:13),
  nrow(env$obs$catch) == 847L,
  nrow(env$obs$index) == 678L,
  all(index_counts[names(expected_counts)] == expected_counts),
  identical(sort(names(index_counts)), sort(names(expected_times))),
  all(vapply(names(expected_times), function(survey) {
    all(env$obs$index$samp_time[env$obs$index$survey == survey] == expected_times[[survey]])
  }, logical(1))),
  identical(as.integer(q_counts[names(expected_times)]), c(2L, 3L, 3L, 2L)),
  all(is.finite(env$obs$pop_covariates$M_assumption)),
  all(env$obs$pop_covariates$M_assumption[env$obs$pop_covariates$age >= 7] == 0.2),
  identical(levels(env$obs$catch$catch_sd_group), c("age3", "age4plus")),
  identical(levels(env$obs$pop_covariates$F_sd_group), c("age3", "age4plus")),
  ncol(prepared$process_sd$F$matrix) == 2L,
  identical(prepared$F_settings$mean_ages, 4:7),
  env$age_plus_group == 13L,
  identical(env$comparison_aggregates, "ssb"),
  grepl("does not share the official ages 8-13+ F state",
    paste(env$background, collapse = " "), fixed = TRUE)
)

cat("JRN-AFWG Northeast Arctic haddock translation preparation passed.\n")
