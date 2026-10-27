root <- file.path("analysis", "comp_assessments")
source_data <- read.csv(file.path(root, "database", "inputs.csv"),
  na.strings = c("", "NA"), stringsAsFactors = FALSE
)
outputs <- read.csv(file.path(root, "database", "outputs.csv"),
  na.strings = c("", "NA"), stringsAsFactors = FALSE
)
assessments <- read.csv(file.path(root, "database", "assessments.csv"),
  na.strings = c("", "NA"), stringsAsFactors = FALSE
)
assumptions <- read.csv(file.path(root, "database", "assumptions.csv"),
  na.strings = c("", "NA"), stringsAsFactors = FALSE
)
id <- "dfo_cod_3ps_2024"
script_path <- file.path(root, "scripts", "translation", "stocks", paste0(id, ".R"))
script_lines <- readLines(script_path, warn = FALSE)
invisible(parse(file = script_path))
stopifnot(
  any(grepl("tinyAM::fit_tam", script_lines, fixed = TRUE)),
  any(grepl("if (do_fit)", script_lines, fixed = TRUE)),
  any(grepl("fit_default_first", script_lines, fixed = TRUE)),
  any(grepl("start_par = start_par", script_lines, fixed = TRUE)),
  any(grepl("cor_rw", script_lines, fixed = TRUE)),
  any(grepl("M_condition", script_lines, fixed = TRUE)),
  any(grepl('process = "iid"', script_lines, fixed = TRUE))
)

source_record <- list(
  assessment = assessments[assessments$assessment_id == id, , drop = FALSE],
  inputs = source_data[source_data$assessment_id == id, , drop = FALSE],
  outputs = outputs[outputs$assessment_id == id, , drop = FALSE],
  assumptions = assumptions[assumptions$assessment_id == id, , drop = FALSE]
)
env <- new.env(parent = globalenv())
env$source <- source_record
env$do_fit <- FALSE
env$silent <- TRUE
env$print_sources <- function(...) character()
env$database_to_tam_obs <- function(assessment_id, inputs, years, ages,
                                    sampling_times, surveys, assumptions) {
  env$converter_call <- list(
    assessment_id = assessment_id, years = years, ages = ages,
    sampling_times = sampling_times, surveys = surveys
  )
  list(pop_covariates = expand.grid(year = years, age = ages))
}
sys.source(script_path, envir = env)
expected_times <- c(
  "DFO RV Spring" = 0.25, "ERHAPS" = 0.17, "GEAC" = 0.875,
  "Sentinel Gillnet" = 0.5, "Sentinel Line-trawl" = 0.5
)
stopifnot(
  identical(env$converter_call$years, 1959:2024),
  identical(env$converter_call$ages, 2:14),
  identical(env$converter_call$surveys, names(expected_times)),
  identical(env$converter_call$sampling_times, expected_times),
  nrow(env$obs$pop_covariates) == 66L * 13L,
  all(is.finite(env$obs$pop_covariates$M_condition)),
  identical(env$N_settings$process, "iid"),
  all(env$obs$pop_covariates$M_condition_age_2_5[
    env$obs$pop_covariates$age >= 6
  ] == 0),
  all(env$obs$pop_covariates$M_condition_age_6_plus[
    env$obs$pop_covariates$age <= 5
  ] == 0),
  all(env$obs$pop_covariates$M_condition_age_2_5[
    env$obs$pop_covariates$age <= 5
  ] == env$obs$pop_covariates$M_condition[
    env$obs$pop_covariates$age <= 5
  ]),
  all(env$obs$pop_covariates$M_condition_age_6_plus[
    env$obs$pop_covariates$age >= 6
  ] == env$obs$pop_covariates$M_condition[
    env$obs$pop_covariates$age >= 6
  ]),
  setequal(unique(env$comparison_outputs$year), 2020:2024),
  !any(env$comparison_outputs$year == 2025),
  identical(env$comparison_aggregates, "ssb"),
  is.null(env$fit)
)
condition_design <- stats::model.matrix(
  env$M_settings$mu_form,
  data = env$obs$pop_covariates
)
stopifnot(identical(
  grep("^M_condition_age_", colnames(condition_design), value = TRUE),
  c("M_condition_age_2_5", "M_condition_age_6_plus")
))

start_years <- 2020:2024
starts <- list(
  log_r = setNames(numeric(length(1960:2024)), as.character(1960:2024)),
  log_f = matrix(0, length(1959:2024), length(2:14),
    dimnames = list(as.character(1959:2024), as.character(2:14))
  )
)
source_starts <- env$source_start_par(starts, data.frame(
  measure = c(rep("recruitment", 5), rep("Fbar", 5)),
  year = c(start_years, start_years),
  age = c(rep(2, 5), rep(NA, 5)),
  value = c(17.4, 20, 34.3, 36.8, 36.6, 0.04, 0.02, 0.03, 0.03, 0.04),
  unit = c(rep("million fish", 5), rep("per year", 5))
))
stopifnot(
  isTRUE(all.equal(
    unname(source_starts$log_r[as.character(start_years)]),
    log(c(17.4, 20, 34.3, 36.8, 36.6) * 1000)
  )),
  all(source_starts$log_r[as.character(1960:2019)] == 0),
  all(source_starts$log_f[as.character(start_years), as.character(5:8)] ==
    matrix(log(c(0.04, 0.02, 0.03, 0.03, 0.04)), 5L, 4L)),
  all(source_starts$log_f[as.character(start_years), as.character(c(2:4, 9:14))] == 0)
)
cat("3Ps translation recipe years, ages, IID N process, grouped condition effects, comparison window, and fit guard passed.\n")
