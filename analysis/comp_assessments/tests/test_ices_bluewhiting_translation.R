source("analysis/comp_assessments/tests/helper_stock.R")
pkgload::load_all(".", quiet = TRUE)
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_bluewhiting_northeast_atlantic_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
before <- source_data
translated <- .test_stock(source_data)
source_q <- source_data$outputs[
  source_data$outputs$type == "catchability" & source_data$outputs$measure == "q", , drop = FALSE
]
expected_q <- vapply(c("1", "2", "3", "4", "5-8"), function(group) {
  rows <- if (group == "5-8") source_q[source_q$age_group %in% group, , drop = FALSE] else
    source_q[as.character(source_q$age) %in% group, , drop = FALSE]
  unique(rows$value) * 1000
}, numeric(1))
names(expected_q) <- names(translated$start_par$log_q)

dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

stopifnot(
  tinyAM::check_tam_data(translated$obs),
  identical(translated$years, 1981:2026),
  identical(translated$ages, 1:10),
  nrow(translated$obs$catch) == 460L,
  nrow(translated$obs$index) == 168L,
  all(translated$obs$index$samp_time == 0.245),
  identical(dat$N_settings$process, "off"),
  identical(dim(translated$start_par$log_f), dim(par$log_f)),
  identical(dim(translated$start_par$log_n0), dim(par$log_n0)),
  identical(names(translated$start_par$log_q), names(par$log_q)),
  identical(names(translated$start_par$log_sd_catch), names(par$log_sd_catch)),
  identical(names(translated$start_par$log_sd_index), names(par$log_sd_index)),
  length(unique(translated$obs$index$q_key)) == 5L,
  setequal(levels(translated$obs$index$sd_group), c("1", "2", "3", "4-6", "7-8")),
  isTRUE(all.equal(exp(translated$start_par$log_q), expected_q)),
  identical(source_data, before),
  all(is.finite(translated$start_par$log_q))
)

cat("Blue whiting translation structure passed.\n")
