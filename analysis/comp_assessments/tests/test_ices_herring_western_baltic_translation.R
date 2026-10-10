source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_herring_western_baltic_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
translated <- .test_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

input <- source_data$inputs
outputs <- source_data$outputs
source_surface <- function(type, measure, scale = 1) {
  rows <- outputs[outputs$type == type & outputs$measure == measure &
                    !is.na(outputs$year) & !is.na(outputs$age), ]
  surface <- matrix(
    NA_real_, length(translated$years), length(translated$ages),
    dimnames = list(as.character(translated$years),
                    as.character(translated$ages))
  )
  index <- cbind(match(as.integer(rows$year), translated$years),
                 match(as.integer(rows$age), translated$ages))
  surface[index] <- as.numeric(rows$value) * scale
  surface
}
source_N <- source_surface("population", "numbers_at_age", 1000)
source_F <- source_surface("mortality", "fishing_mortality_at_age")
source_M <- source_surface("mortality", "natural_mortality_at_age")

catch_source <- input[input$type == "catch" &
                        input$measure == "numbers_at_age", ]
catch_source_key <- paste(as.integer(catch_source$year),
                         as.integer(catch_source$age))
catch_obs_key <- paste(translated$obs$catch$year,
                       translated$obs$catch$age)
precision <- input[input$measure == "relative_precision_weight", ]
precision$year <- as.integer(precision$year)
precision$age <- as.integer(precision$age)
precision$value <- as.numeric(precision$value)

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1991:2025),
  identical(translated$ages, 0:8),
  nrow(translated$obs$catch) == 315L,
  nrow(translated$obs$index) == 372L,
  nrow(translated$obs$weight) == 315L,
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "HERAS"
  ] == 0.625),
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "GERAS"
  ] == 0.8),
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "N20"
  ] == 0.4),
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "IBTS/BITSQ1"
  ] == 0.136365),
  length(unique(translated$obs$index$q_key)) == 8L,
  length(unique(translated$obs$index$sd_block)) == 4L,
  length(unique(translated$obs$catch$sd_block)) == 3L,
  !any(translated$obs$index$survey == "HERAS" &
         translated$obs$index$year == 1999),
  all(translated$obs$index$relative_sd[
    translated$obs$index$survey != "IBTS/BITSQ1"
  ] == 1),
  max(abs(
    translated$obs$index$relative_sd[
      translated$obs$index$survey == "IBTS/BITSQ1"
    ] -
      1 / sqrt(precision$value[match(
        paste(translated$obs$index$year[
          translated$obs$index$survey == "IBTS/BITSQ1"
        ], translated$obs$index$age[
          translated$obs$index$survey == "IBTS/BITSQ1"
        ]),
        paste(precision$year, precision$age)
      )])
  )) < 1e-12,
  !anyNA(match(catch_obs_key, catch_source_key)),
  isTRUE(all.equal(
    translated$obs$catch$obs,
    as.numeric(catch_source$value[match(catch_obs_key, catch_source_key)]) * 1000
  )),
  all(is.finite(source_N)), all(is.finite(source_F)), all(is.finite(source_M)),
  isTRUE(all.equal(translated$start_par$log_r0, log(source_N[1, "0"]))),
  isTRUE(all.equal(
    unname(translated$start_par$log_r),
    unname(log(source_N[-1, "0"]))
  )),
  isTRUE(all.equal(
    unname(translated$start_par$log_n),
    unname(log(source_N[-1, as.character(1:8), drop = FALSE]))
  )),
  isTRUE(all.equal(
    unname(translated$start_par$log_f),
    unname(log(source_F))
  )),
  identical(dim(translated$start_par$log_f), dim(par$log_f)),
  identical(dim(translated$start_par$log_n), dim(par$log_n)),
  identical(names(translated$start_par$log_q), names(par$log_q)),
  all(is.finite(translated$start_par$log_q)),
  identical(dat$N_settings$process, "rw"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$F_settings$process, "rw"),
  identical(dat$F_settings$mean_ages, 2:5),
  identical(dat$M_settings$process, "off"),
  grepl("AR(1)", paste(translated$background, collapse = " "),
        fixed = TRUE)
)

cat("Western Baltic herring translation structure passed.\n")
