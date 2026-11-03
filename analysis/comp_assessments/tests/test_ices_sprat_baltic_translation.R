source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_sprat_baltic_2026"
source_data <- read_assessment(assessment_id, read_database())
translated <- .test_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- translated$start_par

inputs <- source_data$inputs
outputs <- source_data$outputs
catch_1974 <- subset(inputs, type == "catch" & year == 1974 & age == 1)
bass_2001 <- subset(inputs, type == "index" & survey == "BASS_May_SD24_26_28" &
  year == 2001 & age == 1)
m_1974 <- subset(inputs, type == "M" & year == 1974 & age == 1)
ssb_2025 <- subset(outputs, measure == "SSB" & year == 2025)
recruitment_2026 <- subset(outputs, measure == "recruitment" & year == 2026)
total_numbers <- subset(outputs, measure == "total_numbers")
n_at_age <- subset(outputs, measure == "numbers_at_age")
expected_total_numbers <- tapply(
  as.numeric(n_at_age$value), n_at_age$year, sum
)

stopifnot(
  tinyAM::check_obs(translated$obs),
  all(translated$obs$maturity$spawn_time == .4),
  identical(translated$years, 1974:2025),
  identical(translated$ages, 1:8),
  source_data$assessment$terminal_year == 2025,
  source_data$assessment$estimate_terminal_year == 2026,
  nrow(subset(inputs, type == "catch")) == 52L * 8L,
  nrow(subset(inputs, type == "weight")) == 52L * 8L,
  nrow(subset(inputs, type == "catch_weight")) == 52L * 8L,
  nrow(subset(inputs, type == "M")) == 52L * 8L,
  nrow(subset(inputs, type == "maturity")) == 8L,
  nrow(subset(inputs, type == "biology")) == 16L,
  nrow(subset(inputs, type == "index")) == 465L,
  nrow(translated$obs$index) == 464L,
  !any(translated$obs$index$year == 2026),
  nlevels(translated$obs$index$q_key) == 19L,
  all(is.finite(translated$obs$pop_covariates$M_assumption)),
  identical(as.numeric(catch_1974$value), 2854471),
  identical(as.numeric(bass_2001$value), 8225),
  identical(as.numeric(m_1974$value), 0.75),
  identical(as.numeric(ssb_2025$value), 620553),
  identical(as.numeric(recruitment_2026$value), 123196),
  nrow(total_numbers) == 53L,
  identical(as.integer(total_numbers$year), 1974:2026),
  identical(total_numbers$source_type,
    rep("derived_source_output", 53L)),
  all(total_numbers$unit == "million fish"),
  all(is.na(total_numbers$age)),
  isTRUE(all.equal(
    as.numeric(total_numbers$value),
    as.numeric(expected_total_numbers[as.character(total_numbers$year)])
  )),
  any(grepl("2026 is an intermediate-year estimate", total_numbers$notes)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$F_settings$process, "rw"),
  identical(dat$F_settings$mean_ages, 3:5),
  identical(dim(par$log_f), c(52L, 8L)),
  is.null(par$log_n0),
  all(is.finite(par$log_f)),
  all(is.finite(par$log_q)),
  isTRUE(all.equal(par$log_r0, log(76636 * 1000)))
)

cat("Baltic sprat translation structure passed.\n")
