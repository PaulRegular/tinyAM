root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_plaice_north_sea_2026"
source_data <- read_assessment(assessment_id, read_database())
recipe <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     paste0(assessment_id, ".R")), envir = recipe)
translated <- recipe$translate_stock(source_data)
dat <- do.call(tinyAM::make_dat, c(
  list(obs = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

inputs <- source_data$inputs
outputs <- source_data$outputs
stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1957:2025),
  identical(translated$ages, 1:10),
  identical(source_data$assessment$terminal_year, 2025L),
  identical(source_data$assessment$estimate_terminal_year, 2026L),
  nrow(translated$obs$catch) == 690L,
  nrow(translated$obs$index) == 870L,
  identical(as.integer(table(translated$obs$index$survey)),
           c(300L, 88L, 152L, 180L, 150L)),
  setequal(unique(translated$obs$index$survey),
           c("BTS-IBTS Q3", "BTS-Isis", "IBTS Q1", "SNS1", "SNS2")),
  !any(translated$obs$index$survey == "SNS2" &
         translated$obs$index$year == 2003L),
  all(translated$obs$index$samp_time[translated$obs$index$survey == "IBTS Q1"] == .125),
  all(translated$obs$index$samp_time[translated$obs$index$survey != "IBTS Q1"] == .75),
  all(is.finite(translated$obs$weight$M_assumption)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$F_settings$process, "rw"),
  identical(dat$F_settings$mean_ages, 2:6),
  identical(dim(translated$start_par$log_f), dim(par$log_f)),
  identical(dim(translated$start_par$log_n), dim(par$log_n)),
  identical(names(translated$start_par$log_q), names(par$log_q)),
  all(is.finite(translated$start_par$log_f)),
  all(is.finite(translated$start_par$log_n)),
  all(is.finite(translated$start_par$log_q)),
  sum(inputs$type == "M") == 10L,
  sum(inputs$type == "maturity") == 10L,
  sum(outputs$measure == "numbers_at_age") == 690L,
  sum(outputs$measure == "fishing_mortality_at_age") == 690L,
  all(outputs$lwr[outputs$year <= 2025 &
                    outputs$measure %in% c("recruitment", "SSB", "Fbar")] > 0)
)

cat("North Sea plaice translation structure passed.\n")
