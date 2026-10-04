root <- "analysis/comp_assessments"
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "audit_assumptions.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "dfo_cod_2j3kl_2025.R"), envir = stock)
database <- read_database()
source_data <- read_assessment("dfo_cod_2j3kl_2025", database)
audit <- audit_assumptions("dfo_cod_2j3kl_2025", source_data$assumptions)
source_year <- source_data$inputs$year
source_basis <- source_data$inputs$year_basis
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
readiness <- read.csv(file.path(root, "results", "fit_readiness.csv"),
                      stringsAsFactors = FALSE)
dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))

source_maturity <- source_data$inputs[
  source_data$inputs$type == "maturity" &
    source_data$inputs$measure == "maturity_at_age" &
    source_data$inputs$year_basis == "calendar_year" &
    source_data$inputs$year %in% translated$years &
    source_data$inputs$age %in% translated$ages,
  , drop = FALSE
]
source_key <- paste(source_maturity$year, source_maturity$age)
expected_maturity <- as.numeric(source_maturity$value)[
  match(paste(obs$maturity$year, obs$maturity$age), source_key)
]

stopifnot(
  identical(translated$years, 1968:2024),
  identical(translated$ages, 2:14),
  nrow(source_maturity) == 57L * 13L,
  all(source_maturity$type == "maturity"),
  all(source_maturity$measure == "maturity_at_age"),
  all(source_maturity$year_basis == "calendar_year"),
  !any(source_data$inputs$type == "maturity_cohort" |
         source_data$inputs$measure == "maturity_cohort"),
  nrow(obs$catch) == 57L * 13L,
  nrow(obs$weight) == 57L * 13L,
  nrow(obs$maturity) == 57L * 13L,
  !any(source_data$inputs$type == "maturity" &
         !is.na(source_data$inputs$year_basis) &
         source_data$inputs$year_basis == "birth_cohort"),
  !any(grepl("cohort", source_maturity$notes, ignore.case = TRUE)),
  nrow(obs$index) == 507L,
  sum(obs$index$obs == 0) == 106L,
  isTRUE(all.equal(obs$maturity$obs, expected_maturity)),
  all(obs$weight$M_assumption == 0.312),
  setequal(levels(obs$index$q_key), c("age_2", "age_3", "age_4", "age_5", "age_6_plus")),
  identical(source_data$inputs$year, source_year),
  identical(source_data$inputs$year_basis, source_basis),
  readiness$maturity_rows[readiness$assessment_id == "dfo_cod_2j3kl_2025"] == 1065L,
  readiness$maturity_full_year_age_grid[readiness$assessment_id == "dfo_cod_2j3kl_2025"],
  !("maturity_cohort_rows" %in% names(readiness)),
  audit$tinyam_support[audit$setting == "maturity_at_age"] == "supported",
  dat$N_settings$process == "iid",
  dat$N_settings$init == "free",
  dat$F_settings$process == "ar1",
  dat$M_settings$process == "off",
  all(obs$index$samp_time == 0.75)
)

cat("Northern cod translation tests passed.\n")
