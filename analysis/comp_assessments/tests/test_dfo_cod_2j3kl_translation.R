source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "audit_assumptions.R"))

database <- read_database()
source_data <- read_assessment("dfo_cod_2j3kl_2025", database)
audit <- audit_assumptions("dfo_cod_2j3kl_2025", source_data$assumptions)
source_year <- source_data$inputs$year
source_basis <- source_data$inputs$year_basis
translated <- .test_stock(source_data)
obs <- translated$obs
all_maturity <- source_data$inputs[
  source_data$inputs$type == "maturity" &
    source_data$inputs$measure == "maturity_at_age",
  , drop = FALSE
]

tinyAM::check_obs(obs)
readiness <- read.csv(file.path(root, "results", "observation_readiness.csv"),
                      stringsAsFactors = FALSE)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = translated$years, ages = translated$ages),
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
  nrow(all_maturity) == 71L * 15L,
  setequal(all_maturity$year, 1954:2024),
  setequal(all_maturity$age, 0:14),
  all(all_maturity$year_basis == "calendar_year"),
  !anyDuplicated(paste(all_maturity$year, all_maturity$age)),
  !any(source_data$inputs$type == "maturity_cohort" |
         source_data$inputs$measure == "maturity_cohort"),
  nrow(obs$catch) == 57L * 13L,
  nrow(obs$weight) == 57L * 13L,
  nrow(obs$maturity) == 57L * 13L,
  !any(source_data$inputs$type == "maturity" &
         !is.na(source_data$inputs$year_basis) &
         source_data$inputs$year_basis == "birth_cohort"),
  !any(grepl("cohort", source_maturity$notes, ignore.case = TRUE)),
  nrow(obs$index) == 507L + 7L * 12L,
  isTRUE(all.equal(obs$maturity$obs, expected_maturity)),
  all(obs$weight$M_assumption == median(source_data$outputs$value[
    source_data$outputs$measure == "Mbar"])),
  setequal(unique(as.character(obs$index$q_key[obs$index$survey == "DFO fall RV survey"])),
           c("RV age 2", "RV age 3", "RV age 4", "RV age 5", "RV age 6+")),
  identical(source_data$inputs$year, source_year),
  identical(source_data$inputs$year_basis, source_basis),
  readiness$maturity_rows[readiness$assessment_id == "dfo_cod_2j3kl_2025"] == 1065L,
  readiness$maturity_full_year_age_grid[readiness$assessment_id == "dfo_cod_2j3kl_2025"],
  !any(grepl("maturity_cohort", names(readiness))),
  audit$tinyam_support[audit$setting == "maturity_at_age"] == "supported",
  dat$N_settings$process == "off",
  dat$N_settings$init == "exp",
  dat$F_settings$process == "rw",
  dat$M_settings$process == "ar1",
  min(dat$M_settings$years) == 1984,
  all(obs$index$samp_time[obs$index$survey == "DFO fall RV survey"] == 0.75),
  all(obs$index$smith_sound_year[obs$index$survey != "DFO fall RV survey"] == 0),
  !anyNA(dat$q_modmat)
)

smith <- obs$index[obs$index$survey == "Smith Sound acoustic survey", ]
samples <- source_data$inputs[source_data$inputs$survey %in% "Smith Sound acoustic survey" &
                               source_data$inputs$unit %in% "fish_sampled", ]
for (year in unique(smith$year)) {
  cp <- samples[samples$year == year & samples$age %in% 1:13, ]
  weights <- source_data$inputs[source_data$inputs$type == "weight" &
                                 source_data$inputs$year == year, ]
  w <- weights$value[match(cp$age, weights$age)]
  biomass <- source_data$inputs$value[source_data$inputs$survey %in% "Smith Sound acoustic survey" &
                                      source_data$inputs$year == year &
                                      source_data$inputs$measure == "total_biomass"]
  expected <- biomass * 1000 * cp$value / sum(cp$value * w)
  fitted_rows <- smith[smith$year == year, ]
  stopifnot(isTRUE(all.equal(fitted_rows$obs, expected[match(fitted_rows$age, cp$age)])))
}
juvenile_trial <- .test_stock(source_data, juveniles = TRUE)
juvenile_obs <- juvenile_trial$obs$index
stopifnot(identical(juvenile_trial$ages, 0:14),
          all(juvenile_obs$samp_time[juvenile_obs$survey == "Fleming juvenile survey"] == 9.5/12),
          all(juvenile_obs$samp_time[juvenile_obs$survey == "Newman Sound juvenile survey"] == 9/12),
          setequal(unique(as.character(juvenile_obs$q_key[grepl("juvenile", juvenile_obs$survey)])),
                   c("juvenile age 0", "juvenile age 1")),
          identical(source_data$inputs$year, source_year))

for (measure in c("numbers_at_age", "biomass_at_age", "mature_biomass_at_age",
                  "total_mortality_at_age", "natural_mortality_at_age", "fishing_mortality_at_age")) {
  rows <- source_data$outputs[source_data$outputs$measure == measure, ]
  stopifnot(setequal(rows$age, 0:14), !anyDuplicated(paste(rows$year, rows$age)),
            all(is.na(rows$se)), all(is.na(rows$lwr)), all(is.na(rows$upr)))
}

cat("Northern cod translation tests passed.\n")
