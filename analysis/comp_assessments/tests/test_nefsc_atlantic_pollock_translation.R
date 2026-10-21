source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/read_database.R")
database <- read_database()
source_data <- read_assessment("nefsc_pollock_5yz_2024", database)
translated <- .test_stock(source_data)

years <- translated$years
ages <- translated$ages
inputs <- source_data$inputs
obs <- translated$obs
dat <- translated$dat

landings <- inputs[inputs$measure == "landings_numbers_at_age", , drop = FALSE]
discards <- inputs[inputs$measure == "discard_numbers_at_age", , drop = FALSE]
combined <- stats::aggregate(value ~ year + age, rbind(landings, discards), sum)
catch_index <- match(paste(combined$year, combined$age), paste(obs$catch$year, obs$catch$age))
spring <- "NEFSC Bigelow spring trawl survey"
spring_year <- 2009
age_values <- inputs[inputs$type == "index" & inputs$measure == "biomass_at_age" &
  inputs$survey == spring & inputs$year == spring_year, , drop = FALSE]
aggregate_index <- inputs[inputs$type == "index" & inputs$measure == "relative_biomass_index" &
  inputs$survey == spring & inputs$year == spring_year, , drop = FALSE]
survey_weights <- inputs[inputs$type == "weight" & inputs$measure == "weight_at_age" &
  !is.na(inputs$survey) & inputs$survey == spring &
  inputs$year == spring_year, , drop = FALSE]
age_one <- which(obs$index$survey == spring & obs$index$year == spring_year & obs$index$age == 1)
expected_index <- aggregate_index$value *
  (age_values$value[age_values$age == 1] / sum(age_values$value)) /
  survey_weights$value[survey_weights$age == 1]

stopifnot(
  tinyAM::check_obs(obs),
  identical(years, 1970:2023),
  identical(ages, 1:9),
  translated$age_plus_group == 9L,
  nrow(obs$catch) == 54L * 9L,
  nrow(obs$weight) == 54L * 9L,
  nrow(obs$maturity) == 54L * 9L,
  nrow(obs$index) == 105L * 9L,
  isTRUE(all.equal(obs$catch$obs[catch_index], combined$value)),
  setequal(unique(obs$index$survey), c(
    "NEFSC Albatross spring trawl survey",
    "NEFSC Bigelow spring trawl survey",
    "NEFSC Albatross fall trawl survey",
    "NEFSC Bigelow fall trawl survey"
  )),
  setequal(unique(obs$index$samp_time), c(.25, .75)),
  isTRUE(all.equal(obs$index$obs[age_one], expected_index / 1000)),
  all(is.finite(obs$index$relative_sd) & obs$index$relative_sd > 0),
  all(obs$pop_covariates$M_assumption == .2),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$F_settings$process, "ar1"),
  identical(dat$F_settings$mean_ages, 5:7),
  identical(dat$M_settings$process, "off"),
  all(is.finite(tinyAM::make_par(dat)$log_f)),
  is.null(translated$fit_call$start_par),
  grepl("composition sampling uncertainty is not propagated", paste(translated$background, collapse = " "))
)

cat("Atlantic pollock tinyAM observations, transformations, settings and default initialization passed.\n")
