root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "afsc_pollock_goa_2024.R"), envir = stock)
database <- read_database()
source_data <- read_assessment("afsc_pollock_goa_2024", database)
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
stopifnot(any(is.finite(obs$catch$obs[obs$catch$age == 1])),
          any(is.finite(obs$catch$obs[obs$catch$age == 2])),
          all(is.na(obs$catch$obs[obs$catch$year < 1975])))
shelikof <- obs$index$survey == "Shelikof winter acoustic"
stopifnot(!any(shelikof & obs$index$age == 3),
          any(is.finite(obs$index$obs[shelikof & obs$index$age >= 4])))

dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
stopifnot(identical(dat$years, translated$years),
          identical(dat$ages, translated$ages))

stopifnot(identical(deparse(translated$settings$catch_settings$sd_form), "~1"),
          qr(dat$q_modmat)$rank == ncol(dat$q_modmat),
          !any(grepl("^adfg_q_", names(obs$index))),
          nrow(translated$catch_reporting$weights) == 55 * 10,
          nrow(translated$catch_reporting$totals) == 55)

# Both trawls rise with age; both acoustics decline, with independent steps.
q <- drop(dat$q_mono_modmat %*% rep(0.1, ncol(dat$q_mono_modmat)))
for (survey in unique(dat$obs$index$survey)) {
  rows <- dat$obs$index$survey == survey
  by_age <- tapply(q[rows], dat$obs$index$age[rows], mean)
  direction <- if (survey %in% c("NMFS bottom trawl", "ADF&G crab/groundfish trawl")) 1 else -1
  stopifnot(all(direction * diff(by_age) > 0))
}
stopifnot(length(unique(dat$q_mono_steps$by_level)) == 4)

catch_numbers <- source_data$inputs[source_data$inputs$type == "catch" &
                                      source_data$inputs$measure == "numbers_at_age", ]
stopifnot(nrow(catch_numbers) == 49 * 15,
          setequal(catch_numbers$year, 1975:2023),
          setequal(catch_numbers$age, 1:15),
          all(catch_numbers$source_type == "official_table"),
          catch_numbers$value[catch_numbers$year == 2023 & catch_numbers$age == 1] == 0.43,
          catch_numbers$value[catch_numbers$year == 2023 & catch_numbers$age == 2] == 8.57,
          obs$catch$obs[obs$catch$year == 2023 & obs$catch$age == 1] == 430000,
          obs$catch$obs[obs$catch$year == 2023 & obs$catch$age == 2] == 8570000,
          obs$catch$obs[obs$catch$year == 2023 & obs$catch$age == 10] == 8890000)

totals <- source_data$inputs[source_data$inputs$type == "catch" &
                              source_data$inputs$measure == "total_biomass", ]
stopifnot(identical(translated$catch_reporting$totals$yield, totals$value * 1000))
cat("GOA pollock translation tests passed.\n")
