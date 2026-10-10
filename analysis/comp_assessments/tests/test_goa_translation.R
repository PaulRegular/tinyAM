source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

database <- read_database()
source_data <- read_assessment("afsc_pollock_goa_2024", database)
translated <- .test_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
stopifnot(any(is.finite(obs$catch$obs[obs$catch$age == 1])),
          any(is.finite(obs$catch$obs[obs$catch$age == 2])),
          all(is.na(obs$catch$obs[obs$catch$year < 1975])))
shelikof <- obs$index$survey == "Shelikof winter acoustic"
stopifnot(!any(shelikof & obs$index$age == 3),
          any(is.finite(obs$index$obs[shelikof & obs$index$age >= 4])))

dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
stopifnot(identical(dat$years, translated$years),
          identical(dat$ages, translated$ages))

stopifnot(identical(deparse(translated$settings$catch_settings$sd_form), "~1"),
          qr(dat$q_modmat)$rank == ncol(dat$q_modmat),
          !any(grepl("^adfg_q_", names(obs$index))),
          nrow(translated$catch_reporting$weights) == 55 * 10,
          nrow(translated$catch_reporting$totals) == 55)

# Retain paired-age q blocks, with deterministic survival and bounded q.
stopifnot(identical(translated$settings$N_settings$process, "off"),
          identical(translated$settings$N_settings$init, "exp"),
          identical(translated$settings$F_settings$process, "rw"),
          identical(translated$settings$index_settings$q_link, "logit"),
          identical(deparse(translated$settings$index_settings$q_form),
                    "~0 + q_key + environmental_effect + adfg_year"),
          ncol(dat$q_mono_modmat) == 0L)
par <- tinyAM::make_par(dat)
stopifnot("logit_q" %in% names(par), !"log_n" %in% names(par))
reported <- RTMB::MakeADFun(function(p) tinyAM::nll_fun(p, dat), par,
                           silent = TRUE)$report()
stopifnot(all(is.finite(reported$N)),
          all(exp(reported$log_q_obs) == 0.5))
for (survey in unique(obs$index$survey)) {
  z <- obs$index[obs$index$survey == survey, ]
  expected <- c("1-2", "3-4", "5-6", "7-8", "9-10")[ceiling(z$age / 2)]
  stopifnot(identical(as.character(z$q_age_block), expected),
            all(vapply(split(z$q_key, z$q_age_block, drop = TRUE),
                       function(x) length(unique(x)) == 1L, logical(1))))
}

# Audit the full reconstruction before the pooled Shelikof bin is omitted.
translation_inputs <- source_data$inputs[
  !(source_data$inputs$type == "catch" &
      source_data$inputs$measure == "proportion_at_age"), ]
full <- database_to_tam_obs(
  "afsc_pollock_goa_2024", translation_inputs,
  years = translated$years, ages = translated$ages,
  weight_survey = "", index_weight_source = "source",
  maturity_multiplier = 0.5, assumptions = source_data$assumptions
)
weights <- source_data$inputs[source_data$inputs$type == "weight" &
                                source_data$inputs$measure == "weight_at_age" &
                                !is.na(source_data$inputs$survey), ]
for (survey in unique(full$index$survey)) {
  z <- full$index[full$index$survey == survey, ]
  w <- weights[weights$survey == survey, ]
  index <- match(paste(z$year, z$age), paste(w$year, w$age))
  stopifnot(!anyNA(index), all(w$unit == "kg"),
            all(w$basis == "kg_per_fish"))
  reconstructed <- tapply(z$obs * w$value[index], z$year, sum)
  totals <- source_data$inputs[source_data$inputs$type == "index" &
                                source_data$inputs$measure == "total_biomass" &
                                !is.na(source_data$inputs$survey) &
                                source_data$inputs$survey == survey, ]
  stopifnot(all(totals$unit == "million t"))
  expected <- totals$value[match(as.integer(names(reconstructed)), totals$year)] * 1e9
  stopifnot(all(abs(reconstructed / expected - 1) < 1e-12))
  retained <- z[!(survey == "Shelikof winter acoustic" & z$age == 3), ]
  actual <- obs$index[obs$index$survey == survey, ]
  stopifnot(isTRUE(all.equal(actual$obs, retained$obs,
                            check.attributes = FALSE)))
}
# Independent rounded SAFE Table 1.11 check (millions of fish, not billions).
rows <- full$index$survey == "Shelikof winter acoustic" &
          full$index$year == 1992 & full$index$age %in% 4:9
stopifnot(all(abs(full$index$obs[rows] / 1e6 -
                   c(188.1, 368.0, 84.1, 85.0, 171.2, 32.7)) < 0.05))

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
