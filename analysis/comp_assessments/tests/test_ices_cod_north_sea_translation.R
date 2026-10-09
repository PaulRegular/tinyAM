root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "ices_cod_north_sea_2025.R"), envir = stock)
database <- read_committed_database()
source_data <- read_committed_assessment("ices_cod_north_sea_2025", database)
source_inputs <- source_data$inputs
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = translated$years, ages = translated$ages),
  translated$settings
))

weights <- source_inputs[source_inputs$type == "weight" &
                           source_inputs$measure == "weight_at_age" &
                           source_inputs$year %in% translated$years &
                           source_inputs$age %in% translated$ages, ]
weight_average <- stats::aggregate(as.numeric(value) ~ year + age, weights, mean)
maturity <- source_inputs[source_inputs$type == "maturity" &
                            source_inputs$measure == "maturity_at_age" &
                            source_inputs$year %in% translated$years &
                            source_inputs$age %in% translated$ages, ]
maturity_average <- stats::aggregate(as.numeric(value) ~ year + age, maturity, mean)
comparison <- translated$comparison_outputs
source_n <- source_data$outputs[source_data$outputs$measure == "numbers_at_age", ]
sum_n <- stats::aggregate(value ~ year + age, source_n, sum)
source_ssb <- source_data$outputs[source_data$outputs$measure == "SSB", ]
sum_ssb <- stats::aggregate(value ~ year, source_ssb, sum)
source_f <- source_data$outputs[
  source_data$outputs$measure == "fishing_mortality_at_age", ]
source_f_n <- source_n[match(
  paste(source_f$year, source_f$age, source_f$region),
  paste(source_n$year, source_n$age, source_n$region)
), ]
source_f_weighted <- data.frame(
  year = source_f$year,
  age = source_f$age,
  value = source_f$value * source_f_n$value
)
sum_f_weighted <- stats::aggregate(value ~ year + age, source_f_weighted, sum)
sum_f_weights <- stats::aggregate(value ~ year + age, source_f_n, sum)
mean_f <- sum_f_weighted$value / sum_f_weights$value
source_catch <- source_inputs[source_inputs$type == "catch" &
                                source_inputs$measure == "numbers_at_age" &
                                source_inputs$year %in% translated$years &
                                source_inputs$age %in% translated$ages, ]

stopifnot(
  identical(translated$years, 1983:2022),
  identical(translated$ages, 1:7),
  nrow(obs$catch) == 40L * 7L,
  nrow(obs$weight) == 40L * 7L,
  nrow(obs$maturity) == 40L * 7L,
  isTRUE(all.equal(
    obs$catch$obs,
    as.numeric(source_catch$value)[match(
      paste(obs$catch$year, obs$catch$age),
      paste(source_catch$year, source_catch$age)
    )] * 1000
  )),
  setequal(unique(obs$index$survey), "Survey_Q34"),
  all(obs$index$samp_time == 0.75),
  isTRUE(all.equal(obs$weight$obs,
    weight_average$`as.numeric(value)`[match(paste(obs$weight$year, obs$weight$age),
                                             paste(weight_average$year, weight_average$age))])),
  isTRUE(all.equal(obs$maturity$obs,
    maturity_average$`as.numeric(value)`[match(paste(obs$maturity$year, obs$maturity$age),
                                                paste(maturity_average$year, maturity_average$age))])),
  identical(source_data$inputs, source_inputs),
  all(is.na(comparison$region)),
  isTRUE(all.equal(
    comparison$value[comparison$measure == "numbers_at_age"],
    sum_n$value[match(
      paste(comparison$year[comparison$measure == "numbers_at_age"],
            comparison$age[comparison$measure == "numbers_at_age"]),
      paste(sum_n$year, sum_n$age)
    )]
  )),
  isTRUE(all.equal(
    comparison$value[comparison$measure == "SSB"],
    sum_ssb$value[match(comparison$year[comparison$measure == "SSB"],
                        sum_ssb$year)]
  )),
  isTRUE(all.equal(
    comparison$value[comparison$measure == "fishing_mortality_at_age"],
    mean_f[match(
      paste(comparison$year[comparison$measure == "fishing_mortality_at_age"],
            comparison$age[comparison$measure == "fishing_mortality_at_age"]),
      paste(sum_f_weighted$year, sum_f_weighted$age)
    )]
  )),
  dat$N_settings$process == "rw",
  dat$F_settings$process == "rw",
  dat$M_settings$process == "off"
)

cat("North Sea cod translation tests passed.\n")
