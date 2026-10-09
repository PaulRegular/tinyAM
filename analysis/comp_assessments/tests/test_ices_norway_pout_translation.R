pkgload::load_all(".", quiet = TRUE)
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))

database <- read_database()
assessment_id <- "ices_norway_pout_north_sea_2026_benchmark"
source_data <- read_assessment(assessment_id, database)
recipe <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     paste0(assessment_id, ".R")), envir = recipe)
translated <- recipe$translate_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))

inputs <- source_data$inputs
outputs <- source_data$outputs
inputs$year <- suppressWarnings(as.integer(inputs$year))
inputs$age <- suppressWarnings(as.integer(inputs$age))
inputs$value <- suppressWarnings(as.numeric(inputs$value))
outputs$year <- as.integer(outputs$year)
outputs$age <- as.integer(outputs$age)
outputs$value <- as.numeric(outputs$value)

raw_catch <- inputs[
  inputs$type == "catch" & inputs$measure == "numbers_at_age" &
    inputs$year %in% translated$years, ]
raw_catch$fish <- raw_catch$value * 1e6
expected_catch <- stats::aggregate(
  fish ~ year + age, raw_catch, sum
)
actual_catch <- translated$obs$catch[c("year", "age", "obs")]
catch_match <- match(paste(expected_catch$year, expected_catch$age),
                     paste(actual_catch$year, actual_catch$age))

quarterly_m <- inputs[
  inputs$type == "M" & inputs$measure == "natural_mortality_at_age" &
    inputs$year %in% translated$years, ]
quarterly_m$quarter <- as.integer(sub("^Q", "", quarterly_m$season))
expected_m <- do.call(rbind, lapply(split(
  quarterly_m, paste(quarterly_m$year, quarterly_m$age, sep = "_")
), function(x) {
  quarter <- if (x$age[[1]] == 0L) 3:4 else 1:4
  data.frame(year = x$year[[1]], age = x$age[[1]],
             M = sum(x$value[x$quarter %in% quarter]) * 0.25)
}))
actual_m <- translated$obs$weight[c("year", "age", "M_assumption")]
m_match <- match(paste(expected_m$year, expected_m$age),
                 paste(actual_m$year, actual_m$age))

source_weight <- inputs[inputs$type == "weight" &
                          inputs$measure == "weight_at_age" &
                          is.na(inputs$year), ]
source_weight$quarter <- as.integer(sub("^Q", "", source_weight$season))
source_weight <- source_weight[source_weight$quarter ==
                                  ifelse(source_weight$age == 0L, 3L, 1L), ]
expected_weight <- source_weight$value[match(
  translated$obs$weight$age, source_weight$age
)]

source_n <- outputs[outputs$measure == "numbers_at_age", ]
source_n_q1 <- source_n[source_n$season == "Q1" &
                          source_n$year %in% translated$years, ]
source_n_q3_age0 <- source_n[source_n$season == "Q3" &
                               source_n$age == 0L &
                               source_n$year %in% translated$years, ]
source_f <- outputs[outputs$measure == "fishing_mortality_at_age" &
                      outputs$year %in% translated$years, ]
source_f$quarter <- as.integer(sub("^Q", "", source_f$season))
expected_f <- do.call(rbind, lapply(split(
  source_f, paste(source_f$year, source_f$age, sep = "_")
), function(x) {
  quarter <- if (x$age[[1]] == 0L) 3:4 else 1:4
  data.frame(year = x$year[[1]], age = x$age[[1]],
             F = sum(x$value[x$quarter %in% quarter]) * 0.25)
}))
expected_f <- expected_f[order(expected_f$year, expected_f$age), ]

comparison <- translated$comparison_outputs
ssb_q1 <- outputs[outputs$measure == "SSB" & outputs$season == "Q1" &
                    outputs$year %in% translated$years, ]
translated_ssb <- comparison[comparison$measure == "SSB", ]
common_ssb <- comparison[
  comparison$measure == "mature_biomass_at_age", ]
common_ssb <- stats::aggregate(value ~ year, common_ssb, sum)

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1984:2024),
  identical(translated$ages, 0:3),
  nrow(translated$obs$catch) == 41L * 4L,
  nrow(translated$obs$index) == 387L,
  all(translated$obs$index$samp_time %in% c(0.125, 0.625)),
  !any(translated$obs$index$survey == "IBTS other countries Q3" &
         translated$obs$index$age > 1L),
  nlevels(translated$obs$index$q_key) == 12L,
  all(translated$obs$index$q_key[
    translated$obs$index$survey == "EGFS Q3" &
      translated$obs$index$age %in% 2:3
  ] == "egfs_q3_age2_3"),
  all(translated$obs$index$q_key[
    translated$obs$index$survey == "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)" &
      translated$obs$index$age %in% 2:3
  ] == "sgfs_q3_age2_3"),
  length(unique(translated$obs$index$sd_block[
    translated$obs$index$age == 0L &
      translated$obs$index$survey %in% c(
        "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)",
        "SGFS Q3 age 0 (2013 onward)"
      )
  ])) == 1L,
  isTRUE(all.equal(actual_catch$obs[catch_match], expected_catch$fish)),
  isTRUE(all.equal(actual_m$M_assumption[m_match], expected_m$M)),
  isTRUE(all.equal(translated$obs$weight$obs, expected_weight)),
  isTRUE(all.equal(translated$start_par$log_r0,
                   log(as.numeric(source_n_q3_age0$value[1]) * 1e6))),
  max(abs(exp(translated$start_par$log_f) -
            as.matrix(xtabs(F ~ year + age, expected_f)))) < 1e-12,
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$F_settings$process, "ar1"),
  identical(dat$F_settings$mean_ages, 1:2),
  identical(dat$M_settings$process, "off"),
  identical(dim(translated$start_par$log_f), c(41L, 4L)),
  identical(dim(translated$start_par$log_f), dim(tinyAM::make_par(dat)$log_f)),
  !any(comparison$measure == "numbers_at_age" & comparison$age == 0L),
  !any(comparison$year == 2025L),
  identical(translated_ssb$year, ssb_q1$year),
  isTRUE(all.equal(as.numeric(translated_ssb$value), as.numeric(ssb_q1$value))),
  isTRUE(all.equal(common_ssb$value, as.numeric(ssb_q1$value))),
  nrow(comparison[comparison$measure == "biomass_at_age", ]) == 41L * 4L,
  all(is.na(comparison$se[comparison$source_type == "derived_common_definition"]))
)

cat("North Sea Norway pout translation structure passed.\n")
