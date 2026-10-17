## Observations ----
years <- 1978:2023
ages <- 2:11
inputs <- source$inputs

catch <- inputs[
  inputs$type == "catch" & inputs$measure == "numbers_at_age" &
    inputs$season == "spring", ,
  drop = FALSE
]
gear_catch <- catch[c("year", "age", "fleet", "value")]
catch <- aggregate(value ~ year + age, catch, sum)
catch_rows <- data.frame(
  assessment_id = source$assessment$assessment_id,
  type = "catch", measure = "numbers_at_age", basis = "numbers",
  fleet = "", survey = "", sex = "", region = "", season = "spring",
  year = catch$year, year_basis = "calendar_year", age = catch$age,
  value = catch$value, unit = "thousand fish", sampling_time = NA_real_,
  source_type = "translation_assumption",
  source_reference = "DFO Research Document 2024/058, Tables 4 and 8",
  transformation = "Summed fixed- and mobile-gear catch-at-age for each year and age.",
  notes = "Combined stream for tinyAM; the source assessment models its catch observations differently.",
  stringsAsFactors = FALSE
)

gear_weights <- inputs[
  inputs$type == "catch_weight" & inputs$measure == "weight_at_age" &
    inputs$season == "spring" & inputs$value > 0, ,
  drop = FALSE
]
gear_weights <- merge(
  gear_weights[c("year", "age", "fleet", "value")], gear_catch,
  by = c("year", "age", "fleet"), suffixes = c("_weight", "_catch")
)
combined_weights <- do.call(rbind, lapply(
  split(gear_weights, interaction(gear_weights$year, gear_weights$age)),
  function(x) {
    if (!nrow(x)) {
      return(NULL)
    }
    value <- if (sum(x$value_catch) > 0) {
      stats::weighted.mean(x$value_weight, x$value_catch)
    } else {
      mean(x$value_weight)
    }
    data.frame(year = x$year[1], age = x$age[1], value = value)
  }
))

weight_grid <- expand.grid(year = years, age = ages)
weight_grid <- merge(
  weight_grid, combined_weights,
  by = c("year", "age"), all.x = TRUE, sort = FALSE
)
weight_grid <- weight_grid[order(weight_grid$year, weight_grid$age), ]

for (age in ages) {
  rows <- weight_grid$age == age
  observed <- rows & is.finite(weight_grid$value)
  weight_grid$value[rows] <- stats::approx(
    weight_grid$year[observed], weight_grid$value[observed],
    xout = weight_grid$year[rows], rule = 2
  )$y
}

source_weight <- matrix(
  weight_grid$value,
  nrow = length(years), byrow = TRUE,
  dimnames = list(year = years, age = ages)
)
stock_weight <- source_weight

for (i in 2:length(years)) {
  for (j in 2:length(ages)) {
    stock_weight[i, j] <- sqrt(
      source_weight[i - 1, j - 1] * source_weight[i, j]
    )
  }
}

weight_rows <- expand.grid(year = years, age = ages)
weight_rows <- data.frame(
  assessment_id = source$assessment$assessment_id,
  type = "weight", measure = "weight_at_age", basis = "kg_per_fish",
  fleet = "", survey = "", sex = "", region = "", season = "spring",
  year = weight_rows$year, year_basis = "calendar_year",
  age = weight_rows$age,
  value = as.vector(stock_weight), unit = "kg per fish",
  sampling_time = NA_real_, source_type = "translation_assumption",
  source_reference = "DFO Research Document 2024/058, Tables 5 and 9; DFO Research Document 2022/068, p. 6",
  transformation = paste(
    "Catch-number-weighted arithmetic mean of available fixed/mobile weights",
    "by year-age (unweighted mean if both catches are zero); missing",
    "values interpolated within age. Beginning-year age a uses the geometric",
    "mean of age a-1 in year t-1 and age a in year t; age 2 and 1978 use",
    "same-year weights because earlier ages/years are unavailable."
  ),
  notes = "Fit-only stock-weight approximation; the official fitted weight surface and gear-combination rule are not published.",
  stringsAsFactors = FALSE
)

add_composition_fields <- function(x) {
  x$observation_id <- ""
  x$length_bin <- NA_real_
  x$length_bin_lower <- NA_real_
  x$length_bin_upper <- NA_real_
  x$sample_size <- NA_real_
  x$age_error <- NA_real_
  x$partition <- NA_real_
  x[names(inputs)]
}

catch_rows <- add_composition_fields(catch_rows)
weight_rows <- add_composition_fields(weight_rows)

fit_inputs <- inputs[
  !(inputs$type == "catch" & inputs$measure == "numbers_at_age") &
    !(inputs$type == "weight" & inputs$measure == "weight_at_age"), ,
  drop = FALSE
]
fit_inputs <- rbind(fit_inputs, catch_rows, weight_rows)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  fit_inputs,
  years = years,
  ages = ages,
  sampling_times = c(
    "Spring fixed-gear CPUE" = 0.25,
    "4Tmno acoustic survey" = 0.75
  ),
  surveys = c("Spring fixed-gear CPUE", "4Tmno acoustic survey"),
  assumptions = source$assumptions
)

obs$pop_covariates$M_process_center <- 0.2
obs$pop_covariates$age_blocks <- cut_ages(
  obs$pop_covariates$age,
  c(2, 4, 6, 8, 11)
)

obs$pop_covariates$fishery_period <- cut_years(
  obs$pop_covariates$year,
  c(1978, 1990, 2005, max(obs$pop_covariates$year))
)

# Restrict the catchability random walk to CPUE.
# The first CPUE year defines the baseline; its RW effect is zero.
is_cpue <- obs$index$survey == "Spring fixed-gear CPUE"
first_cpue_year <- min(obs$index$year[is_cpue])

obs$index$cpue_active <- as.numeric(is_cpue)
obs$index$cpue_year <- ifelse(
  is_cpue, obs$index$year, first_cpue_year
)


## Background and comparisons ----
age_plus_group <- 11

comparison_outputs <- source$outputs[
  source$outputs$measure %in% c(
    "numbers_at_age", "fishing_mortality_at_age",
    "biomass_at_age", "recruitment", "Fbar",
    "total_numbers", "total_biomass", "SSB",
    "mature_biomass_at_age"
  ) &
    !(source$outputs$measure %in% c("total_numbers", "total_biomass") &
        !is.na(source$outputs$age_group)), ,
  drop = FALSE
]


# Request the additional age-specific SSB comparison in the shared fish/kg calculation.
comparison_scales <- c(ssb_mat = .001)

comparison_definitions <- list(
  ssb = list(
    status = "matched",
    definition = "January 1 mature biomass at ages 4-11+, derived from reported biomass-at-age and maturity",
    reason = "This common-definition comparison is not the official assessment's April 1 SSB. Source uncertainty is unavailable."
  )
)


background <- c(
  "### Southern Gulf spring-spawning herring: detailed assessment to 2023",
  "",
  print_sources(source$assessment),
  "",
  "The canonical record represents DFO's official 2024 assessment. The 2026 assessment is recorded separately as summary-only because its detailed methods and numerical tables have not been recovered. The 2022 methods report provides model context where the 2024 supporting document is silent.",
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",

  "| Years | Population estimates cover 1978–2023; spring fixed-gear CPUE is available for 1990–2021. | Fit 1978–2023 using the available catch, CPUE and acoustic observations. | The fitted population period extends beyond the CPUE series. |",

  "| Ages | Ages 2–11+, with age 11 as the plus group. | Retain ages 2–11+, including the plus group. | The official population age structure is preserved. |",

  "| N | Age-2 recruitment varies annually; older cohorts follow survival with initial-cohort deviations. | Random-walk recruitment, IID older-cohort survival deviations and exponential initial abundance. | The IID survival process adds flexibility absent from the documented SCA survival model; initial states also differ. |",

  "| F | The documented SCA uses logistic fishery selectivity in three periods, with annual variation in fishing mortality. | Correlated F random walks (cor_rw), allowing annual log-F increments to be correlated across ages. | The translation represents persistent time variation but does not impose period-specific logistic selectivity or a shared annual fishing intensity. |",

  "| M | Log-M follows separate random walks for ages 2–6 and 7–11+, with initial-M prior means of 0.2 and fixed increment SD of 0.075. | Two age-blocked AR1 M processes centered on 0.2, with estimated process variation. | Stationary AR1 dynamics provide a mortality anchor but do not reproduce the source random walks, initial-state priors or fixed increment SD. |",

  "| Catch | Fixed- and mobile-gear catch-at-age are reported separately in thousand fish. The official SCA uses total catch and age-composition likelihoods. | Sum both gears into one age-specific catch stream, convert to fish and estimate one lognormal catch-error SD. | The total-plus-composition observation likelihood is not reproduced. |",

  "| CPUE catchability | The 2022 methods describe annual random-walk variation in log catchability, beginning in 1990, with fixed increment SD of 0.1. | Estimate a CPUE baseline and annual log-q random walk with an estimated increment SD. The initial random-walk deviation is zero. | The temporal structure is similar, but its process SD is estimated because the fixed 0.1 specification did not converge in the translation. |",

  "| Survey selectivity | The documented SCA uses age-dependent survey selectivity, including logistic CPUE selectivity and an acoustic selectivity curve. | Estimate separate non-decreasing age-q curves for CPUE and acoustic observations using mono(age, by = survey). | The flexible monotone curves do not impose the source's logistic shapes or historical CPUE selectivity periods. |",

  "| Indices | Spring fixed-gear CPUE and acoustic observations enter age-composition and aggregate-biomass likelihoods. The documented acoustic biomass likelihood has weight 3. | Fit CPUE ages 4–11 and acoustic ages 2–10 as age-specific lognormal indices, with separate survey SDs and approximate timings of 0.25 and 0.75. | Composition, aggregate-biomass likelihoods and the acoustic likelihood weight are not reproduced. |",

  "| Index scale | The acoustic age-index series is documented in thousand fish in DFO 2016/060, Table 16. | Convert acoustic observations to fish through canonical database units and retain a log-q link. | The documented numerical scale is retained without a manual multiplier. |",

  "| Weights and maturity | Maturity is knife-edge between ages 3 and 4. Beginning-year weights combine gear-specific weights and use adjacent age-year geometric means. | Combine fixed/mobile weights using catch-number weighting, interpolate gaps and apply the adjacent-cell geometric mean and maturity schedule. | Gear weighting is inferred; its median difference from weights implied by published biomass/N is 0.073%. Boundary cells remain approximate. |",

  "| Comparison | January 1 N, biomass-at-age, F-at-age, recruitment and Fbar are reported. | Compare age-specific estimates and derived total abundance, total biomass and January 1 mature biomass. | The SSB comparison is January 1 mature biomass, not official April 1 SSB. Source age-specific M and uncertainty are unavailable. |",

  "",
  "The revised translation uses correlated F random walks rather than stationary AR1 F around age-blocked means. Period-specific logistic fishery selectivity was explored but is not retained; the final model does not use fishery_period in its F specification.",
  "",
  "CPUE catchability has a survey-specific estimated baseline and an annual random walk on log q. Its first-year random-walk effect is fixed at zero, while the baseline determines absolute catchability. The random-walk increment SD is estimated rather than fixed at the 0.1 value documented in the 2022 SCA methods.",
  "",
  "The acoustic survey has a separate baseline and monotone age-q curve but no temporal catchability random walk. Both surveys retain their published observations and separate estimated observation SDs.",
  "",
  "The translation is a structural approximation, not a reproduction of the official likelihood. The main remaining differences concern aggregate-biomass and age-composition observations, logistic selectivity, M initial-state priors and process penalties.",
  "",
  "Acoustic values are converted using canonical units. The January 1 mature-biomass comparison is distinct from official April 1 SSB, and the 2026 assessment's updated biomass data are not mixed into the 2024 record."
)



## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()

  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(
      process = "iid",
      init = "exp"
    ),
    F_settings = list(
      process = "cor_rw",
      mu_form = NULL,
      mean_ages = 6:8
    ),
    M_settings = list(
      process = "ar1",
      mu_form = NULL,
      mu_supplied = ~M_process_center,
      age_breaks = c(2, 7, 11),
      first_dev_year = 1978L
    ),
    catch_settings = list(
      sd_form = ~1
    ),
    index_settings = list(
      q_form = ~ 0 + survey +
        mono(age, by = survey) +
        rw(cpue_year, by = cpue_active),
      q_link = "log",
      sd_form = ~ 0 + survey
    ),
    silent = silent
  )
}
