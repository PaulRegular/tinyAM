source("analysis/comp_assessments/tests/helper_stock.R")
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_record <- read_assessment("dfo_herring_4tvn_fall_2024", database)
source_original <- source_record
translated <- .test_stock(source_record)
obs <- translated$obs

fixed_catch <- source_record$inputs[
  source_record$inputs$type == "catch" &
    source_record$inputs$measure == "numbers_at_age" &
    source_record$inputs$season == "fall" &
    source_record$inputs$fleet == "Fixed gear", ,
  drop = FALSE
]
expected_catch <- aggregate(value ~ year + age, fixed_catch, sum)
expected_catch_index <- match(
  paste(obs$catch$year, obs$catch$age),
  paste(expected_catch$year, expected_catch$age)
)
cpue_index <- obs$index$survey %in% c(
  "Fall fixed-gear CPUE - North",
  "Fall fixed-gear CPUE - Middle",
  "Fall fixed-gear CPUE - South"
)

stopifnot(
  identical(translated$years, 1978:2023),
  identical(translated$ages, 2:11),
  check_obs(obs),
  nrow(obs$catch) == length(1978:2023) * length(2:11),
  all(is.finite(obs$catch$obs)),
  !anyNA(expected_catch_index),
  isTRUE(all.equal(obs$catch$obs, expected_catch$value[expected_catch_index])),
  all(is.finite(obs$weight$obs) & obs$weight$obs > 0),
  all(obs$maturity$obs[obs$maturity$age < 4] == 0),
  all(obs$maturity$obs[obs$maturity$age >= 4] == 1),
  all(obs$maturity$spawn_time == 212 / 365),
  identical(sort(unique(obs$index$survey)), sort(c(
    "Fall fixed-gear CPUE - North",
    "Fall fixed-gear CPUE - Middle",
    "Fall fixed-gear CPUE - South",
    "4Tmno acoustic survey",
    "DFO multispecies RV survey"
  ))),
  all(obs$index$samp_time[grepl("Fall fixed-gear CPUE|4Tmno", obs$index$survey)] == 0.8),
  all(obs$index$samp_time[obs$index$survey == "DFO multispecies RV survey"] == 0.75),
  identical(range(obs$index$age[cpue_index]), c(4, 11)),
  identical(range(obs$index$age[obs$index$survey == "4Tmno acoustic survey"]), c(2, 10)),
  identical(range(obs$index$age[obs$index$survey == "DFO multispecies RV survey"]), c(2, 11)),
  all(obs$index$cpue_active[grepl("Fall fixed-gear CPUE", obs$index$survey)] == 1),
  all(obs$index$cpue_active[!grepl("Fall fixed-gear CPUE", obs$index$survey)] == 0),
  isTRUE(all.equal(
    obs$catch$obs[obs$catch$year == 1978 & obs$catch$age == 4],
    sum(source_record$inputs$value[
      source_record$inputs$type == "catch" &
        source_record$inputs$measure == "numbers_at_age" &
        source_record$inputs$season == "fall" &
        source_record$inputs$fleet == "Fixed gear" &
        source_record$inputs$year == 1978 &
        source_record$inputs$age == 4
    ])
  )),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$F_settings$process, "cor_rw"),
  identical(translated$settings$M_settings$process, "ar1"),
  !"start_par" %in% names(translated$fit_call),
  all(c("numbers_at_age", "fishing_mortality_at_age", "SSB") %in%
    unique(translated$comparison_outputs$measure)),
  "ssb" %in% translated$comparison_aggregates,
  all(is.na(translated$comparison_outputs$age_group[
    translated$comparison_outputs$measure == "SSB"
  ])),
  translated$comparison_definitions$biomass_at_age$status == "approximate",
  translated$comparison_definitions$ssb$status == "matched",
  any(grepl("duplicat", translated$background, ignore.case = TRUE)),
  identical(source_record, source_original)
)

cat("Southern Gulf fall herring translation tests passed.\n")
