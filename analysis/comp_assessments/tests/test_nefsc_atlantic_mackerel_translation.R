root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "nefsc_atlantic_mackerel_2018"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
recipe <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     paste0(assessment_id, ".R")), envir = recipe)
translated <- recipe$translate_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

inputs <- source_data$inputs
outputs <- source_data$outputs
summary_outputs <- database$outputs[database$outputs$assessment_id == "nefsc_atlantic_mackerel_2025", ]
egg_ssb <- database$inputs[database$inputs$assessment_id == assessment_id & database$inputs$measure == "total_biomass", ]
expected_catch <- inputs$value[
  inputs$type == "catch" & inputs$measure == "numbers_at_age" &
    inputs$year == 1968 & inputs$age == 1
] * 1000
expected_index <- inputs$value[
  inputs$type == "index" & inputs$measure == "numbers_at_age" &
    inputs$survey == "NEFSC Albatross spring trawl" &
    inputs$year == 1974 & inputs$age == 3
]
expected_n <- outputs$value[
  outputs$type == "population" & outputs$measure == "numbers_at_age" &
    outputs$year == 2016 & outputs$age == 1
]
expected_recruitment <- outputs$value[
  outputs$type == "recruitment" & outputs$measure == "recruitment" &
    outputs$year == 2016
]

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1968:2016),
  identical(translated$ages, 1:10),
  translated$age_plus_group == 10L,
  nrow(translated$obs$catch) == 49L * 10L,
  nrow(translated$obs$weight) == 49L * 10L,
  nrow(translated$obs$maturity) == 49L * 10L,
  nrow(translated$obs$index) == 320L,
  setequal(unique(translated$obs$index$survey),
           c("NEFSC Albatross spring trawl", "NEFSC Bigelow spring trawl")),
  setequal(unique(translated$obs$index$age), 3:10),
  all(translated$obs$index$age[translated$obs$index$survey ==
                                 "NEFSC Bigelow spring trawl"] <= 7),
  isTRUE(all.equal(translated$obs$catch$obs[
    translated$obs$catch$year == 1968 & translated$obs$catch$age == 1
  ], expected_catch)),
  isTRUE(all.equal(translated$obs$index$obs[
    translated$obs$index$survey == "NEFSC Albatross spring trawl" &
      translated$obs$index$year == 1974 & translated$obs$index$age == 3
  ], expected_index)),
  all(translated$obs$weight$M_assumption == 0.2),
  all(abs(translated$obs$index$samp_time - 95.7 / 365) < 1e-12),
  identical(dat$N_settings$process, "iid"),
  identical(dat$F_settings$process, "ar1"),
  identical(dat$F_settings$mean_ages, 6:10),
  identical(dat$M_settings$process, "off"),
  all(is.finite(par$log_n)),
  nrow(source_data$inputs[source_data$inputs$measure == "total_biomass", ]) == 17L,
  nrow(outputs[outputs$measure == "natural_mortality_at_age", ]) == 49L * 10L,
  isTRUE(all.equal(expected_n, expected_recruitment)),
  inputs$value[inputs$type == "catch" & inputs$year == 1968 & inputs$age == 1] == 161471,
  inputs$value[inputs$type == "weight" & inputs$year == 1968 & inputs$age == 1] == 0.148,
  inputs$value[inputs$type == "maturity" & inputs$year == 1968 & inputs$age == 10] == 0.999,
  inputs$value[inputs$type == "index" & inputs$survey == "NEFSC Albatross spring trawl" & inputs$year == 1974 & inputs$age == 3] == 1.09,
  expected_n == 455.43,
  outputs$value[outputs$measure == "SSB" & outputs$year == 2016] == 43519,
  outputs$value[outputs$measure == "Fbar" & outputs$year == 2016] == 0.468,
  outputs$age[outputs$measure == "recruitment" & outputs$year == 2016] == 1,
  egg_ssb$value[egg_ssb$year == 1979] == 1131094,
  nrow(summary_outputs) == 30L,
  summary_outputs$value[summary_outputs$measure == "SSB" & summary_outputs$year == 2024] == 94702,
  summary_outputs$lwr[summary_outputs$measure == "SSB" & summary_outputs$year == 2024] == 52539,
  summary_outputs$upr[summary_outputs$measure == "SSB" & summary_outputs$year == 2024] == 170702,
  summary_outputs$value[summary_outputs$measure == "Fbar" & summary_outputs$year == 2024] == 0.04,
  summary_outputs$value[summary_outputs$measure == "recruitment" & summary_outputs$year == 2024] == 1292885,
  grepl("omitted from the fit", paste(translated$background, collapse = " "))
)

cat("Atlantic mackerel translation structure passed.
")
