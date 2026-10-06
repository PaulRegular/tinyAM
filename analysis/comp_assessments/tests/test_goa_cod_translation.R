root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

source_data <- read_assessment("afsc_cod_goa_2026", read_database())
inputs_before <- source_data$inputs
stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "afsc_cod_goa_2026.R"), stock)
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
stopifnot(
  identical(source_data$inputs, inputs_before),
  identical(translated$years, 2007:2025),
  identical(translated$ages, 1:10),
  nrow(obs$catch) == 19L * 10L,
  nrow(obs$index) == 10L * 10L,
  all(is.finite(obs$catch$obs)),
  all(is.finite(obs$index$obs)),
  identical(translated$settings$F_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$M_settings$process, "off"),
  all(is.finite(obs$weight$M_assumption)),
  all(obs$weight$M_assumption[obs$weight$year %in% 2014:2016] == 0.84),
  all(obs$weight$M_assumption[!obs$weight$year %in% 2014:2016] == 0.50),
  all(translated$comparison_outputs$unit == "kg"),
  identical(sort(unique(translated$comparison_outputs$year)), 2007:2025)
)

index_totals <- source_data$inputs[
  source_data$inputs$type == "index" &
    source_data$inputs$survey == "NMFS bottom-trawl survey" &
    source_data$inputs$measure == "total_numbers" &
    source_data$inputs$year %in% translated$years, , drop = FALSE
]
index_sum <- aggregate(obs ~ year, obs$index, sum)
expected_index <- index_totals$value * 1000
stopifnot(isTRUE(all.equal(
  index_sum$obs[match(index_totals$year, index_sum$year)],
  expected_index, tolerance = 1e-8
)))

dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
stopifnot(identical(dat$years, translated$years),
          identical(dat$ages, translated$ages))

cat("GOA Pacific cod translation tests passed.\n")
