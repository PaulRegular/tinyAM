root <- "analysis/comp_assessments"
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path(root, "R", "database_to_tam_ref.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "dfo_herring_4tvn_spring_2024.R"), envir = stock)
source_data <- read_assessment("dfo_herring_4tvn_spring_2024", read_database())
original <- source_data
translated <- stock$translate_stock(source_data)
obs <- translated$obs

stopifnot(
  identical(translated$years, 1978:2023),
  identical(translated$ages, 2:11),
  tinyAM::check_obs(obs),
  nrow(obs$catch) == 460L,
  nrow(obs$index) == 526L,
  nrow(obs$weight) == 460L,
  nrow(obs$maturity) == 460L,
  setequal(unique(obs$index$survey),
           c("Spring fixed-gear CPUE", "4Tmno acoustic survey")),
  all(obs$index$samp_time[obs$index$survey == "Spring fixed-gear CPUE"] == 0.25),
  all(obs$index$samp_time[obs$index$survey == "4Tmno acoustic survey"] == 0.75),
  all(obs$maturity$obs[obs$maturity$age %in% 2:3] == 0),
  all(obs$maturity$obs[obs$maturity$age >= 4] == 1),
  all(obs$weight$M_process_center == 0.2),
  !"M_assumption" %in% names(obs$weight),
  identical(attr(obs, "translation")$M$status, "estimated_in_source"),
  identical(source_data, original),
  !any(translated$comparison_outputs$measure == "natural_mortality_at_age"),
  all(c("total_numbers", "total_biomass", "SSB", "mature_biomass_at_age") %in%
        translated$comparison_outputs$measure),
  all(is.na(translated$comparison_outputs$age_group[
    translated$comparison_outputs$measure %in% c("total_numbers", "total_biomass")
  ]))
)

acoustic_source <- source_data$inputs[
  source_data$inputs$type == "index" &
    source_data$inputs$measure == "numbers_at_age" &
    source_data$inputs$survey == "4Tmno acoustic survey", , drop = FALSE
]
acoustic_observed <- obs$index[
  obs$index$survey == "4Tmno acoustic survey", , drop = FALSE
]
acoustic_match <- merge(
  acoustic_observed[c("year", "age", "obs")],
  acoustic_source[c("year", "age", "value")], by = c("year", "age")
)
stopifnot(
  nrow(acoustic_source) == 270L,
  all(range(as.numeric(acoustic_source$year)) == c(1994, 2023)),
  all(range(as.numeric(acoustic_source$age)) == c(2, 10)),
  all(is.finite(acoustic_source$value)),
  all(acoustic_source$unit == "number (index scale; multiplier not stated)"),
  all(acoustic_source$basis == "index_scale"),
  all(acoustic_source$source_type == "official_table"),
  all(grepl("DFO Research Document 2024/058 Table 15",
            acoustic_source$source_reference, fixed = TRUE)),
  acoustic_source$value[acoustic_source$year == 1994 & acoustic_source$age == 3] == 231932,
  acoustic_source$value[acoustic_source$year == 2023 & acoustic_source$age == 10] == 140,
  nrow(acoustic_observed) == 270L,
  nrow(acoustic_match) == 270L,
  isTRUE(all.equal(acoustic_match$obs, acoustic_match$value)),
  !anyNA(obs$index$cpue_period),
  all(as.character(acoustic_observed$cpue_period) == "baseline"),
  identical(translated$settings$index_settings$q_link, "log")
)

latest <- read_assessment("dfo_herring_4tvn_spring_2026", read_database())
latest_biomass <- latest$inputs[
  latest$inputs$type == "index" &
    latest$inputs$measure == "total_biomass" &
    latest$inputs$survey == "4Tmno acoustic survey", , drop = FALSE
]
stopifnot(
  nrow(latest_biomass) == 32L,
  all(range(as.numeric(latest_biomass$year)) == c(1994, 2025)),
  all(latest_biomass$unit == "tonnes"),
  all(latest_biomass$source_type == "official_machine_readable"),
  all(grepl("herring_historic_biomass_2025.csv",
            latest_biomass$source_reference, fixed = TRUE)),
  latest_biomass$value[latest_biomass$year == 2023] == 7358,
  latest_biomass$value[latest_biomass$year == 2025] == 25607,
  !any(source_data$inputs$measure == "total_biomass" &
         source_data$inputs$survey == "4Tmno acoustic survey")
)

source_catch <- source_data$inputs[
  source_data$inputs$type == "catch" &
    source_data$inputs$measure == "numbers_at_age" &
    source_data$inputs$season == "spring", , drop = FALSE
]
source_catch <- aggregate(value ~ year + age, source_catch, sum)
observed_catch <- obs$catch[c("year", "age", "obs")]
catch_comparison <- merge(observed_catch, source_catch, by = c("year", "age"))
stopifnot(
  nrow(catch_comparison) == nrow(source_catch),
  isTRUE(all.equal(catch_comparison$obs, catch_comparison$value * 1000))
)

dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)
stopifnot(
  as.character(dat$M_settings$age_blocks["2"]) == "2-6",
  as.character(dat$M_settings$age_blocks["7"]) == "7-11",
  identical(dim(par$log_f), c(46L, 10L)),
  identical(dim(par$log_m), c(46L, 2L)),
  sum(dat$fill_missing_map) == 0L,
  isTRUE(all.equal(rownames(par$log_m), as.character(translated$years))),
  all(is.finite(translated$start_par$log_r)),
  all(is.finite(translated$start_par$log_f)),
  all(is.finite(translated$start_par$log_n)),
  identical(names(translated$start_par), names(par)),
  all(is.finite(translated$start_par$log_q)),
  identical(translated$settings$N_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$M_settings$process, "ar1"),
  identical(translated$settings$M_settings$first_dev_year, 1978L)
)

reference <- database_to_tam_ref(
  "dfo_herring_4tvn_spring_2024", translated$comparison_outputs,
  obs = obs, years = translated$years, ages = translated$ages,
  terminal_year = 2023, age_plus_group = 11,
  comparison_scales = translated$comparison_scales
)
stopifnot(
  is.null(reference$pop$M),
  !"M" %in% names(reference$pop),
  nrow(reference$pop$abundance) == 46L,
  nrow(reference$pop$biomass) == 46L,
  nrow(reference$pop$ssb) == 46L,
  all(reference$pop$ssb$source_type == "derived_source_output"),
  all(is.na(reference$pop$ssb$se)),
  all(is.na(reference$pop$ssb$lwr)),
  all(is.na(reference$pop$ssb$upr)),
  all(grepl("January 1", reference$pop$ssb$notes, fixed = TRUE)),
  grepl("not the accepted assessment's April 1 SSB",
        translated$comparison_definitions$ssb$reason, fixed = TRUE)
)

annual_numbers <- aggregate(est ~ year, reference$pop$N, sum)
annual_mature <- aggregate(est ~ year, reference$pop$ssb_mat, sum)
stopifnot(
  isTRUE(all.equal(reference$pop$abundance$est,
                   annual_numbers$est[match(reference$pop$abundance$year, annual_numbers$year)])),
  isTRUE(all.equal(reference$pop$ssb$est,
                   annual_mature$est[match(reference$pop$ssb$year, annual_mature$year)]))
)

matched_fit <- list(
  dat = list(years = translated$years, ages = translated$ages),
  pop = list(ssb_mat = transform(reference$pop$ssb_mat, est = est * 1000))
)
differences <- .assessment_percent_differences(matched_fit, reference)
matched <- differences[differences$metric == "ssb_mat", ]
stopifnot(nrow(matched) == 460L, all(matched$comparison_status == "matched"),
          all(matched$unit == "t"),
          all(matched$source == matched$tinyAM),
          all(matched$percent_difference[!is.na(matched$percent_difference)] == 0))

# Verify the gear mixture separately from the adjacent age/year geometric mean.
gear_weights <- source_data$inputs[source_data$inputs$type == "catch_weight" &
                                    source_data$inputs$value > 0, ]
gear_catch <- source_data$inputs[source_data$inputs$type == "catch" &
                                  source_data$inputs$measure == "numbers_at_age", ]
cell_weight <- function(year, age) {
  w <- gear_weights[gear_weights$year == year & gear_weights$age == age, ]
  n <- gear_catch$value[match(paste(w$year, w$age, w$fleet),
                             paste(gear_catch$year, gear_catch$age, gear_catch$fleet))]
  stats::weighted.mean(w$value, n)
}
expected_weight <- sqrt(cell_weight(1978, 4) * cell_weight(1979, 5))
stopifnot(abs(obs$weight$obs[obs$weight$year == 1979 & obs$weight$age == 5] -
               expected_weight) < 1e-12)
implicit <- merge(reference$pop$N[c("year", "age", "est")],
                  reference$pop$biomass_at_age[c("year", "age", "est")],
                  by = c("year", "age"), suffixes = c("_N", "_B"))
implicit <- merge(implicit, obs$weight[c("year", "age", "obs")],
                  by = c("year", "age"))
stopifnot(stats::median(abs(100 * (implicit$obs - implicit$est_B / implicit$est_N) /
                             (implicit$est_B / implicit$est_N))) < 0.1)

cat("Southern Gulf spring herring translation tests passed.\n")
