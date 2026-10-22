pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

assessment_id <- "ices_sole_eastern_channel_7d_2026"
database <- read_database()
source_record <- read_assessment(assessment_id, database)
original_inputs <- source_record$inputs
original_outputs <- source_record$outputs
stock <- .assessment_stock(source_record, do_fit = FALSE, silent = TRUE)
prepared <- stock$env

stopifnot(
  is.null(stock$error),
  is.null(prepared$fit),
  isTRUE(tinyAM::check_obs(prepared$obs)),
  identical(prepared$years, 1982:2025),
  identical(prepared$ages, 1:11),
  identical(prepared$age_plus_group, 11L),
  nrow(prepared$obs$catch) == 44L * 11L,
  sum(prepared$obs$catch$obs == 0, na.rm = TRUE) == 1L,
  nrow(prepared$obs$index) == 228L,
  setequal(unique(prepared$obs$index$survey),
    c("UK(E&W)-BTS-Q3", "UK(E&W)-YFS", "FR-YFS")),
  sum(prepared$obs$index$survey == "UK(E&W)-BTS-Q3") == 171L,
  sum(prepared$obs$index$survey == "UK(E&W)-YFS") == 20L,
  sum(prepared$obs$index$survey == "FR-YFS") == 37L,
  !any(prepared$obs$index$survey == "UK(E&W)-BTS-Q3" &
    prepared$obs$index$year == 2022L),
  !any(prepared$obs$index$survey == "UK(E&W)-BTS-Q3" &
    prepared$obs$index$year >= 2010L & prepared$obs$index$age < 4L),
  all(prepared$obs$index$samp_time == 0.625),
  length(unique(prepared$obs$index$q_key[
    prepared$obs$index$survey == "UK(E&W)-BTS-Q3" &
      prepared$obs$index$age %in% 5:6
  ])) == 1L,
  length(unique(prepared$obs$index$q_key[
    prepared$obs$index$survey == "UK(E&W)-BTS-Q3" &
      prepared$obs$index$age %in% 1:4
  ])) == 4L,
  length(unique(prepared$obs$index$sd_block[
    prepared$obs$index$survey == "UK(E&W)-BTS-Q3" &
      prepared$obs$index$age %in% 3:6
  ])) == 1L,
  all(prepared$obs$pop_covariates$M_assumption == 0.1),
  identical(prepared$N_settings$process, "iid"),
  identical(prepared$N_settings$init, "exp"),
  identical(prepared$rec_settings$process, "rw"),
  identical(prepared$F_settings$process, "rw"),
  identical(prepared$F_settings$mean_ages, 3:7),
  identical(prepared$M_settings$process, "off"),
  identical(deparse(prepared$catch_settings$sd_form), deparse(~0 + sd_block)),
  identical(deparse(prepared$index_settings$q_form), deparse(~0 + q_key)),
  identical(deparse(prepared$index_settings$sd_form), deparse(~0 + sd_block)),
  identical(source_record$inputs, original_inputs),
  identical(source_record$outputs, original_outputs),
  identical(source_record$inputs, read_assessment(assessment_id, database)$inputs),
  all(prepared$comparison_outputs$source_type[
    prepared$comparison_outputs$measure == "mature_biomass_at_age"
  ] == "derived_source_output")
)

ssb_age <- prepared$comparison_outputs[
  prepared$comparison_outputs$measure == "mature_biomass_at_age", , drop = FALSE
]
ssb_total <- source_record$outputs[
  source_record$outputs$measure == "SSB" &
    source_record$outputs$year <= 2025, c("year", "value"), drop = FALSE
]
ssb_sum <- stats::aggregate(value ~ year, ssb_age, sum)
ssb_check <- merge(ssb_total, ssb_sum, by = "year", suffixes = c("_reported", "_derived"))
stopifnot(
  nrow(ssb_age) == 44L * 11L,
  all(is.na(ssb_age$se)),
  max(abs(ssb_check$value_reported - ssb_check$value_derived) /
    ssb_check$value_reported) < 0.0002,
  any(abs(ssb_check$value_reported[ssb_check$year >= 2021L] -
    ssb_check$value_derived[ssb_check$year >= 2021L]) > 0)
)

cat("Eastern Channel sole tinyAM preparation and common-definition SSB checks passed.\n")
