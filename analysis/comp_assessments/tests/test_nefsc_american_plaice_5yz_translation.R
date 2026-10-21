source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

record <- read_assessment("nefsc_american_plaice_5yz_2022_run29f4")
stock <- .test_stock(record)

stopifnot(
  is.null(stock$fit),
  identical(stock$years, 1980:2019),
  identical(stock$ages, 1:11),
  stock$age_plus_group == 11L,
  !"start_par" %in% names(stock$fit_call),
  all(stock$obs$maturity$spawn_time == 0.25),
  all(stock$obs$pop_covariates$M_assumption == 0.3),
  length(unique(stock$obs$index$q_key)) == 4L * 11L,
  all(is.finite(stock$obs$index$relative_sd)),
  all(stock$obs$index$relative_sd > 0),
  all(stock$obs$index$samp_time %in% c(0.25, 0.75)),
  any(stock$obs$index$samp_time[stock$obs$index$survey == "NEFSC spring Bigelow"] == 0.75),
  any(stock$obs$index$samp_time[stock$obs$index$survey == "NEFSC fall Albatross"] == 0.25),
  identical(stock$settings$N_settings$process, "iid"),
  identical(stock$settings$N_settings$init, "free"),
  identical(stock$settings$rec_settings$process, "iid"),
  identical(stock$settings$F_settings$process, "ar1"),
  identical(stock$settings$M_settings$process, "off"),
  identical(all.vars(stock$settings$index_settings$q_form), "q_key"),
  stock$comparison_definitions$F_bar$status == "matched",
  stock$comparison_definitions$M$status == "approximate",
  stock$comparison_definitions$M_bar$status == "approximate"
)

cat("NEFSC American plaice 5YZ fit-disabled preparation and source mapping passed.\n")
