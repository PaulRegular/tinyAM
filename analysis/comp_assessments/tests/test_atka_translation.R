source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
record <- read_assessment("afsc_atka_mackerel_bsai_2024")
stock <- .test_stock(record)
stopifnot(is.null(stock$fit), identical(stock$years, 1990:2023),
  identical(stock$ages, 2:11),
  all(stock$obs$maturity$spawn_time == 7 / 12),
  all(stock$obs$pop_covariates$M_assumption == .3),
  max(abs(stock$obs$index$samp_time - 6.5 / 12)) < 1e-12,
  all(stock$obs$index$year <= 2022),
  all(stock$obs$index$relative_sd > 0),
  all(stock$obs$index$q_block[stock$obs$index$age >= 10] == "10"),
  !"start_par" %in% names(stock$fit_call))
ref <- database_to_tam_ref("afsc_atka_mackerel_bsai_2024", stock$comparison_outputs,
  obs = stock$obs, years = stock$years, ages = stock$ages, age_plus_group = 11,
  assumptions = record$assumptions)
stopifnot(nrow(ref$pop$recruitment) == length(stock$years),
  all(ref$pop$recruitment$est ==
  record$outputs$value[record$outputs$measure == "numbers_at_age" &
    record$outputs$age == 2 & record$outputs$year %in% stock$years]),
  all(ref$pop$recruitment$unit == "million fish"),
  sum(stock$comparison_outputs$measure == "mature_biomass_at_age") == 48 * 11)
cat("Atka preparation, native units, time, age restriction and recruitment alignment passed.\n")
