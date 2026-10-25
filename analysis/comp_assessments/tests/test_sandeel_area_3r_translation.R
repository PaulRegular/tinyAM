source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source_record <- read_assessment("ices_sandeel_area_3r_2026")
before <- source_record
x <- .test_stock(source_record)
obs <- x$obs; dat <- x$dat
stopifnot(identical(source_record, before), is.null(x$fit),
  identical(x$years, 1987:2025), identical(x$ages, 1:4), tinyAM::check_obs(obs),
  nrow(obs$catch) == 156L, nrow(obs$index) == 88L,
  all(obs$catch$unit == "thousand fish"),
  obs$catch$obs[obs$catch$year == 1987 & obs$catch$age == 1] == (33760 + 65) * 1000,
  all(obs$catch$obs[obs$catch$year == 2025] == 0),
  sum(obs$index$survey == "Dredge") == 20L,
  sum(obs$index$survey == "Acoustic") == 68L,
  all(obs$index$samp_time[obs$index$survey == "Dredge"] == .9375),
  all(obs$index$samp_time[obs$index$survey == "Acoustic"] == .3),
  all(obs$maturity$spawn_time == 0),
  dat$N_settings$process == "off", dat$N_settings$init == "free",
  dat$rec_settings$process == "iid", dat$F_settings$process == "cor_rw",
  is.null(dat$F_settings$mu_form), dat$M_settings$process == "off",
  ncol(dat$q_modmat) == 4L, ncol(dat$sd_index_modmat) == 2L,
  ncol(dat$sd_catch_modmat) == 2L,
  is.finite(tinyAM::nll_fun(tinyAM::make_par(dat), dat)))
acoustic <- obs$index[obs$index$survey == "Acoustic", ]
stopifnot(setequal(acoustic$age, 1:4),
  length(unique(acoustic$q_group[acoustic$age >= 3])) == 1L,
  length(unique(acoustic$q_group)) == 3L,
  length(unique(acoustic$sd_group)) == 1L,
  !anyDuplicated(paste(acoustic$year, acoustic$age)))

inputs <- source_record$inputs
key <- function(z) paste(z$year, z$age)
m <- inputs[inputs$type == "M" & inputs$age %in% x$ages, ]
annual_m <- aggregate(value ~ year + age, m, sum)
stopifnot(isTRUE(all.equal(obs$pop_covariates$M_assumption,
  annual_m$value[match(key(obs$pop_covariates), key(annual_m))])),
  obs$pop_covariates$M_assumption[obs$pop_covariates$year == 1987 &
    obs$pop_covariates$age == 1] == 1.08,
  abs(obs$pop_covariates$M_assumption[obs$pop_covariates$year == 2025 &
    obs$pop_covariates$age == 1] - 1.52) < 1e-12)
stock_w <- inputs[inputs$type == "weight" & inputs$season == "1", ]
stopifnot(isTRUE(all.equal(obs$weight$obs,
  stock_w$value[match(key(obs$weight), key(stock_w))] / 1000)))
direct <- inputs[inputs$type == "index" & inputs$age >= 1, ]
stopifnot(identical(obs$index$obs, direct$value[match(
  paste(obs$index$survey, key(obs$index)), paste(direct$survey, key(direct)))]))

c <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" & inputs$age >= 1, ]
w <- inputs[inputs$type == "catch_weight" & inputs$age >= 1, ]
i <- match(paste(key(c), c$season), paste(key(w), w$season))
c$biomass <- c$value * 1e6 * w$value[i] / 1000
annual <- aggregate(biomass ~ year, c, sum)
cw <- x$catch_reporting$weights
stopifnot(isTRUE(all.equal(x$catch_reporting$totals$yield, annual$biomass)),
  all(x$catch_reporting$totals$yield[x$catch_reporting$totals$year == 2025] == 0),
  abs(cw$weight[cw$year == 2025 & cw$age == 1] - (7.49 + 8.92) / 2000) < 1e-12,
  all(is.finite(cw$weight)), all(cw$weight > 0),
  all(x$comparison_outputs$season[x$comparison_outputs$measure == "fishing_mortality_at_age"] == "annual"))

ref <- database_to_tam_ref(source_record$assessment$assessment_id,
  x$comparison_outputs, obs = obs, years = x$years, ages = x$ages,
  age_plus_group = 4, assumptions = source_record$assumptions,
  comparison_aggregates = x$comparison_aggregates,
  comparison_definitions = x$comparison_definitions)
n <- source_record$outputs[source_record$outputs$measure == "numbers_at_age" &
  source_record$outputs$age == 1, ]
stopifnot(all(ref$pop$recruitment$est == n$value[match(ref$pop$recruitment$year, n$year)]),
  all(is.na(ref$pop$recruitment$se)), all(ref$pop$N$age >= 1),
  all(attr(ref, "native_pop")$recruitment$age == 0),
  all(ref$pop$F$est[ref$pop$F$year == 2025] == 0),
  any(grepl("Finite log-F states cannot represent structural zero F", x$background, fixed = TRUE)))
cat("Area 3r annual preparation, unchanged survey values, q/SD sharing, supplied M, zero catch and same-year reference checks passed.\n")
