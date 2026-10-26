source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
s <- read_assessment("ices_sandeel_area_4_2025")
before <- s
x <- .test_stock(s)
obs <- x$obs; dat <- x$dat
stopifnot(identical(s, before), is.null(x$fit), tinyAM::check_obs(obs),
  identical(x$years, 1993:2024), identical(x$ages, 1:4),
  nrow(obs$catch) == 128L, nrow(obs$index) == 42L,
  all(obs$catch$unit == "thousand fish"), all(obs$index$age == 1),
  all(obs$maturity$spawn_time == 0),
  dat$N_settings$process == "off", dat$N_settings$init == "free",
  dat$rec_settings$process == "iid", dat$F_settings$process == "rw",
  dat$M_settings$process == "off", ncol(dat$q_modmat) == 4L,
  ncol(dat$sd_index_modmat) == 4L, ncol(dat$sd_catch_modmat) == 2L,
  is.finite(tinyAM::nll_fun(tinyAM::make_par(dat), dat)))
juvenile <- obs$index$source_age == 0
stopifnot(sum(juvenile) == 21L, sum(!juvenile) == 21L,
  all(obs$index$year[juvenile] == obs$index$source_year[juvenile] + 1L),
  all(obs$index$year[!juvenile] == obs$index$source_year[!juvenile]),
  all(obs$index$samp_time[juvenile] == 0), all(obs$index$samp_time[!juvenile] == .9375),
  !any(obs$index$source_year[juvenile] == 2024),
  !any(obs$index$survey == "Dredge (source age 1)" & obs$index$source_year == 2003))
direct <- s$inputs[s$inputs$type == "index", ]
key <- function(survey, year, age) paste(survey, year, age)
source_survey <- sub(" \\(source age [01]\\)$", "", as.character(obs$index$survey))
stopifnot(identical(obs$index$obs, direct$value[match(
  key(source_survey, obs$index$source_year, obs$index$source_age),
  key(direct$survey, direct$year, direct$age))]))

# The unchanged juvenile observation has the same conditional log likelihood.
M0 <- unique(s$inputs$value[s$inputs$type == "M" & s$inputs$age == 0])
R0 <- c(100, 1000, 10000); q0 <- .7; sd <- .3
average_survival <- exp(-.75 * M0) * (-expm1(-.25 * M0)) / (.25 * M0)
N1 <- R0 * exp(-M0)
q_next_january <- q0 * exp(M0) * average_survival
original <- q0 * R0 * average_survival
mapped <- q_next_january * N1
stopifnot(isTRUE(all.equal(original, mapped)),
  isTRUE(all.equal(sum(dnorm(log(original * 1.1), log(original), sd, log = TRUE)),
    sum(dnorm(log(original * 1.1), log(mapped), sd, log = TRUE)))))

inputs <- s$inputs
m <- inputs[inputs$type == "M" & inputs$age >= 1, ]
annual_m <- aggregate(value ~ year + age, m, sum)
key_pop <- function(z) paste(z$year, z$age)
stopifnot(isTRUE(all.equal(obs$pop_covariates$M_assumption,
  annual_m$value[match(key_pop(obs$pop_covariates), key_pop(annual_m))])))
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" & inputs$age >= 1, ]
annual_catch <- aggregate(value ~ year + age, catch, sum)
stopifnot(isTRUE(all.equal(obs$catch$obs,
  annual_catch$value[match(key_pop(obs$catch), key_pop(annual_catch))] * 1000)),
  all(obs$catch$obs[obs$catch$year == 2024] == 0),
  setequal(setdiff(x$years, x$catch_reporting$totals$year), c(2006:2008, 2010:2011, 2013)),
  x$catch_reporting$totals$yield[x$catch_reporting$totals$year == 2024] == 0,
  all(is.finite(x$catch_reporting$weights$weight)), all(x$catch_reporting$weights$weight > 0))
ref <- database_to_tam_ref(s$assessment$assessment_id, x$comparison_outputs, obs = obs,
  years = x$years, ages = x$ages, age_plus_group = 4, assumptions = s$assumptions,
  comparison_aggregates = x$comparison_aggregates, comparison_definitions = x$comparison_definitions)
n <- s$outputs[s$outputs$measure == "numbers_at_age" & s$outputs$age == 1, ]
stopifnot(all(ref$pop$recruitment$est == n$value[match(ref$pop$recruitment$year, n$year)]),
  all(is.na(ref$pop$recruitment$se)), all(attr(ref, "native_pop")$recruitment$age == 0),
  all(ref$pop$F$est[ref$pop$F$year %in% c(2009, 2024)] == 0))
cat("Area 4 annual preparation, unchanged juvenile values/likelihood, provenance, missing catches and same-year reference passed.\n")
