source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source_record <- read_assessment("ices_sandeel_area_2r_2026")
before <- source_record
x <- .test_stock(source_record)
obs <- x$obs
dat <- x$dat
stopifnot(identical(source_record, before), is.null(x$fit),
  identical(x$years, 1983:2025), identical(x$ages, 1:4),
  tinyAM::check_obs(obs), nrow(obs$catch) == 172L, nrow(obs$index) == 22L,
  all(obs$catch$unit == "thousand fish"),
  obs$catch$obs[obs$catch$year == 1983 & obs$catch$age == 2] == (6191 + 877) * 1000,
  all(obs$index$age >= 1), sum(obs$index$survey == "Dredge") == 22L,
  all(obs$index$samp_time[obs$index$survey == "Dredge"] == .9375),
  all(obs$maturity$spawn_time == 0),
  dat$N_settings$process == "off", dat$N_settings$init == "free",
  dat$rec_settings$process == "iid", dat$F_settings$process == "iid",
  dat$M_settings$process == "off", ncol(dat$F_modmat) == 13L,
  ncol(dat$q_modmat) == 1L, ncol(dat$sd_index_modmat) == 1L,
  ncol(dat$sd_catch_modmat) == 2L,
  is.finite(tinyAM::nll_fun(tinyAM::make_par(dat), dat)))
expected <- c("1983" = "early.1", "1998" = "early.1", "1999" = "middle1.1",
  "2004" = "middle1.1", "2005" = "middle2.1", "2009" = "middle2.1",
  "2010" = "recent.1", "2025" = "recent.1")
young <- obs$pop_covariates[obs$pop_covariates$age == 1, ]
stopifnot(identical(as.character(young$F_mean_group[
  match(as.integer(names(expected)), young$year)]), unname(expected)))

inputs <- source_record$inputs
key <- function(z) paste(z$year, z$age)
m <- inputs[inputs$type == "M" & inputs$age %in% x$ages, ]
annual_m <- aggregate(value ~ year + age, m, sum)
stopifnot(isTRUE(all.equal(obs$pop_covariates$M_assumption,
  annual_m$value[match(key(obs$pop_covariates), key(annual_m))])),
  obs$pop_covariates$M_assumption[obs$pop_covariates$year == 1983 &
    obs$pop_covariates$age == 1] == .98)
stock_w <- inputs[inputs$type == "weight" & inputs$season == "1", ]
stopifnot(isTRUE(all.equal(obs$weight$obs,
  stock_w$value[match(key(obs$weight), key(stock_w))] / 1000)))
direct <- inputs[inputs$type == "index" & inputs$age >= 1, ]
stopifnot(identical(obs$index$obs, direct$value[match(
  paste(obs$index$survey, key(obs$index)), paste(direct$survey, key(direct)))]))

# Catch reporting retains observed seasonal weights, not fitted-state weights.
c <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" & inputs$age >= 1, ]
w <- inputs[inputs$type == "catch_weight" & inputs$age >= 1, ]
i <- match(paste(key(c), c$season), paste(key(w), w$season))
c$biomass <- c$value * 1e6 * w$value[i] / 1000
annual <- aggregate(biomass ~ year, c, sum)
stopifnot(isTRUE(all.equal(x$catch_reporting$totals$yield, annual$biomass)),
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
  all(attr(ref, "native_pop")$recruitment$age == 0))
cat("Sandeel annual catch/M aggregation, January biology, original indices, settings and same-year recruitment passed.\n")
source("analysis/comp_assessments/R/plot_stock_order.R")
metadata <- assessment_plot_stocks(source_record$assessment, source_record$stocks,
  source_record$assessment$assessment_id)
stopifnot(nrow(metadata) == 1L, metadata$scientific_name == "Ammodytes spp.",
  as.character(metadata$basin) == "Northeast Atlantic")
