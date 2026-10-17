source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("afsc_alaska_plaice_bsai_2024")
s <- .test_stock(x)
stopifnot(is.null(s$fit), identical(s$years,2000:2023), identical(s$ages,5:25),
  !"start_par" %in% names(s$fit_call), ncol(s$dat$q_modmat) == 0,
  all(s$obs$maturity$spawn_time == 0),
  all(s$obs$index$samp_time == .5),
  all(s$obs$pop_covariates$M_assumption == .13),
  identical(as.integer(sort(unique(s$obs$catch$year[!is.na(s$obs$catch$obs)]))),c(2000L,2002L,2003L)))
p <- x$inputs[x$inputs$measure == "proportion_at_age", ]
for (id in unique(p$observation_id)) {
  z <- p[p$observation_id == id, ]
  total <- x$inputs[x$inputs$measure == "total_biomass" &
    x$inputs$type == z$type[1] & x$inputs$year == z$year[1], ]
  w <- s$catch_weights$value[match(s$key(z),s$key(s$catch_weights))]/1000
  expected <- total$value * 1000 / sum(z$value*w)
  if (z$type[1] == "index") expected <- expected / exp(.182322)
  actual <- sum(s$observations$value[s$observations$observation_id == id])
  stopifnot(abs(actual/expected-1) < 1e-5)
}
pooled <- s$comparison_outputs[s$comparison_outputs$measure == "numbers_at_age", ]
stopifnot(nrow(pooled) == 50*26, all(pooled$sex == "combined"), all(is.na(pooled$se)))
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=25,assumptions=x$assumptions)
stopifnot(nrow(ref$pop$recruitment) == 24,all(ref$pop$recruitment$unit == "million fish"))
cat("Alaska plaice pooling, fixed q calibration, fit-disabled preparation and reference alignment passed.\n")
