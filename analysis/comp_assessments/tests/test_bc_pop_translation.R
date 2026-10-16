source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("dfo_pop_bc_2023")
original <- x
s <- .test_stock(x)
stopifnot(identical(x,original),is.null(s$fit),identical(s$years,2003:2022),identical(s$ages,5:60),
  !"start_par" %in% names(s$fit_call),s$settings$rec_settings$sd==.9,
  is.null(s$settings$M_settings$mu_supplied),length(unique(s$obs$index$survey))==3,
  all(abs(s$obs$index$samp_time-5/12)<1e-12),all(s$obs$maturity$spawn_time==0),
  all(is.finite(s$obs$weight$obs)),all(s$obs$maturity$obs>=0 & s$obs$maturity$obs<=1),
  s$comparison_definitions$ssb$status=="non_equivalent",
  all(is.na(s$comparison_outputs$lwr)),all(is.na(s$comparison_outputs$upr)))
for(year in s$years) {
  raw <- s$observations[s$observations$type=="catch" & s$observations$year==year,]
  pooled <- s$obs$catch[s$obs$catch$year==year,]
  if(length(unique(raw$fleet))==3) {
    expected <- tapply(raw$value,raw$age,sum)[as.character(s$ages)]
    stopifnot(max(abs(pooled$obs-expected))<1e-7)
  } else stopifnot(all(is.na(pooled$obs)))
}
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=60,assumptions=x$assumptions,
  comparison_definitions=s$comparison_definitions)
stopifnot(all(ref$pop$M$est==.05229),all(is.na(ref$pop$N$est)))
cat("BC POP female conversion, complete-area pooling, missing-area traps, native timing and disabled-fit preparation passed.\n")
