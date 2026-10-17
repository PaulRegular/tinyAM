source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("afsc_yellowfin_sole_ebs_2024")
original <- x
s <- .test_stock(x)
stopifnot(identical(x,original),is.null(s$fit),identical(s$years,1982:2023),
  identical(s$ages,2:20),!"start_par" %in% names(s$fit_call),
  all(s$obs$maturity$spawn_time==4/12),all(s$obs$index$samp_time==.5),
  all(s$obs$pop_covariates$M_assumption==.12),
  all(is.finite(as.matrix(s$obs$index[c("temperature","start_date","interaction","log_sd")]))),
  all(s$observations$sex=="female"),all(s$comparison_outputs$sex=="female"),
  !any(s$comparison_outputs$measure=="natural_mortality_at_age"))
print_error <- function(x) .5*10^(floor(log10(pmax(abs(x),1e-300)))-5)
for(year in 1982:2023) {
  z <- s$comparison_outputs[s$comparison_outputs$measure=="mature_biomass_at_age" &
    s$comparison_outputs$year==year,]
  native <- x$outputs[x$outputs$measure=="SSB" & x$outputs$year==year,]
  nn <- x$outputs[x$outputs$measure=="numbers_at_age" & x$outputs$sex=="female" & x$outputs$year==year,]
  ff <- x$outputs[x$outputs$measure=="fishing_mortality_at_age" & x$outputs$sex=="female" & x$outputs$year==year,]
  coef <- z$value/nn$value
  bound <- 1000*print_error(native$value) +
    sum(coef*((nn$value+print_error(nn$value))*exp(print_error(ff$value)*4/12)-nn$value))
  stopifnot(abs(sum(z$value)-native$value*1000)<=bound+1e-4)
}
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=20,assumptions=x$assumptions,
  comparison_definitions=s$comparison_definitions)
stopifnot(nrow(ref$pop$recruitment)==42,all(ref$pop$recruitment$unit=="billion fish"),
  all(ref$pop$N$sex=="female"))
cat("Yellowfin sole female conversion, timing, fixed M, SSB identity and fit-disabled preparation passed.\n")
