source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("afsc_northern_rock_sole_ebs_2024")
original <- x
s <- .test_stock(x)
stopifnot(identical(x,original),is.null(s$fit),identical(s$years,1982:2023),
  identical(s$ages,5:20),!"start_par" %in% names(s$fit_call),
  all(s$obs$maturity$spawn_time==1/12),all(s$obs$index$samp_time==.5),
  all(s$observations$sex=="female"),all(s$comparison_outputs$sex=="female"),
  all(is.finite(s$obs$index$log_sd)),all(s$obs$index$log_sd>0),
  is.null(s$settings$M_settings$mu_supplied),
  identical(as.character(s$settings$M_settings$mu_form),c("~","1")))
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=20,assumptions=x$assumptions,
  comparison_definitions=s$comparison_definitions)
stopifnot(nrow(ref$pop$recruitment)==42,all(ref$pop$recruitment$unit=="million fish"),
  all(ref$pop$N$sex=="female"),all(ref$pop$M$est==.192),all(ref$pop$M$se==.003))
cat("Northern rock sole female observations, estimated M, units, timing and fit-disabled preparation passed.\n")
