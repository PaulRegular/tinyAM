source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("joint_hake_coastal_2025")
original <- x
s <- .test_stock(x)
stopifnot(identical(x,original),is.null(s$fit),identical(s$years,1975:2024),
  identical(s$ages,2:15),!"start_par" %in% names(s$fit_call),
  all(s$obs$maturity$spawn_time==0),all(s$obs$index$samp_time==.5),
  s$settings$rec_settings$sd==1.4,length(s$dat$M_terms)==1L,
  abs(unique(s$dat$log_mu_supplied_m)+1.60944)<1e-12,
  identical(levels(s$obs$index$q_block),as.character(2:6)),
  all(s$obs$index$q_block[s$obs$index$age>=6]=="6"),
  all(abs(tapply(s$inputs$value[s$inputs$measure=="proportion_at_age"],
    s$inputs$observation_id[s$inputs$measure=="proportion_at_age"],sum)-1)<1e-10))
i <- x$inputs
for(measure in c("weight_at_age","maturity_fecundity_at_age")) {
  z <- i[i$measure==measure & i$age>=15 & i$year %in% s$years,]
  stopifnot(all(vapply(split(z$value,paste(z$type,z$survey,z$fleet,z$year)),
    function(v)max(v)-min(v)<1e-10,logical(1))))
}
w <- s$obs$weight
m <- s$obs$maturity
source_product <- i[i$measure=="maturity_fecundity_at_age",]
expected <- source_product$value[match(paste(w$year,w$age),paste(source_product$year,source_product$age))]*.5
stopifnot(max(abs(w$obs*m$obs-expected))<1e-10)
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=15,assumptions=x$assumptions,
  comparison_definitions=s$comparison_definitions)
stopifnot(nrow(ref$pop$recruitment)==50,all(ref$pop$recruitment$unit=="million fish"),
  is.null(ref$pop$F) || all(is.na(ref$pop$F$est)))
defaults <- tinyAM::make_par(s$dat)
partial <- s$retry_start
e <- new.env(parent=environment(partial))
e$source <- x
zero <- e$source$outputs$measure=="numbers_at_age" &
  e$source$outputs$year==1976 & e$source$outputs$age==2
stopifnot(sum(zero)==1)
e$source$outputs$value[zero] <- 0
environment(partial) <- e
start <- partial()
merged <- tinyAM:::.merge_start_par(defaults,start)
missing <- setdiff(names(defaults$log_r),names(start$log_r))
stopifnot(length(missing)>0,all(is.finite(unlist(start))),
  identical(merged$log_r[missing],defaults$log_r[missing]),
  identical(merged$log_f,defaults$log_f),identical(merged$log_q,defaults$log_q),
  all(merged$log_r[names(start$log_r)]==start$log_r))
cat("Coastal hake native-percent conversion, fecundity/sex scaling, q sharing, plus biology and fit-disabled preparation passed.\n")
