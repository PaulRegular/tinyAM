source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
x <- read_assessment("afsc_cod_ebs_2024")
original <- x
s <- .test_stock(x)
stopifnot(identical(x,original),is.null(s$fit),identical(s$years,2000:2023),
  identical(s$ages,1:12),!"start_par" %in% names(s$fit_call),
  all(s$obs$index$samp_time==.5),all(s$obs$maturity$spawn_time==0),
  all(s$obs$pop_covariates$M_assumption==.3866),
  s$settings$rec_settings$sd==.6908,all(is.finite(s$obs$weight$obs)),
  all(s$obs$weight$obs>0),all(s$obs$maturity$obs>=0 & s$obs$maturity$obs<=.5),
  all(abs(rowSums(s$key)-1)<1e-12))
raw <- x$inputs
total <- raw[raw$measure=="total_numbers",]
active <- raw[raw$type=="index" & raw$measure=="proportion_at_age" & grepl("^Active",raw$notes),]
for(year in setdiff(2000:2023,2020)) {
  p <- active[active$year==year,]; p$value <- p$value/sum(p$value)
  expected <- p$value[p$age %in% s$ages]*total$value[total$year==year]*1000
  got <- s$survey$value[s$survey$year==year & s$survey$age %in% s$ages]
  stopifnot(max(abs(got-expected))<1e-7)
}
stopifnot(all(is.na(s$obs$index$obs[s$obs$index$year==2020])))
length_rows <- raw[raw$type=="catch" & raw$measure=="proportion_at_length" & raw$year %in% s$years,]
coverage <- vapply(split(length_rows,length_rows$year),function(z)
  sum(z$value[z$length_bin %in% s$lengths])/sum(z$value),numeric(1))
stopifnot(min(coverage)>.994,all(coverage<1),
  all(s$comparison_definitions$F$status=="approximate"),
  s$comparison_definitions$ssb$status=="approximate")
stopifnot(s$dat$F_settings$process == "rw")
# The reviewed shared-level alternative must not add free projection means.
candidate_obs <- s$obs
candidate_obs$pop_covariates$F_year <- factor(candidate_obs$pop_covariates$year)
candidate_settings <- s$settings
candidate_settings$F_settings <- list(process="iid",mu_form=~F_year+factor(age))
candidate <- do.call(tinyAM::prepare_tam,c(
  list(data=candidate_obs,years=s$years,ages=s$ages),candidate_settings))
stopifnot(ncol(candidate$F_modmat) == length(s$years) + length(s$ages) - 1L,
  qr(candidate$F_modmat)$rank == ncol(candidate$F_modmat))
projected <- do.call(tinyAM::prepare_tam,c(
  list(data=candidate_obs,years=s$years,ages=s$ages,proj_settings=list(n_proj=2)),candidate_settings))
stopifnot(ncol(projected$F_modmat) == ncol(candidate$F_modmat))
for(age in s$ages) {
  rows <- projected$pop_data$age == age & projected$pop_data$year >= max(s$years)
  stopifnot(all(apply(projected$F_modmat[rows,,drop=FALSE],2,function(x)length(unique(x))) == 1L))
}
ref <- database_to_tam_ref(x$assessment$assessment_id,s$comparison_outputs,
  obs=s$obs,years=s$years,ages=s$ages,age_plus_group=12,assumptions=x$assumptions,
  comparison_definitions=s$comparison_definitions)
stopifnot(all(is.finite(ref$pop$N$est)),all(is.finite(ref$pop$F$est)))
cat("EBS cod empirical age allocation, omitted length mass, number-index scaling, fixed M, proxy biology and disabled-fit preparation passed.\n")
