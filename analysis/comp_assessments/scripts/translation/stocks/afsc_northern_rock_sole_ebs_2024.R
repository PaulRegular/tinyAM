source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1982:2023
ages <- 5:20
inputs <- source$inputs
compositions <- inputs[inputs$measure=="proportion_at_age" & inputs$year %in% years, ]
observations <- lapply(split(compositions,compositions$observation_id),function(p) {
  total <- inputs[inputs$type==p$type[1] & inputs$measure=="total_biomass" & inputs$year==p$year[1], ]
  if(!nrow(total)) return(NULL)
  weight <- inputs[inputs$type==if(p$type[1]=="catch")"catch_weight" else "weight", ]
  weight <- weight[weight$year==p$year[1] &
    if(p$type[1]=="catch") !is.na(weight$fleet) else !is.na(weight$survey), ]
  w <- weight$value[match(paste(p$sex,p$age),paste(weight$sex,weight$age))]
  stopifnot(nrow(total)==1,all(is.finite(w)),sum(p$value*w)>0)
  p$value <- total$value*1e6*p$value/sum(p$value*w)
  p <- p[p$sex=="female", ]
  p$measure <- "numbers_at_age"
  p$basis <- "numbers"
  p$unit <- "fish"
  p$sampling_time <- if(p$type[1]=="index") .5 else NA_real_
  p$notes <- "Translation only: female numbers from total biomass and joint-sex compositions with original sex-specific kg weights; both sexes enter the denominator before retaining females."
  p
})
observations <- do.call(rbind,observations)
biology <- inputs[inputs$sex=="female" & inputs$type %in% c("weight","maturity"), ]
obs <- database_to_tam_obs(source$assessment$assessment_id,rbind(observations,biology),
  years=years,ages=ages,assumptions=source$assumptions)
obs$maturity$spawn_time <- 1/12
sd <- inputs[inputs$measure=="index_sd", ]
biomass <- inputs[inputs$type=="index" & inputs$measure=="total_biomass", ]
cv <- sd$value/biomass$value[match(sd$year,biomass$year)]
obs$index$log_sd <- sqrt(log1p(cv^2))[match(obs$index$year,sd$year)]

## Background and comparisons ----
age_plus_group <- 20
comparison_outputs <- source$outputs[source$outputs$sex=="female", ]
ssb_age <- comparison_outputs[comparison_outputs$measure=="numbers_at_age", ]
w <- inputs[inputs$type=="weight" & inputs$sex=="female" & is.na(inputs$survey), ]
mat <- inputs[inputs$type=="maturity", ]
f <- comparison_outputs[comparison_outputs$measure=="fishing_mortality_at_age", ]
key <- function(x) paste(x$year,x$age)
ssb_age$value <- ssb_age$value*1000*w$value[match(key(ssb_age),key(w))] *
  mat$value[match(key(ssb_age),key(mat))] * exp(-(.192+f$value[match(key(ssb_age),key(f))])/12)
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Female February 1 mature biomass from rounded official N, derived F, rounded estimated female M and original biology; million fish times kg times 1000 gives tonnes. No uncertainty reconstructed."
ssb_age$se <- ssb_age$lwr <- ssb_age$upr <- NA_real_
comparison_outputs <- rbind(comparison_outputs,ssb_age)
comparison_definitions <- list(
  N=list(status="matched",definition="Female begin-year abundance at age",reason="Male states remain canonical, outside this translation."),
  F=list(status="approximate",definition="Female F at age from reported curve parameters",reason="Source F reconstructed from rounded parameters and deviations, with no joint covariance."),
  F_bar=list(status="approximate",definition="N-weighted mean of reconstructed female F",reason="Source F comes from rounded curve parameters and deviations, not an exact native F surface."),
  abundance=list(status="matched",definition="Total female abundance at modeled ages",reason="Not combined-sex stock abundance."),
  recruitment=list(status="matched",definition="Female age-5 entry in the same calendar year",reason="Comparison uses official female N at age 5, not combined-sex age-1 recruitment."),
  ssb=list(status="approximate",definition="Female February 1 SSB at modeled ages",reason="Calculated from rounded official states and parameters without covariance."))
background <- c(
  "### Eastern Bering Sea northern rock sole: official 2024 Model 24.2", "",print_sources(source$assessment),"",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1975–2024, ages 1–20+, both sexes. | 1982–2023, female ages 5–20+. | Use years with survey biomass and age samples, and ages with positive recorded catch weights. Totals are converted using both sexes and all sample ages before restriction. |",
  "| N and recruitment | Deterministic sex-specific survival, free initial ages, Ricker relationship with lag one over 1978–2018 and fixed log SD 0.6. | Deterministic female survival, free initial ages and IID age-5 recruitment around a constant with estimated SD. | Native juvenile survival and the age-1 stock–recruit relationship are not reproduced. |",
  "| F | Annual intensity and sex-specific logistic selectivity, annual slope/a50 deviations and male offset; ages 17+ shared. | Independent random walks in female F at age. | Source common-intensity, curve and sharing constraints are omitted. |",
  "| M | Both sexes' constant M estimated with priors; female estimate 0.192. | Estimate one female constant without the source prior. | Male M is outside the fitted population; no source-estimated M is turned into a fixed input. |",
  "| Catch and survey | Biomass and separate joint-sex age likelihoods; supplied composition weights and survey SEs. | Female numbers reconstructed from original totals, joint compositions and kg weights; estimated catch SD and supplied biomass survey log SD. | Males contribute to conversion but are not fitted. Composition covariance and ageing error are omitted. |",
  "| q and timing | Estimated scalar q with prior; female logistic selectivity. Biomass time 0.5, age proportions based on time-0 N; environmental effects inactive. | Free scalar and logistic q at time 0.5; no environmental effect. | Source q prior and different biomass/composition timings are not reproduced. |",
  "| Biology | Empirical female weights and maturity, February 1 female SSB. | Original female biology and time, without an extra half-sex multiplier. | Native combined-sex abundance, biomass and recruitment are not compared as female totals. Common SSB uses rounded official N/F/M and has no invented uncertainty. |",
  "", "Default starting values are tried first; a numerical failure permits one retry from aligned official female N and derived F. Missing source states retain package defaults. Data and assumptions do not change between attempts."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data=obs,years=years,ages=ages,
  N_settings=list(process="off",init="free"),
  rec_settings=list(process="iid",mu_form=~1),
  F_settings=list(process="rw"),
  M_settings=list(process="off",mu_form=~1,mu_supplied=NULL),
  catch_settings=list(sd_form=~1),
  index_settings=list(q_form=~logistic(age),sd_form=~0,sd_supplied=~log_sd),
  silent=silent,...
)
fit <- NULL
if(do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model,retry_start=function() {
    n <- .translation_start_surface(comparison_outputs,"population","numbers_at_age",years,ages,1e6)
    f <- .translation_start_surface(comparison_outputs,"mortality","fishing_mortality_at_age",years,ages)
    list(log_r0=log(n[1,1]),log_r=log(n[-1,1]),log_n0=log(n[1,-1]),log_f=log(f))
  },attempt_file=file.path("analysis/comp_assessments/results/cache",source$assessment$assessment_id,"initialization_attempts.rds"))
}
