source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1982:2023
ages <- 2:20
inputs <- source$inputs
compositions <- inputs[inputs$measure == "proportion_at_age" &
  inputs$sex %in% c("female","male") & inputs$year %in% years, ]
observations <- lapply(split(compositions,compositions$observation_id),function(p) {
  total <- inputs[inputs$type==p$type[1] & inputs$measure=="total_biomass" &
    inputs$year==p$year[1], ]
  # Some source age-sample years have no corresponding biomass index.
  if(!nrow(total)) return(NULL)
  weight <- inputs[inputs$type==if(p$type[1]=="catch")"catch_weight" else "weight", ]
  weight <- weight[weight$year==p$year[1] &
    if(p$type[1]=="catch") !is.na(weight$fleet) else !is.na(weight$survey), ]
  w <- weight$value[match(paste(p$sex,p$age),paste(weight$sex,weight$age))]/1000
  stopifnot(nrow(total)==1,all(is.finite(w)),sum(p$value*w)>0)
  p$value <- total$value*1e6*p$value/sum(p$value*w)
  p <- p[p$sex=="female", ]
  p$measure <- "numbers_at_age"
  p$basis <- "numbers"
  p$unit <- "fish"
  p$sampling_time <- if(p$type[1]=="index") .5 else NA_real_
  p$notes <- "Translation only: female component of total biomass times joint-sex age proportions divided by their sex-specific weighted mean. Both sexes enter the denominator; male observations are then omitted."
  p
})
observations <- do.call(rbind,observations)
biology <- inputs[inputs$sex=="female" & inputs$type %in% c("weight","maturity","M"), ]
obs <- database_to_tam_obs(source$assessment$assessment_id,rbind(observations,biology),
  years=years,ages=ages,assumptions=source$assumptions)
obs$maturity$spawn_time <- 4/12
for(name in c("temperature","start_date","interaction")) {
  z <- inputs[inputs$type=="covariate" & inputs$observation_id==name, ]
  obs$index[[name]] <- (z$value-mean(z$value))[match(obs$index$year,z$year)]
}
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
ssb_age$value <- ssb_age$value*1e6*w$value[match(key(ssb_age),key(w))]/1000 *
  mat$value[match(key(ssb_age),key(mat))] * exp(-(.12+f$value[match(key(ssb_age),key(f))])*4/12)
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Female May 1 mature biomass from native sex-specific N/F, original weights and maturity, fixed female M. Billion fish times grams times 1000 gives tonnes. No uncertainty reconstructed."
comparison_outputs <- rbind(comparison_outputs,ssb_age)
comparison_definitions <- list(
  N=list(status="matched",definition="Female begin-year abundance at age",reason="Male states remain canonical but are not represented in this translation."),
  abundance=list(status="matched",definition="Total female abundance at modeled ages",reason="Not combined-sex stock abundance."),
  recruitment=list(status="matched",definition="Female age-2 entry in the same calendar year",reason="Native recruitment is combined-sex age 1; comparison uses official female N at age 2."),
  ssb=list(status="approximate",definition="Female May 1 SSB at modeled ages",reason="Calculated from rounded official state reports without covariance."))
background <- c(
  "### Eastern Bering Sea yellowfin sole: official 2024 Model 23.0", "",print_sources(source$assessment),"",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1954–2024, ages 1–20+, both sexes. | 1982–2023, female ages 2–20+. | Use the VAST survey period with matching age samples. The source has zero age-1 catch weights; age 1 is omitted rather than assigned invented weights. This is one stock case, not two sex cases. |",
  "| N and recruitment | Deterministic sex-specific survival, free initial ages, fitted Ricker relationship over 1978–2018 and fixed log SD 0.6. | Deterministic female survival, free initial ages and IID age-2 recruitment around a constant mean with estimated SD. | Juvenile survival and the age-1 Ricker relationship are not reproduced. |",
  "| F | Annual intensity and ascending fishery selectivity with time deviations; ages 17+ shared. | Independent random walks in female F at age. | Source common-intensity and age-sharing constraints are omitted. |",
  "| M | Female M fixed 0.12; male M estimated. | Female M fixed 0.12. | Estimated male mortality is retained in the database but is not compared with female fixed M. |",
  "| Catch | Total biomass and joint-sex age compositions. | Female age numbers reconstructed from total biomass, joint proportions and sex-specific catch weights. | Male cells contribute to biomass conversion but are omitted from fitting. Separate composition likelihood and weights are not reproduced. |",
  "| Survey and q | VAST biomass at time 0.5; joint age proportions predicted from time-0 N. Estimated scalar log q depends on temperature, start-date and interaction. | Female age indices at time 0.5, ascending logistic q with the same three centered supplied covariates; supplied biomass log SD. | The source's different timing for biomass and proportions is approximate; additional composition uncertainty, source q priors and sex structure are omitted. |",
  "| Biology and SSB | Original empirical sex-specific weights and female maturity; May 1 female SSB. | Same original female weights, maturity and timing; no extra sex multiplier. | Comparisons cover female modeled ages. Combined-sex biomass and recruitment remain native database outputs, not directly comparable female metrics. |",
  "", "Default starts are tried first. One retry from aligned official female N and F is allowed only after numerical failure. Source observations and assumptions remain unchanged between attempts."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data=obs,years=years,ages=ages,
  N_settings=list(process="off",init="free"),
  rec_settings=list(process="iid",mu_form=~1),
  F_settings=list(process="rw"),
  M_settings=list(process="off",mu_form=NULL,mu_supplied=~M_assumption),
  catch_settings=list(sd_form=~1),
  index_settings=list(q_form=~temperature+start_date+interaction+logistic(age),
    sd_form=~0,sd_supplied=~log_sd),
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
