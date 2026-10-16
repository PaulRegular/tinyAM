source("analysis/comp_assessments/R/fit_default_first.R",local=TRUE)

## Observations ----
years <- 2003:2022
ages <- 5:60
inputs <- source$inputs
mid_weight <- inputs[inputs$type=="weight" & inputs$season=="mid_year",]
compositions <- inputs[inputs$measure=="proportion_at_age" & inputs$year %in% years,]
observations <- lapply(split(compositions,compositions$observation_id),function(p) {
  total <- inputs[inputs$type==p$type[1] & inputs$measure=="total_biomass" & inputs$year==p$year[1],]
  total <- if(p$type[1]=="catch")total[total$fleet==p$fleet[1],] else total[total$survey==p$survey[1],]
  if(!nrow(total))return(NULL)
  w <- mid_weight$value[match(paste(p$sex,p$age),paste(mid_weight$sex,mid_weight$age))]
  stopifnot(nrow(total)==1,all(is.finite(w)),sum(p$value*w)>0)
  p$value <- total$value*1000*p$value/sum(p$value*w)
  p <- p[p$sex=="female",]
  p$measure <- "numbers_at_age";p$basis <- "numbers";p$unit <- "fish"
  if(p$type[1]=="index")p$sampling_time <- total$sampling_time
  p$notes <- "Translation: original area biomass times joint-sex number composition, divided by joint mean fixed source weight; both sexes/all ages enter before retaining female ages."
  p
})
observations <- do.call(rbind,observations)
catch <- observations[observations$type=="catch",]
complete <- vapply(split(catch,catch$year),function(p)length(unique(p$fleet))==3,logical(1))
catch <- catch[catch$year %in% as.numeric(names(complete)[complete]),]
catch <- aggregate(value~year+age,catch,sum)
catch_rows <- observations[rep(which(observations$type=="catch")[1],nrow(catch)),]
catch_rows$year <- catch$year;catch_rows$age <- catch$age;catch_rows$value <- catch$value
catch_rows$fleet <- "Pooled_trawl";catch_rows$region <- "BC coast";catch_rows$sampling_time <- NA_real_
catch_rows$notes <- "Translation: sum female numbers after separate area conversions; all three area compositions required. Missing area coverage gives missing coastwide catch, not a partial total."
biology <- inputs[(inputs$type=="weight" & inputs$season=="begin_year" | inputs$type=="maturity") & inputs$sex=="female",]
biology$season <- NA_character_
obs <- database_to_tam_obs(source$assessment$assessment_id,
  rbind(catch_rows,observations[observations$type=="index",],biology),
  years=years,ages=ages,assumptions=source$assumptions)
obs$maturity$spawn_time <- 0

## Comparisons and background ----
age_plus_group <- 60
comparison_outputs <- source$outputs[source$outputs$sex=="female" & source$outputs$region=="BC coast",]
comparison_outputs$lwr <- comparison_outputs$upr <- NA_real_
comparison_definitions <- list(
  M=list(status="matched",definition="Female constant M shared across source areas",reason="Official posterior median versus tinyAM likelihood estimate. Native 90% credible limits stay canonical; omitted here to avoid labelling them as 95%."),
  ssb=list(status="non_equivalent",definition="Female mature biomass versus unresolved source reproductive output",reason="Only 2024 boundary summary recovered; native control fecundity/first-mature-age settings differ from weight-times-maturity. No historical posterior N or SSB series available."))
background <- c(
  "### British Columbia Pacific ocean perch: official three-area 2023 R21v3","",print_sources(source$assessment),"",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Coverage | 1935–2023, ages 0–60+, female/male in three nonmixing areas 5ABC, 3CD, 5DE. | 2003–2022, female ages 5–60+, pooled coastwide. | Keep one case. Start with modern synoptic surveys and exclude provisional repeated 2023 catch. Male data enter number conversion, then females are retained. |",
  "| Catch | Separate area catch biomass and joint-sex age compositions, supplied N, Francis weights and ageing error. | Convert each area separately with fixed source sex weights, retain females, then sum only complete three-area years. | No missing area is treated as zero or a coastwide total. Composition precision/covariance and ageing error are omitted. |",
  "| Survey | Area-specific biomass and joint-sex age samples; different age selectivity curves. | Female numbers from original totals/compositions at native time 5/12. Free survey intercepts and survey-specific natural cubic age splines (4 df). | A pooled population assumes constant area shares within each survey's q. Source spatial recruitment/selectivity constraints and priors are not reproduced. |",
  "| N and recruitment | Deterministic sex/area survival, early deviations, Beverton–Holt with estimated steepness, fixed log SD 0.9 and spatial recruitment allocation. | Deterministic female survival, free initial ages, constant-mean IID age-5 recruitment with log SD 0.9. | Spatial allocation, stock–recruit feedback and bias adjustments differ. |",
  "| F and M | Area-specific F/selectivity and separate estimated constant M by sex, M prior mean 0.06 SD 0.018. | Independent female F random walks and one estimated female M without source prior/bounds. | Male M is outside the fitted population. Source-estimated M is not a fixed input. |",
  "| Biology | Fixed sex-specific growth, length variation/allometry and female maturity; no empirical weight matrix. | Recover fixed weights with matching SS3 at zero optimization evaluations; use female beginning-year population weights and original maturity at January 1. | Generated native N/F/M are never reference estimates or observation adjustments. No additional half-sex multiplier. |",
  "| Reference | Bayesian medians with 90% intervals; reported 2024 boundary SSB, sex-specific M, no recovered annual posterior N/F/R/SSB files. | Compare female M; show unavailable metrics clearly. Native reported SSB stays as context and is non-equivalent. | Public control uses fecundity option 2, coefficients 1/0 and first mature age read as 9; reproductive output differs from mature biomass. This source-definition issue remains unresolved. Native 90% limits are not relabelled 95%. |",
  "", "Default package starts first. Only a numerical failure permits one retry initialized at the official female M median; unavailable N/F starts remain defaults. Zero source age proportions remain zero and are treated as missing by the tinyAM lognormal observation likelihood."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data=obs,years=years,ages=ages,
  N_settings=list(process="off",init="free"),
  rec_settings=list(process="iid",mu_form=~1,sd=.9),
  F_settings=list(process="rw"),
  M_settings=list(process="off",mu_form=~1,mu_supplied=NULL),
  catch_settings=list(sd_form=~1),
  index_settings=list(q_form=~survey*splines::ns(age,df=4),sd_form=~survey),silent=silent,...
)
fit <- NULL
if(do_fit) {
  fit_stage <- "fit";fit_started <- Sys.time()
  fit <- fit_default_first(fit_model,retry_start=function()list(mu_m=log(.05229)),
    attempt_file=file.path("analysis/comp_assessments/results/cache",source$assessment$assessment_id,"initialization_attempts.rds"))
}
