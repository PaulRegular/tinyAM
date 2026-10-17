source("analysis/comp_assessments/R/fit_default_first.R", local=TRUE)

## Observations and biology ----
years <- 2000:2023
ages <- 1:12
inputs <- source$inputs
key_rows <- inputs[inputs$type=="catch" & inputs$measure=="conditional_proportion_at_age",]
lengths <- sort(unique(key_rows$length_bin_lower))
key <- matrix(0,length(lengths),13,dimnames=list(lengths,0:12))
for(length in lengths) {
  z <- key_rows[key_rows$length_bin_lower==length,]
  z$value <- z$value/ave(z$value,z$observation_id,FUN=sum)
  count <- tapply(z$value*z$sample_size,z$age,sum)
  key[as.character(length),] <- count[as.character(0:12)]/sum(count)
}
length_obs <- inputs[inputs$type=="catch" & inputs$measure=="proportion_at_length" & inputs$year %in% years,]
catch <- lapply(split(length_obs,length_obs$year),function(z) {
  z <- z[order(z$length_bin),]
  p <- z$value/sum(z$value)
  covered <- z$length_bin %in% lengths
  weight <- 5.40706e-6*(z$length_bin+2.5)^3.19601
  biomass <- inputs$value[inputs$type=="catch" & inputs$measure=="total_biomass" & inputs$year==z$year[1]]
  stopifnot(length(biomass)==1,sum(p*weight)>0)
  numbers <- biomass*1000/sum(p*weight)
  age_p <- as.numeric(p[covered] %*% key[as.character(z$length_bin[covered]),,drop=FALSE])
  out <- z[rep(1,13),]
  out$age <- 0:12;out$value <- numbers*age_p
  out$measure <- "numbers_at_age";out$basis <- "numbers";out$unit <- "fish"
  out$notes <- paste("Translation: pooled supplied-N-weighted empirical fishery age-length key; uncovered number fraction",
    signif(1-sum(p[covered]),5),"left unallocated. Total numbers use original fixed weight-length relation at 5-cm bin midpoints.")
  out
})
catch <- do.call(rbind,catch)
survey <- inputs[inputs$type=="index" & inputs$measure=="proportion_at_age" &
  grepl("^Active",inputs$notes) & inputs$year %in% years,]
survey$value <- survey$value/ave(survey$value,survey$observation_id,FUN=sum)
total <- inputs[inputs$type=="index" & inputs$measure=="total_numbers",]
survey$value <- survey$value*total$value[match(survey$year,total$year)]*1000
survey$measure <- "numbers_at_age";survey$basis <- "numbers";survey$unit <- "fish"
survey$notes <- "Translation: active marginal survey number proportions times original VAST total abundance; no fitted parameter adjusts observations."

outputs <- source$outputs
state_key <- function(z)paste(z$year,z$age)
n <- outputs[outputs$measure=="numbers_at_age" & outputs$year %in% years,]
w <- outputs[outputs$measure=="weight_at_age" & outputs$year %in% years,]
w <- w[match(state_key(n),state_key(w)),]
mat <- plogis(.132*((w$value/5.40706e-6)^(1/3.19601)-58))
mat[n$age<1] <- 0
biology <- lapply(split(seq_len(nrow(n)),paste(n$year,pmin(n$age,12))),function(j) {
  out <- inputs[rep(which(inputs$type=="M")[1],2),]
  out$year <- n$year[j[1]];out$year_basis <- "calendar_year";out$age <- min(n$age[j])
  out$type <- c("weight","maturity");out$measure <- c("weight_at_age","maturity_at_age")
  out$value <- c(weighted.mean(w$value[j],n$value[j]),
    weighted.mean(mat[j],n$value[j]*w$value[j]))
  out$unit <- c("kg/fish","proportion");out$basis <- c("kg_per_fish","proportion")
  out$source_type <- "derived_source_output"
  out$notes <- "Translation only: fitted population weights; maturity proxy from weight-equivalent length, not official integrated maturity. Age-12+ biology weighted by official N; no fitted states used to allocate source observations."
  out
})
biology <- do.call(rbind,biology)
obs <- database_to_tam_obs(source$assessment$assessment_id,
  rbind(catch,survey,biology,inputs[inputs$type=="M",]),years=years,ages=ages,
  maturity_multiplier=.5,assumptions=source$assumptions)
obs$maturity$spawn_time <- 0

## Comparisons and background ----
age_plus_group <- 12
comparison_outputs <- outputs
ssb_age <- n
ssb_age$value <- n$value*w$value*.5*mat
ssb_age$type <- "biomass";ssb_age$measure <- "mature_biomass_at_age";ssb_age$unit <- "t"
ssb_age$se <- ssb_age$lwr <- ssb_age$upr <- NA_real_
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Common-definition female January biomass using weight-derived maturity proxy; not official length-integrated female SSB. Thousand fish times kg gives tonnes."
comparison_outputs <- rbind(comparison_outputs,ssb_age)
comparison_definitions <- list(
  F=list(status="approximate",definition="F from beginning/midyear official N, pooled at age 12+",reason="Rounded N limits precision; tail F is N-weighted, without covariance."),
  F_bar=list(status="approximate",definition="Population-weighted mean of recovered F",reason="Derived from rounded beginning/midyear N; not native apical F."),
  ssb=list(status="approximate",definition="Common-age female biomass with shared maturity proxy",reason="Original maturity integrates over length; native female SSB and reported SD remain separate."),
  recruitment=list(status="matched",definition="Same-calendar-year age-1 entry from official N",reason="Native recruitment is age 0."))
background <- c(
  "### Eastern Bering Sea Pacific cod: official 2024 Model 24.1","",print_sources(source$assessment),"",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Coverage | 1977–2024, ages 0–20+, combined sex. | 2000–2023, ages 1–12+. | Use active survey ages; 2020 survey ages are unavailable. Original totals are converted before age restriction. |",
  "| Catch | Catch tonnes and length compositions; fishery conditional ages inactive. | Allocate original length proportions with a pooled empirical key from those biological samples, weighted by supplied N. | This adds a fixed age allocation outside the official likelihood. Uncovered length bins, 0.003–0.54% of number proportions, remain unallocated. No fitted N or q adjusts observations. |",
  "| Weight conversion | Fixed weight–length relation and estimated growth. | Catch total numbers use kg=5.40706e-6 × (lower length edge + 2.5 cm)^3.19601. | Five-cm midpoint approximation; within-bin length distribution and key uncertainty omitted. |",
  "| Survey and q | July 1 VAST total numbers plus active marginal ages; length selectivity, ageing error, estimated scalar q and extra SD. | Original total numbers times active age proportions at time 0.5; free q at each age and estimated survey log SD. | No biomass conversion. Shared precision, length selectivity, age error and q bounds are omitted. |",
  "| N and recruitment | Deterministic survival; steepness fixed 1, fixed recruitment log SD 0.6908, early regime and bias adjustments. | Deterministic survival, free initial ages, constant-mean IID age-1 recruitment with log SD 0.6908. | Juvenile survival, early-regime and bias adjustments differ. |",
  "| F and M | Time-varying length selectivity and annual F; fixed M 0.3866. | Independent age-specific F random walks; the same fixed M. | Source selectivity sharing and penalties are omitted; fixed M is an input check only. |",
  "| Biology and SSB | Estimated growth/weights, length-integrated maturity, January female SSB. | Plug-in fitted weights; maturity proxy at weight-equivalent length, female fraction 0.5, age-12+ N-weighted biology. | Original empirical stock biology is unavailable. Common-definition SSB is approximate; official female SSB remains separately reported. Worksheet combined mature biomass is not mistaken for female SSB. |",
  "", "Native totals and compositions remain canonical. Source zero observations remain zero; tinyAM treats lognormal zero observations as missing. Default starts first; only a numerical failure allows one source-N/F retry. No assumptions change between attempts."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data=obs,years=years,ages=ages,
  N_settings=list(process="off",init="free"),
  rec_settings=list(process="iid",mu_form=~1,sd=.6908),
  F_settings=list(process="rw"),
  M_settings=list(process="off",mu_supplied=~M_assumption),
  catch_settings=list(sd_form=~1),
  index_settings=list(q_form=~factor(age),sd_form=~1),silent=silent,...
)
fit <- NULL
if(do_fit) {
  fit_stage <- "fit";fit_started <- Sys.time()
  fit <- fit_default_first(fit_model,retry_start=function() {
    states <- n;states$age <- pmin(states$age,12)
    n_sum <- aggregate(value~year+age,states,sum)
    nn <- .translation_start_surface(transform(n_sum,type="population",measure="numbers_at_age"),
      "population","numbers_at_age",years,ages,1)
    f <- outputs[outputs$measure=="fishing_mortality_at_age" & outputs$year %in% years & outputs$age>=1,]
    f$value <- f$value*n$value[match(state_key(f),state_key(n))]
    f$age <- pmin(f$age,12);f <- aggregate(value~year+age,f,sum)
    f$value <- f$value/n_sum$value[match(state_key(f),state_key(n_sum))]
    ff <- .translation_start_surface(transform(f,type="mortality",measure="fishing_mortality_at_age"),
      "mortality","fishing_mortality_at_age",years,ages)
    list(log_r0=log(nn[1,1]),log_r=log(nn[-1,1]),log_n0=log(nn[1,-1]),log_f=log(ff))
  },attempt_file=file.path("analysis/comp_assessments/results/cache",source$assessment$assessment_id,"initialization_attempts.rds"))
}
