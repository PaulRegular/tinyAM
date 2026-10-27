source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1975:2024
ages <- 2:15
inputs <- source$inputs
inputs <- inputs[!is.na(inputs$year) & inputs$year %in% years, ]
p <- inputs$measure=="frequency_at_age"
totals <- ave(inputs$value[p],inputs$observation_id[p],FUN=sum)
inputs$value[p] <- inputs$value[p]/totals
inputs$measure[p] <- "proportion_at_age"
inputs$basis[p] <- "proportion_numbers"
inputs$unit[p] <- "proportion"
inputs$notes[p] <- paste(inputs$notes[p],"Translation: native percentages normalized jointly as in SS3; source percentages remain canonical.")
w <- inputs[inputs$type=="weight" & is.na(inputs$survey), ]
mat <- inputs[inputs$measure=="maturity_fecundity_at_age", ]
key <- function(x) paste(x$year,x$age)
mat$value <- mat$value/w$value[match(key(mat),key(w))]
mat$value[mat$age<2] <- 0
mat$type <- "maturity"
mat$measure <- "maturity_at_age"
mat$basis <- mat$unit <- "proportion"
mat$notes <- "Translation: original maturity-fecundity weight divided by original population weight; source first mature age 2 applied. Female fraction 0.5 supplied separately."
inputs <- rbind(inputs[!inputs$type %in% c("biology") & inputs$measure!="log_index_sd", ],mat)
obs <- database_to_tam_obs(source$assessment$assessment_id,inputs,years=years,ages=ages,
  maturity_multiplier=.5,assumptions=source$assumptions)
obs$maturity$spawn_time <- 0
obs$pop_covariates$M_prior_group <- factor("stock")
obs$index$q_block <- factor(pmin(obs$index$age,6))

## Background and comparisons ----
age_plus_group <- 15
comparison_outputs <- source$outputs
ssb_age <- source$outputs[source$outputs$measure=="numbers_at_age" & source$outputs$year %in% years, ]
ssb_age$value <- ssb_age$value*1000*w$value[match(key(ssb_age),key(w))] *
  .5*mat$value[match(key(ssb_age),key(mat))]
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Female begin-year mature biomass from marginal posterior median N and supplied deterministic biology; a sum of marginal medians need not equal the median of total SSB. No combined interval reconstructed."
comparison_outputs <- rbind(comparison_outputs,ssb_age)
comparison_definitions <- list(
  N=list(status="approximate",definition="Posterior median N at age, pooling source ages 15-20+",
    reason="Ages 2-14 retain marginal medians; the 15+ sum of marginal medians is not necessarily the posterior median of pooled N."),
  abundance=list(status="approximate",definition="Sum of marginal median N over modeled ages",
    reason="Joint posterior draws unavailable; not a native posterior median of total abundance."),
  biomass=list(status="approximate",definition="Sum of marginal biomass medians over modeled ages",
    reason="A sum of age-specific posterior medians need not be the posterior median of aggregate biomass."),
  biomass_at_age=list(status="approximate",definition="Marginal median biomass at age, pooling source ages 15-20+",
    reason="The 15+ sum of marginal medians is not necessarily the median of pooled biomass."),
  ssb=list(status="approximate",definition="Modeled-age female begin-year mature biomass from marginal N medians",
    reason="Native posterior median SSB and its 95% credible interval remain separate; derived sums have no reconstructed uncertainty."),
  recruitment=list(status="matched",definition="Age-2 entry from official N in the same calendar year",
    reason="Native recruitment is age 0; comparison uses marginal posterior median N at age 2."))
background <- c(
  "### Joint U.S.–Canada coastal hake: official 2025 base model", "",print_sources(source$assessment),"",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1966–2024, modeled ages 0–20+, observed ages 1–15+. | 1975–2024, ages 2–15+. | Begin with the first fishery age samples and the acoustic age-2+ definition. Totals are converted before filtering ages. Original weights and fecundity are identical at ages 15–20+, allowing the observation plus group. The 2026 update supplies newer projections, not a different historical fit. |",
  "| N and recruitment | Deterministic survival; early cohort deviations; Beverton–Holt with estimated R0/steepness, fixed log SD 1.4 and terminal recruitment tied to the mean curve. | Deterministic survival, free initial ages, age-2 recruitment with a two-year-lagged Beverton–Holt mean and IID deviations at fixed log SD 1.4. | Age-2 entry includes juvenile survival; early-cohort priors and terminal fixing are not reproduced. |",
  "| F | Annual intensity and nonparametric fishery selectivity, with time variation at ages 2–6 and shared older increments. | Independent temporal F random walks at all ages. | The source selectivity normalization, age sharing and fixed deviation SD are not reproduced. |",
  "| M | One constant estimated with log-M penalty centered at -1.60944, SD 0.1, and bounds 0.05–0.4. | One shared Gaussian log-M effect with the same center and fixed SD. | The density contributes once per stock. M is estimated with Laplace integration; source bounds and full Bayesian posterior integration are not reproduced. |",
  "| Catch and index | Catch biomass plus percentage age compositions; acoustic age-2+ biomass at time 0.5; separate Dirichlet-multinomial likelihoods and estimated extra survey SD. | Numbers reconstructed from original totals, normalized native percentages and original fleet weights; estimated catch and survey log SDs. | Ageing error, composition covariance and the source precision model are omitted. No fitted parameter adjusts an observation. |",
  "| q | Analytically profiled scalar and free nonparametric age increments through age 6, older increments zero. | Free q by age 2–6, shared at ages 6+. | Preserve the source sharing and allow nonmonotonic shapes; profiling, bounds and source uncertainty differ. q is estimated, not fixed to one. |",
  "| Biology | Original annual empirical weights and maturity-fecundity products; female fraction 0.5; first mature age 2; January 1 female SSB. | Divide the supplied product by population weight, apply the same mature-age cutoff, female fraction and time. | Native posterior median SSB is retained separately from common-age sums based on marginal N medians. No intervals are invented for sums or plus groups. |",
  "", "Official points are Bayesian posterior medians; tinyAM reports likelihood-based estimates and conditional random-effect modes. Joint posterior F-at-age is unavailable, and exploitation fractions or annual F intensity are not substituted for it. Default starts are tried first; one conditional retry uses aligned N, leaving other starts at package defaults. The prior center is a source assumption, not the fitted official M estimate."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data=obs,years=years,ages=ages,
  N_settings=list(process="off",init="free"),
  rec_settings=list(process="iid",mu_form=~bh(ssb,lag=2),sd=1.4),
  F_settings=list(process="rw"),
  M_settings=list(process="off",mu_form=~0+iid(M_prior_group,sd=.1),
                  mu_supplied=~I(exp(-1.60944))),
  catch_settings=list(sd_form=~1),
  index_settings=list(q_form=~0+q_block,sd_form=~1),
  silent=silent,...
)
retry_start <- function() {
  n <- source$outputs[source$outputs$measure=="numbers_at_age", ]
  n$age <- pmin(n$age,15)
  n <- aggregate(value~year+age,n,sum)
  n <- n[n$year %in% years & n$age %in% ages & is.finite(n$value) & n$value>0, ]
  # Named partial vectors let the package retain defaults for rounded zeros.
  r <- n[n$age==min(ages) & n$year>min(years), ]
  n0 <- n[n$year==min(years) & n$age>min(ages), ]
  start <- list(log_r=setNames(log(r$value*1e3),as.character(r$year)),
    log_n0=setNames(log(n0$value*1e3),as.character(n0$age)))
  r0 <- n$value[n$year==min(years) & n$age==min(ages)]
  if(length(r0)) start$log_r0 <- log(r0*1e3)
  start
}
fit <- NULL
if(do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model,retry_start=retry_start,
    attempt_file=file.path("analysis/comp_assessments/results/cache",source$assessment$assessment_id,"initialization_attempts.rds"))
}
