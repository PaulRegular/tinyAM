source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1975:2024
ages <- 2:15
inputs <- source$inputs
inputs <- inputs[!is.na(inputs$year) & inputs$year %in% years, ]
p <- inputs$measure == "frequency_at_age"
totals <- ave(inputs$value[p], inputs$observation_id[p], FUN = sum)
inputs$value[p] <- inputs$value[p] / totals
inputs$measure[p] <- "proportion_at_age"
inputs$basis[p] <- "proportion_numbers"
inputs$unit[p] <- "proportion"
inputs$notes[p] <- paste(inputs$notes[p], "Translation: native percentages normalized jointly as in SS3; source percentages remain canonical.")
w <- inputs[inputs$type == "weight" & is.na(inputs$survey), ]
mat <- inputs[inputs$measure == "maturity_fecundity_at_age", ]
key <- function(x) paste(x$year, x$age)
mat$value <- mat$value / w$value[match(key(mat), key(w))]
mat$value[mat$age < 2] <- 0
mat$type <- "maturity"
mat$measure <- "maturity_at_age"
mat$basis <- mat$unit <- "proportion"
mat$notes <- "Translation: original maturity-fecundity weight divided by original population weight; source first mature age 2 applied. Female fraction 0.5 supplied separately."
inputs <- rbind(
  inputs[inputs$type != "biology", , drop = FALSE],
  mat
)
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages,
  maturity_multiplier = .5, assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0
obs$pop_covariates$M_prior_group <- factor("stock")
obs$index$q_block <- factor(pmin(obs$index$age, 6))

## Background and comparisons ----
age_plus_group <- 15
comparison_outputs <- source$outputs
ssb_age <- source$outputs[source$outputs$measure == "numbers_at_age" & source$outputs$year %in% years, ]
ssb_age$value <- ssb_age$value * 1000 * w$value[match(key(ssb_age), key(w))] *
  .5 * mat$value[match(key(ssb_age), key(mat))]
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Female begin-year mature biomass from marginal posterior median N and supplied deterministic biology; a sum of marginal medians need not equal the median of total SSB. No combined interval reconstructed."
comparison_outputs <- rbind(comparison_outputs, ssb_age)
comparison_definitions <- list(
  N = list(
    status = "approximate", definition = "Posterior median N at age, pooling source ages 15-20+",
    reason = "Ages 2-14 retain marginal medians; the 15+ sum of marginal medians is not necessarily the posterior median of pooled N."
  ),
  abundance = list(
    status = "approximate", definition = "Sum of marginal median N over modeled ages",
    reason = "Joint posterior draws unavailable; not a native posterior median of total abundance."
  ),
  biomass = list(
    status = "approximate", definition = "Sum of marginal biomass medians over modeled ages",
    reason = "A sum of age-specific posterior medians need not be the posterior median of aggregate biomass."
  ),
  biomass_at_age = list(
    status = "approximate", definition = "Marginal median biomass at age, pooling source ages 15-20+",
    reason = "The 15+ sum of marginal medians is not necessarily the median of pooled biomass."
  ),
  ssb = list(
    status = "approximate", definition = "Modeled-age female begin-year mature biomass from marginal N medians",
    reason = "Native posterior median SSB and its 95% credible interval remain separate; derived sums have no reconstructed uncertainty."
  ),
  recruitment = list(
    status = "matched", definition = "Age-2 entry from official N in the same calendar year",
    reason = "Native recruitment is age 0; comparison uses marginal posterior median N at age 2."
  )
)
background <- c(
  "### Joint U.S.–Canada coastal hake: official 2025 base model", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM candidate | Principal difference |",
  "|---|---|---|---|",
  "| Years and ages | 1966–2024, ages 0–20+, with observed ages 1–15+. | Fit 1975–2024, ages 2–15+. | Starts with fishery age data and the age-2+ acoustic definition; recruitment and total biomass are not strictly age-matched. Source biology is constant over ages 15–20+, supporting the translated 15+ group. |",
  "| N and recruitment | Deterministic cohort survival; Beverton–Holt recruitment at age 0 with estimated R0/steepness, fixed log SD 1.4, early cohort deviations and terminal rules. | N process off, exponential initial abundance, IID age-2 recruitment around a constant mean with fixed SD 1.4. | Age-2 recruitment includes juvenile survival; the source stock–recruit relationship, early-cohort priors, terminal rules and pre-1975 initialization differ. |",
  "| F | Annual fishing intensity with nonparametric selectivity; deviations at ages 2–6 and shared older increments. | Correlated random walks in log F-at-age (`process = 'cor_rw'`). | Permits coherent age-specific F changes but does not reproduce the source's selectivity normalization, age restrictions or temporal penalties. |",
  "| M | One constant estimated with a log-M penalty centered at -1.60944 (SD 0.1), bounded 0.05–0.4. | A single stock-wide Gaussian log-M effect with the same center and fixed SD; no M process. | Laplace-based fit rather than the source Bayesian integration; native bounds are not reproduced. |",
  "| Catch | Catch biomass and age percentages with a separate Dirichlet–multinomial composition likelihood; catch log SD fixed at 0.01. | Reconstructed catch-at-age in numbers with an estimated common lognormal SD. | Aggregate catch and composition likelihoods, ageing error and within-year composition covariance are simplified. |",
  "| Acoustic index errors | Age-2+ biomass survey at year fraction 0.5; supplied sampling error and estimated extra survey SD combined in quadrature, with a separate age-composition likelihood. | Reconstructed age-specific index; retained supplied `log_index_sd` via `sd_supplied = ~ relative_sd`. With `sd_form = ~1`, tinyAM drops the intercept and uses supplied SDs directly. | No fitted extra survey SD in this candidate; shared index error and composition covariance across ages are not reproduced. |",
  "| Survey q / selectivity | Analytically profiled scalar q (reported posterior median 0.894) and distinct time-invariant nonparametric age selectivity, with older increments fixed. | `q_form = ~ q_block` gives separate effective q at ages 2–5 and shared q at ages 6+, while `q_link = 'logit'` constrains every effective q to (0, 1). | Bounding q was pivotal for obtaining a reasonable abundance scale, but fitted q is near the upper boundary and its intervals cover nearly the allowed range. Effective q-at-age is not separately normalized into scalar q and selectivity, and is not directly equivalent to the profiled source scalar. |",
  "| Biology and SSB | Annual empirical population/fishery/survey weights and maturity–fecundity products; female fraction 0.5; first mature age 2; January 1 spawning. | Divide the source maturity–fecundity product by population weight; retain source weights, cutoff, female fraction and timing. | Common-age sums from marginal posterior medians are approximate and distinct from the source's posterior median aggregate SSB. |",
  "", "This candidate converged and gives a reasonably close approximation of the accepted assessment, as reported by the user. Constraining survey catchability below one was important for resolving the earlier abundance-scale discrepancy. The fitted effective q values remain weakly identified: estimates approach one and intervals extend toward both zero and one. That is a limitation of parameter inference, not evidence that the broad population comparison failed. No new fit statistics have been independently verified for this documentation update.",
  "", "The operational model profiles a scalar q and separately estimates age selectivity; tinyAM estimates bounded effective q by age block. The source reports Bayesian posterior medians, whereas tinyAM provides likelihood-based estimates and conditional random-effect modes. Neither the official F-at-age posterior nor joint posterior draws are available for exact comparisons. The fitted script makes one direct `fit_tam()` call without a source-state initialization retry."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "iid", mu_form = ~1, sd = 1.4),
    F_settings = list(process = "cor_rw"),
    M_settings = list(
      process = "off", mu_form = ~ 0 + iid(M_prior_group, sd = 0.1),
      mu_supplied = ~ I(exp(-1.60944))
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~ q_block,
      sd_form = ~ 1,
      sd_supplied = ~ relative_sd,
      q_link = "logit"
    ),
    silent = silent
  )
}
