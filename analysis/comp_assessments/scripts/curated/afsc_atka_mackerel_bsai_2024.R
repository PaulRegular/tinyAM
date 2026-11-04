# BSAI Atka mackerel, 2024 operational update of AMAK Model 16.0b.
# Working specification: logistic mean F with AR1 deviations; monotone survey q plateau at age 7.
# No official fitted states are used to construct the observations.
## Observations ----
years <- 1990:2023
ages <- 2:11
obs <- database_to_tam_obs(
  source$assessment$assessment_id, source$inputs,
  years = years, ages = ages, maturity_multiplier = 0.5,
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 7 / 12

# Official survey selectivity is nonparametric, time-invariant, and may be dome-shaped.
# This tinyAM fit instead uses mono(q_block): non-decreasing effective q at ages
# 2-7, with ages 7+ coupled to the same q level to stabilize high-age estimates.
# The native assessment couples ages 10+, not ages 7+. q_age is retained here
# for comparison with earlier RW candidates; only q_block enters the current fit.
obs$index$q_age <- pmin(obs$index$age, 7)
obs$index$q_block <- factor(pmin(obs$index$age, 7))

# The report supplies rounded aggregate survey CVs. Convert CV to lognormal SD
# and repeat across reconstructed age observations as a documented approximation.
# These SDs do not capture the original age-composition covariance or effective N.
cv <- source$inputs[
  source$inputs$type == "index" &
    source$inputs$measure == "relative_standard_error" &
    source$inputs$survey == "AI_bottom_trawl", , drop = FALSE
]
if (anyDuplicated(paste(cv$year, cv$survey))) {
  stop("Atka survey has duplicate annual CV observations.", call. = FALSE)
}
cv_row <- match(
  paste(obs$index$year, obs$index$survey),
  paste(cv$year, cv$survey)
)
if (anyNA(cv_row) ||
    any(!is.finite(cv$value[cv_row]) | cv$value[cv_row] <= 0)) {
  stop("Atka survey CVs do not cover the reconstructed index years.", call. = FALSE)
}
obs$index$relative_sd <- sqrt(log1p(cv$value[cv_row]^2))

## Background and comparisons ----
age_plus_group <- 11
comparison_outputs <- source$outputs
ssb_age <- source$outputs[source$outputs$measure == "numbers_at_age", ]
weights <- source$inputs[
  source$inputs$type == "weight" & is.na(source$inputs$survey),
]
maturity <- source$inputs[source$inputs$type == "maturity", ]
f <- source$outputs[
  source$outputs$measure == "fishing_mortality_at_age",
]
f_match <- match(paste(ssb_age$year, ssb_age$age), paste(f$year, f$age))
ssb_age$value <- ssb_age$value * 1000 *
  weights$value[match(ssb_age$age, weights$age)] *
  0.5 * maturity$value[match(ssb_age$age, maturity$age)] *
  exp(-(f$value[f_match] + 0.3) * 7 / 12)
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- paste(
  "Common-definition August 1 female biomass calculated from rounded",
  "official N and F, fixed M, population weights and maturity.",
  "No uncertainty reconstructed."
)
comparison_outputs <- rbind(comparison_outputs, ssb_age)
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "August 1 female SSB at modeled ages from rounded official N/F and recorded biology",
    reason = "Calculated at report precision; native age-1+ SSB and its intervals remain separate."
  )
)

background <- c(
  "### BSAI Atka mackerel: official 2024 Model 16.0b", "", print_sources(source$assessment), "",
  "| Component | Official assessment | Current tinyAM translation | Key differences |",
  "|---|---|---|---|",
  "| Years / ages | 1977–2024, ages 1–11+. | 1990–2023, ages 2–11+. | The published early catch-at-age table contains incomplete rows; catch age 1 is unavailable. The 2024 totals have no matching age compositions and cannot enter this age-specific likelihood without assumptions. |",
  "| Population / recruitment | Deterministic survival; equilibrium initial conditions; age-1 Beverton–Holt with steepness 0.8 and an estimated recruitment SD. | `N process = 'off'`, `init = 'free'`; age-2 IID recruitment around an estimated constant mean, with estimated SD. | The source stock–recruit relationship, age-1 recruitment, early population history and bias adjustments are not reproduced. |",
  "| F | Annual fishing intensity multiplied by age selectivity that can change over time; flexible source penalties and shared selectivity at ages 10+. | `F process = 'ar1'`, `mu_form = ~ logistic(age)`: a rising logistic mean F-at-age with stationary age- and year-correlated log-F residuals. | Residual F can vary across ages and years, but this is not a normalized fishing-intensity × selectivity model, does not force shared F at 10+, and does not reproduce source penalty structure. |",
  "| M | Fixed M = 0.3 at every age. | `process = 'off'`, `mu_supplied = ~ M_assumption`; same M. | Agreement in M is input consistency, not independent parameter recovery. |",
  "| Catch | Total biomass with 5% CV, plus separate multinomial age compositions with supplied effective sample sizes. | Source biomass and age composition converted to catch-at-age using original catch weights; estimated common lognormal catch SD. | Separate biomass and composition likelihoods, covariance, sampling weights and years lacking compositions are not fully represented. |",
  "| Survey / q | One Aleutian Islands bottom-trawl biomass index and age composition; time-invariant nonparametric selectivity (potentially dome-shaped), with separately estimated scalar q, its prior and normalization across ages 4–10; source selectivity shared at ages 10+. | `q_form = ~ 1 + mono(q_block)`, with `q_block = factor(pmin(age, 7))` and `q_link = 'log'`. Effective q is non-decreasing through age 7 and identical from age 7 onward, with no year effect. | The age-7 plateau is a deliberate regularization to avoid extreme high-age q estimates, not a source assumption. The curve cannot descend at older ages and is not decomposed into normalized selectivity × scalar q or constrained to q < 1. |",
  "| Survey uncertainty | Year-specific supplied survey CVs for aggregate biomass plus a separate multinomial composition likelihood. | Rounded published CVs converted to log-observation SD by `sqrt(log1p(CV^2))`, passed as `sd_supplied = ~ relative_sd` with `sd_form = ~ 0`. | Source aggregate index uncertainty is applied to reconstructed age rows treated as independent; original age-composition uncertainty and covariance remain unrepresented. No SD multiplier is estimated. |",
  "| Biology / SSB | Population, catch and survey weights distinct; female maturity and 0.5 sex fraction; survey time 6.5/12, spawning at 7/12. | Retain original biological data and timing; `maturity_multiplier = 0.5`. | Approximate common-age SSB comparison is constructed from rounded official N/F; native aggregate SSB and its uncertainty remain separate. |",
  "", "The selected model is a single direct fit_tam call without a fallback initialization or source-fitted states supplied as starting values. The user reports this specification as a reasonable practical translation, including improved behaviour after coupling q at ages 7+; current optimizer diagnostics and numerical comparison metrics were not provided with the revised script.",
  "", "The earlier objective 511.565 and comparison summaries (SSB MAPD 25.7%, terminal difference -23.6%, correlation 0.966) belong to an older independent-F-RW / independent age-q fit and must NOT be attributed to this AR1-F / monotone-q specification.",
  "", "The source allows more flexible, potentially dome-shaped survey selectivity and time-varying fishery selectivity. Here, monotone effective q regularizes the high-age behaviour while stationary age/year F deviations approximate changing selectivity. The model remains a simplified age-specific lognormal reconstruction rather than the accepted separate aggregate and age-composition likelihoods."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()

  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "iid", mu_form = ~ 1),
    F_settings = list(process = "ar1", mu_form = ~logistic(age)),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~ M_assumption
    ),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(
      q_form = ~ 1 + mono(q_block),
      q_link = "log",
      sd_form = ~ 0,
      sd_supplied = ~ relative_sd
    ),
    silent = silent
  )
}
