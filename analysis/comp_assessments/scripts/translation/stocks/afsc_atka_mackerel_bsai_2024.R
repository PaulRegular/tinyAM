source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1990:2023
ages <- 2:11
obs <- database_to_tam_obs(source$assessment$assessment_id, source$inputs,
  years = years, ages = ages, maturity_multiplier = .5,
  assumptions = source$assumptions)
obs$maturity$spawn_time <- 7 / 12
obs$index$q_block <- factor(pmin(obs$index$age, 10))
cv <- source$inputs[source$inputs$measure == "relative_standard_error", ]
obs$index$relative_sd <- sqrt(log1p(cv$value[match(obs$index$year, cv$year)]^2))

## Background and comparisons ----
age_plus_group <- 11
comparison_outputs <- source$outputs
ssb_age <- source$outputs[source$outputs$measure == "numbers_at_age", ]
weights <- source$inputs[source$inputs$type == "weight" & is.na(source$inputs$survey), ]
maturity <- source$inputs[source$inputs$type == "maturity", ]
f <- source$outputs[source$outputs$measure == "fishing_mortality_at_age", ]
f_match <- match(paste(ssb_age$year, ssb_age$age), paste(f$year, f$age))
ssb_age$value <- ssb_age$value * 1000 * weights$value[match(ssb_age$age, weights$age)] *
  .5 * maturity$value[match(ssb_age$age, maturity$age)] *
  exp(-(f$value[f_match] + .3) * 7 / 12)
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "Common-definition August 1 female biomass calculated from rounded official N and F, fixed M, population weights and maturity. No uncertainty reconstructed."
comparison_outputs <- rbind(comparison_outputs, ssb_age)
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "August 1 female SSB at modeled ages from rounded official N/F and recorded biology",
  reason = "Calculated at report precision; native age-1+ SSB and its intervals remain separate."))
background <- c(
  "### BSAI Atka mackerel: official 2024 Model 16.0b", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1977–2024, ages 1–11+. | 1990–2023, ages 2–11+. | Restrict to the period with complete reported catch-age cells. Age 1 is absent from the catch-age table. The 2024 biomass observations have no matching age compositions, so are omitted rather than given invented ages. |",
  "| N and recruitment | Deterministic survival; equilibrium initial ages; Beverton–Holt recruitment with fixed steepness 0.8 and estimated recruitment SD. | Deterministic survival, free initial ages and IID recruitment around an estimated constant. | Age-2 entry includes prior juvenile survival and is not the native age-1 stock–recruit relationship. Comparable recruitment is official N at age 2 in the same calendar year. |",
  "| F | Annual intensity times flexible age selectivity, with age/time penalties and shared selectivity at ages 10+. | Independent temporal random walks for F at each age. | Source penalties and terminal-age F sharing are not reproduced. |",
  "| M | Fixed 0.3 at all ages. | Same fixed input. | Checked for input consistency, not credited as estimated agreement. |",
  "| Catch | Biomass and separate multinomial age compositions; catch CV 5%. | Numbers reconstructed from original biomass, reported compositions and annual catch weights; estimated lognormal SD. | Composition rounding and age restriction are explicit. The separate source likelihood and supplied composition weights are not reproduced. |",
  "| Index and q | One combined AI bottom-trawl biomass series plus age compositions; time 6.5/12; estimated scalar q and nonparametric selectivity with a prior. | Reconstructed age indices using original survey weights and time; q by age, shared at ages 10+; log SD supplied from the rounded index CV. | The q prior, ageing error and additional correlated composition uncertainty are omitted. No artificial age series is created for survey-only 2024. |",
  "| Biology and SSB | Separate original population, catch and survey weights; female maturity and August 1 female SSB. | Same recorded weights, maturity multiplied by 0.5 and spawning time 7/12. | Common-definition biomass and SSB use ages 2–11+; native age-1+ values remain available separately. |",
  "",
  "Official N is reported in whole millions; the comparison cannot resolve differences smaller than that reporting precision. Source F is reconstructed from two rounded output tables. Derived sums and F have no invented uncertainty.",
  "Default starting values are tried first. A single retry from aligned official N and F is allowed only after failed numerical checks."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ 1),
  F_settings = list(process = "rw"),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(q_form = ~ q_block, sd_form = ~ 0,
                        sd_supplied = ~ relative_sd),
  silent = silent, ...
)
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs, "population", "numbers_at_age",
                                     years, ages, 1e6)
    f <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age",
                                     years, ages)
    list(log_r0 = log(n[1, 1]), log_r = log(n[-1, 1]),
         log_n0 = log(n[1, -1]), log_f = log(f))
  }, attempt_file = file.path("analysis/comp_assessments/results/cache",
    source$assessment$assessment_id, "initialization_attempts.rds"))
}
