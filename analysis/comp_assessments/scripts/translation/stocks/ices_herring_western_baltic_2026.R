## Observations ----
years <- 1991:2025

# Use age 5 as the population plus group so all survey plus groups
# fall at the population plus age. The translation helper aggregates
# older catch/index observations into 5+.
ages <- 0:5

assessment_id <- source$assessment$assessment_id
inputs <- source$inputs

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = c("HERAS", "GERAS", "N20", "IBTS/BITSQ1"),
  assumptions = source$assumptions
)

catch_block <- cut_ages(obs$catch$age, c(0, 1, 2, max(ages)))
levels(catch_block) <- c("sam_sd_5", "sam_sd_6", "sam_sd_0")
obs$catch$sd_block <- factor(as.character(catch_block))
geras <- obs$index$survey == "GERAS"
geras_q <- cut_ages(obs$index$age[geras], 1:3)
levels(geras_q) <- c("sam_q_2", "sam_q_2", "sam_q_7")
obs$index$q_key <- with(obs$index, ifelse(
  survey == "HERAS", ifelse(age == 2, "sam_q_0", "sam_q_1"),
  ifelse(survey == "N20", "sam_q_3", paste0("sam_q_", age + 1L))
))
obs$index$q_key[geras] <- as.character(geras_q)
obs$index$q_key <- factor(obs$index$q_key)
obs$index$sd_block <- factor(ifelse(
  obs$index$survey == "HERAS", "sam_sd_1",
  ifelse(obs$index$survey == "GERAS", "sam_sd_2",
    ifelse(obs$index$survey == "N20", "sam_sd_3", "sam_sd_4")
  )
))
obs$index$relative_sd[!is.finite(obs$index$relative_sd)] <- 1


## Background and comparisons ----

age_plus_group <- max(ages)

comparison_scales <- c(
  N = 0.001, recruitment = 0.001, ssb = 0.001, biomass = 0.001, abundance = 0.001,
  biomass_at_age = 0.001, F = 1, M = 1, F_bar = 1
)


background <- c(
  "### Western Baltic spring-spawning herring: official 2026 SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main simplification |",
  "|---|---|---|---|",
  "| Years and ages | 1991–2025; ages 0–8+, with recruitment at age 0. | Fit ages 0–5+, aggregating source catch and survey observations at ages 5 and older into the new plus group. | The reduced age structure avoids survey plus groups falling below the population plus group, but does not reproduce all age-specific dynamics of ages 5–8+. |",
  "| N | SAM uses a segmented (hockey-stick) stock-recruitment relationship, separate recruitment-process variance, and shared survival-process variance across ages 1–8. | Use IID older-cohort survival deviations, free initial abundance, and a separate recruitment process. | The hockey-stick relationship, reduced age structure, and initial-state treatment remain different. |",
  "| F | SAM estimates age-specific F with AR(1) dependence across ages, and shares the F state for ages 7 and 8. Fbar is ages 2–5. | Use correlated F random walks (cor_rw) over ages 0–5+, with mean_ages = 2:5. | Cross-age dependence is better approximated, but the source age-specific F structure and shared oldest-age state are not reproduced exactly. The translated age-5+ F is not directly equivalent to source age-5 F. |",
  "| M | Fixed, time-invariant age-specific M derived from NSAS herring and profiled at the 2025 benchmark. | Retain the published M values through age 5. | The age-5 M value is applied to the translated 5+ group rather than constructing an aggregate mortality schedule for ages 5–8+. |",
  "| Catch | Total catch-at-age in thousands of fish, with lognormal observation errors grouped for ages 0, 1, and 2–8. | Convert to individual fish, sum ages 5–8+ into 5+, and retain the three catch observation-SD groups. | Catch totals are preserved under aggregation, but the original age-specific observations and error structure are not. |",
  "| Indices | HERAS ages 2–6+, GERAS ages 1–3, N20 age 0, and IBTS/BITS Q1 ages 3–5+. Sampling fractions are 0.625, 0.8, 0.4, and 0.136365, respectively. | Combine HERAS ages 5 and 6+ into one 5+ index; retain the native IBTS/BITS 5+ index. Other survey ages remain unchanged. Preserve sampling fractions, q-sharing groups, and available relative SD factors. | Aggregated HERAS observations do not reproduce the original joint lognormal likelihood or cross-age residual correlation. Five HERAS observations are missing in 1999. |",
  "| Weights and maturity | Annual stock and catch weights, constant maturity ogive, and spawning fractions of 0.168 for F and 0.25 for M. | Retain supplied biology through age 5, without pre-spawning survival adjustments. | Age-5 biological values proxy for the translated 5+ group; older-age biological contributions and spawning-time survival are approximated. |",
  "| SSB and comparison | SAM calculates SSB using its full 0–8+ population and survival to spawning. | Calculate start-year SSB over the reduced 0–5+ model and use a common definition when comparing with official states. | Official abundance must be aggregated consistently into 5+ for comparison. Native SAM SSB and Fbar are not strictly equivalent to reduced-age quantities. |",
  "",
  "The age truncation is an intentional approximation to align the model plus group with the earliest survey plus group. The translation helper sums source catch and index observations above the selected maximum age; it does not automatically construct abundance-weighted biology or fitted-state starting values for the new plus group.",
  "",
  "The N process has been changed to iid and the F process to cor_rw to more closely approximate the official SAM configuration. The official hockey-stick recruitment relationship and age-correlated survey observation errors are not reproduced.",
  "",
  "The IBTS/BITS Q1 observation precision factors are retained through their equivalent SD multipliers. Missing relative SD factors for other indices use a neutral multiplier of one."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_block),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + sd_block, sd_supplied = ~relative_sd),
    silent = silent
  )
}
