## Observations ----
years <- 1979:2024
ages <- 1:12
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)


## Background and comparisons ----

age_plus_group <- 12

comparison_scales <- c(
  ssb = 0.001, recruitment = 0.001, N = 0.001, abundance = 0.001, biomass = 0.001,
  biomass_at_age = 0.001
)


background <- c(
  "### Icelandic haddock: 2025 benchmark assessment", "", print_sources(source$assessment),
  "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
  "| Years | Catch and biological inputs cover 1979–2024; population estimates extend to 2025. | Fit 1979–2024. | The 2025 estimates lie beyond the detailed input period used in tinyAM. |",
  "| Ages | Ages 1–12+, with age 12 as the plus group. | Fit ages 1–12, retaining age 12 as the plus group. | The accepted age range and plus group are retained. |",
  "| N | SAM estimates annual age-1 recruitment and abundance-at-age with process variation. | Random-walk recruitment, IID cohort-survival deviations and free initial abundance. | The source's exact process covariance, variance sharing and initial-state treatment are not reproduced. |",
  "| F | Fishing mortality and selectivity vary over time under the 2025 benchmark SAM model. | Correlated random-walk F (cor_rw), allowing annual log-F increments to be correlated across ages. | This replaces the previous stationary AR1 approximation, but the accepted covariance and parameter-sharing structure are not fully verified or reproduced. |",
  "| M | Natural mortality is fixed at 0.2 per year for all ages. | Supply fixed M = 0.2 for every age and year. | The accepted assumption is retained. |",
  "| Catch | Commercial catch numbers-at-age and catch weights are supplied. | Retain catch-at-age with one estimated lognormal observation SD. | tinyAM does not reproduce the complete accepted catch observation likelihood. |",
  "| Indices | IS-SMB in March and IS-SMH in October provide age-specific indices. | Retain both surveys, with survey-age catchability and survey-specific observation SDs. | Index units and exact native sampling fractions remain unresolved; observation covariance may differ. |",
  "| Weights and maturity | Stock weights and maturity come from the March survey; catch weights come from commercial samples. Pre-1985 stock-weight and maturity inputs use 1985 values. | Retain the supplied annual biological surfaces and early-year substitutions. | No alternative biological inputs are introduced. |",
  "| SSB | Spawning biomass accounts for fractions 0.4 of F and 0.3 of M before spawning. | Calculate beginning-of-year mature biomass from fitted N, weights and maturity. | The comparison uses a common-definition mature biomass rather than the accepted native spawning-time SSB. |",
  "",
  "The revised translation uses IID cohort-survival deviations and correlated F random walks, replacing the previous stationary AR1 F process. Recruitment remains a separate random walk under the default tinyAM settings. The revised F process permits persistent temporal changes in age-specific fishing mortality, with correlation between annual increments at neighbouring ages.",
  "",
  "For tinyAM, survey timings of 0.20 for March and 0.80 for October are month-based approximations. The accepted model's exact sampling fractions were not recovered. Both surveys retain their published native index scales.",
  "",
  "The source's native SSB includes pre-spawning mortality fractions of 0.4 for F and 0.3 for M. tinyAM does not apply these fractions in this translation. The comparison therefore uses common-definition mature biomass calculated from accepted N and the translated weights and maturity.",
  "",
  "The accepted reference biomass is defined by fish length (45 cm and larger). This cannot be reconstructed exactly from the available age-based tinyAM inputs and is not used as a directly comparable total-biomass metric."
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
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + survey),
    silent = silent
  )
}
