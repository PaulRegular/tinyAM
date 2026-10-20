## Observations ----
years <- 1985:2015
ages <- 1:11
age_plus_group <- 11

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  sampling_times = c("DFO Summer RV" = 0.5),
  assumptions = source$assumptions
)

obs$index <- obs$index[
  obs$index$year %in% years & obs$index$age <= 10,
  , drop = FALSE
]
obs$index$q_age <- factor(obs$index$age, levels = 1:10)

## Background and comparisons ----
comparison_outputs <- source$outputs
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Age-4+ biomass from the VPA compared with tinyAM SSB using the framework's static spawning weights and maturity.",
    reason = "The source reports annual age-4+ biomass, while tinyAM uses time-invariant per-recruit biological values because annual stock weights are unavailable."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-1 abundance in the same calendar year.",
    reason = "Recruitment age matches, but tinyAM's initial-state and process assumptions differ from the VPA."
  )
)

background <- c(
  "### 4X5Y haddock: 2016 DFO ADAPT VPA assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Fit 1985–2015, ages 1–11+; population estimates through 2016. | Fit 1985–2015, ages 1–11+, with age 11 as the plus group. | The source model grid is retained; source catch ages 11–16 are summed to 11+. |",
  "| N and recruitment | ADAPT VPA estimates population abundance; age-1 N in 2016 is assigned from a geometric mean excluding the 2013 year class. | Exponential initial abundance with no separate N process; recruitment is age 1. | tinyAM's initial-state construction does not reproduce the VPA terminal-age parameterization. |",
  "| F | VPA estimates annual F-at-age and defines F11+ through a Fratio. | AR(1) F process with age-specific mean states. | This is a practical tinyAM process approximation, not the VPA estimation method or Fratio. |",
  "| M | Fixed at 0.2, except ages 10–11+ at 0.3, 0.6, and 0.9 in successive five-year blocks from 2000 to 2015. | Same fixed annual age-specific M surface. | Numerical M inputs are preserved. |",
  "| Catch | Commercial catch-at-age is an input to the VPA; ages 1–11+ are used. | Catch-at-age is fitted with tinyAM's lognormal observation likelihood. | tinyAM does not treat catch-at-age as an exact VPA input. |",
  "| Survey and q | DFO Summer RV swept-area abundance-at-age, ages 1–10; age-specific q; zero index cells treated as missing. | Same survey index ages and age-specific q; zero values are omitted by the log-scale observation likelihood. | The annual within-summer timing is not reported; 0.5 is used as a midyear approximation. |",
  "| Biology and SSB | Source SSB is the reported age-4+ biomass. | Static spawning weight and maturity from the framework's per-recruit table are used to calculate tinyAM SSB. | The annual population weight surface is unavailable, so SSB agreement is approximate. |",
  "",
  "The 2024 advice is stored as a separate summary-only record. Its 2022 SSB value and management-procedure results are not mixed with this 2016 VPA series. The one unresolved 2015 age-2 survey cell is omitted from the input table and remains documented in the source review."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~ 0 + factor(age)),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_age, sd_form = ~1),
    silent = silent
  )
}
