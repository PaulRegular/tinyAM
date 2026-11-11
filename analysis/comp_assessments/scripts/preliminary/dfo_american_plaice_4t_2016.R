## Observations ----
years <- 1976:2015
ages <- 4:16
age_plus_group <- 16
sampling_times <- c(
  "DFO Research Vessel Trawl Survey" = 0.75,
  "Mobile Gear Sentinel Survey" = 0.625
)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  surveys = names(sampling_times),
  sampling_times = sampling_times,
  assumptions = source$assumptions
)

# The source rows carry a shared 0.75 approximation. The Sentinel survey ran in August.
obs$index$samp_time[obs$index$survey == "Mobile Gear Sentinel Survey"] <- 0.625

obs$pop_covariates$M_assumption <- ifelse(
  obs$pop_covariates$age <= 9, 0.6, 0.3
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$age,
  drop = TRUE, sep = "."
)

comparison_aggregates <- c("abundance", "ssb")
comparison_age_groups <- list(
  F = list(`4-9` = 4:9, `10+` = 10:16),
  M = list(`4-9` = 4:9, `10+` = 10:16)
)
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total January 1 abundance summed over ages 4-16+.",
    reason = "The source total is derived from the reported marginal median N-at-age values; it is not a posterior median of total abundance."
  ),
  recruitment = list(
    status = "matched",
    definition = "January 1 abundance at age 4 in the same calendar year.",
    reason = "Both series use age-4 abundance as recruitment; the source recruitment values are derived from the reported N-at-age medians."
  ),
  ssb = list(
    status = "matched",
    definition = "Beginning-of-year spawning biomass using January 1 N, weight-at-age and maturity-at-age.",
    reason = "The reference uses reported median N-at-age and the same weights and maturity as tinyAM. Summing these contributions gives 59,928 t in 2015. The VPA's reported posterior median SSB is 55,023 t."
  )
)

background <- c(
  "### 4T American plaice: 2016 DFO VPA",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | VPA covers 1976-2015, ages 4-16+, with age 16 as the plus group. | Fit the same years and ages. | The source model grid is retained. |",
  "| N and recruitment | The VPA estimates 12 terminal-age abundance parameters and reconstructs population numbers from catch and M. Recruitment is age 4. | Free initial abundance at age, with an age-4 recruitment random walk. | tinyAM estimates abundance forward and does not reproduce the source VPA recursion or its terminal-state parameterization. |",
  "| F | Fishing mortality is recovered within the VPA from catch and population dynamics. | Age-specific temporal random walks. | tinyAM uses an explicit F process instead of the source VPA catch equation. |",
  "| M | Two log-scale random walks: ages 4-9 and 10+, with fixed increment SD 0.075 and initial-M prior means 0.6 and 0.3. | Fix the first year at 0.6 and 0.3 per year, then fit two random walks starting in 1977 with an estimated common increment SD. | The 1977 process states are free of a penalty tying them to the 1976 rates. |",
  "| Catch and indices | Commercial catch-at-age, including estimated discards, is a VPA input. The RV and mobile sentinel indices use age-specific q and lognormal observation errors. | Catch-at-age is fitted with tinyAM's lognormal likelihood; both indices have separate age-specific q and observation SDs. | tinyAM's likelihoods and population equations differ from the VPA. |",
  "| Survey timing | The RV survey occurs in September and the mobile Sentinel survey in August. The 2003 RV observations were excluded because an uncalibrated vessel was used. | Retain the source index years. Use 0.75 for RV and 0.625 for the August Sentinel midpoint. | The fractions approximate annual survey timing. |",
  "| Biology and SSB | Annual population weights are adjusted to January 1; maturity is age- and year-specific. Figure 23 reports annual VPA SSB. | Retain both source surfaces and use them with source median N-at-age to calculate the SSB comparison. | In 2015 the comparison is 59,928 t and the reported VPA posterior median is 55,023 t. The 2003 age-16+ weight of 0.000 remains unresolved. |",
  "| Outputs | N-at-age is tabulated. Annual F, M and SSB medians appear in figures; 2015 SSB is also reported as 55,023 t. | Compare N-at-age, derived total abundance, age-4 recruitment, common-definition SSB, and digitized grouped F/M where recoverable. | Figure values are approximate; some grouped mortality points cannot be separated, and F uncertainty bands are not visible. |",
  "",
  "The fit uses package starting values. F and M age-group summaries use the source's 4-9 and 10+ group definitions."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "rw", mu_form = NULL),
    M_settings = list(
      process = "rw", mu_form = NULL, mu_supplied = ~M_assumption,
      age_breaks = c(4, 10, 16)
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + q_key),
    silent = silent
  )
}
