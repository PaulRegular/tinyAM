## Observations ----
years <- 1979:2026
ages <- 1:12
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)


## Background and comparisons ----
age_plus_group <- 12L

background <- c(
  "### Icelandic haddock: 2026 assessment update", "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | Catch and catch weights through 2025; spring survey, stock biology and abundance through 2026; F-at-age through 2025. | Fit 1979–2026. | The 2026 year has spring survey observations but no catch or reported F-at-age. |",
  "| Ages | Ages 1–12+, with age 12 as the plus group. | Fit ages 1–12, retaining age 12 as the plus group. | The source age structure is retained. |",
  "| N | SAM reports abundance-at-age through 2026. | Random-walk recruitment, IID cohort-survival deviations and free initial abundance. | The official process covariance, variance sharing and initial-state integration are not fully published. |",
  "| F | Fishing mortality and selectivity vary over time; F-at-age is reported through 2025. | Correlated random-walk F (cor_rw), allowing adjacent-age increments to correlate. | This is a flexible approximation; the precise official process covariance and parameter sharing are not reconstructed. |",
  "| M | Natural mortality is fixed at 0.2 per year for all ages. | Supply fixed M = 0.2 for every age and year. | The reported assumption is retained. |",
  "| Catch | Commercial catch numbers-at-age and weights are available through 2025. | Retain catch-at-age with one estimated lognormal observation SD. | tinyAM does not reproduce the full SAM catch likelihood. The source lists a uniform 4,000-g catch weight in 2026 without catch-at-age; its projection role is unspecified, so it is excluded from the fitted catch-weight series. |",
  "| Surveys | IS-SMB spring survey (1985–2026) and IS-SMH autumn survey (1995–2025) provide age-specific indices. | Retain both surveys, with survey-age catchability and survey-specific observation SDs. | Survey units are unresolved and timings 0.20 and 0.80 are month-based approximations. |",
  "| Biology | Catch weights, stock weights and maturity are provided; pre-1985 stock weights and maturity use the 1985 vectors. | Retain the published annual surfaces and early-year substitutions. | No alternative biological series are introduced. |",
  "| SSB | The official estimate accounts for pre-spawning fractions 0.4 of F and 0.3 of M. | Compare common-definition mature biomass calculated from fitted N and the same weights and maturity. | This is not the official spawning-time SSB. |",
  "| Reference biomass | Biomass of fish at least 45 cm. | Not directly represented. | tinyAM uses age-based biomass, not length-based biomass. |",
  "",
  "The public table's annual F column does not specify its age range or averaging rule and is not imported as Fbar. The age-specific F surface is retained. Summary low/high values for SSB, recruitment and reference biomass are stored without assigning a confidence level because the table does not specify interval coverage.",
  "",
  "The translated model uses the 2025 benchmark framework described in the 2026 report. The exact SAM likelihood covariance and parameter-sharing settings were not fully recovered. Survey timing remains approximate, and the published survey index scale is retained."
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
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey),
    silent = silent
  )
}
