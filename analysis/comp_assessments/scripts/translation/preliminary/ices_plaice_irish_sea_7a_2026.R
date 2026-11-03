## Observations ----
years <- 1981:2025
ages <- 1:8
surveys <- "UK BTS-Q3"

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)
obs$index$q_key <- factor(obs$index$age, levels = 1:7)

N_settings <- list(process = "rw", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "cor_rw", mu_form = NULL, mean_ages = 3:6)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~1)

age_plus_group <- 8L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  F_bar = list(
    status = "non_equivalent",
    definition = "SAM reports the arithmetic mean of F at ages 3-6; tinyAM reports abundance-weighted F.",
    reason = "The source Fbar and tinyAM F_bar use different averaging weights."
  ),
  abundance = list(
    status = "matched",
    definition = "Total abundance is the sum of N at ages 1-8+, in thousand fish.",
    reason = ""
  ),
  ssb = list(
    status = "matched",
    definition = "Start-of-year SSB uses the same stock weights, maturity, and zero pre-spawning mortality.",
    reason = ""
  )
)

background <- c(
  "### Irish Sea plaice: Official 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "The Official SAM assessment fits 1981–2025 at ages 1–8+, with age 8 as the plus group. The native fitted object supplies N-at-age, F-at-age, fixed M, SSB, recruitment, and other fitted summaries. tinyAM uses the same historical years, ages, fixed M, weights, and maturity.",
  "",
  "| Component | Official assessment | tinyAM translation | Difference |",
  "|---|---|---|---|",
  "| Catch | One commercial dead-catch series at age; two year-age cells are missing. Discard survival is 40%. | Keep the fitted catch series, including missing cells. | The cached input data are pre-adjustment; canonical catch is the final fitted series after the TAF survival adjustment. |",
  "| Surveys | Age-specific UK BTS Q3 index at ages 1–7; aggregate NIGFS-WIBTS Q1/Q4 SSB indices. | Fit UK BTS Q3 with q by age. | The two NIGFS indices have no age composition, so they are omitted from tinyAM. The UK BTS scale multiplier is unresolved. |",
  "| N, F, and M | SAM uses a random-walk N process, age-correlated random-walk F, age-7/8 shared F state, and fixed age-specific M. | Use tinyAM random walks for N and recruitment, correlated-RW F, and the same fixed M. | tinyAM approximates SAM's process covariance and state sharing. Fixed M is an input, not an estimated match. |",
  "| Biology and SSB | Fixed maturity, discard-adjusted stock weights, and zero pre-spawning F and M. | Use the same biological surfaces and SSB timing. | No biological values are substituted. |",
  "| Comparisons | The fitted object supplies N/F at age and aggregate SSB, recruitment, total biomass, and Fbar. | Compare N/F at age, abundance, SSB, and recruitment. | Native arithmetic Fbar is shown as non-equivalent to tinyAM's abundance-weighted F_bar. |",
  "",
  "The Q1 and Q4 aggregate biomass indices remain in the canonical database, but tinyAM's current index likelihood requires age-specific values and no age allocation is inferred."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = N_settings,
    rec_settings = rec_settings,
    F_settings = F_settings,
    M_settings = M_settings,
    catch_settings = catch_settings,
    index_settings = index_settings,
    silent = silent
  )
}
