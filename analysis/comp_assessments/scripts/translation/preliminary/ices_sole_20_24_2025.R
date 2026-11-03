## Observations ----
years <- 1984:2024
ages <- 1:9
assessment_id <- source$assessment$assessment_id[[1L]]
survey_names <- c("Q4Solesurvey1-9+", "Gillnet", "Trawl")
translation_inputs <- source$inputs[!(
  source$inputs$type == "index" &
    source$inputs$measure == "numbers_at_age" &
    is.na(source$inputs$value)
), , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  surveys = survey_names,
  sampling_times = c(
    "Q4Solesurvey1-9+" = 0.8754645,
    Gillnet = 0.56,
    Trawl = 0.875
  ),
  assumptions = source$assumptions
)

obs$index$q_block <- ifelse(
  obs$index$survey == "Q4Solesurvey1-9+",
  ifelse(obs$index$age >= 4, "4+", as.character(obs$index$age)),
  ifelse(obs$index$age >= 5, "5+", as.character(obs$index$age))
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_block,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- factor(obs$index$survey, levels = survey_names)
obs$maturity$spawn_time <- 0

## Model and comparison ----
age_plus_group <- 9L
comparison_outputs <- source$outputs
comparison_definitions <- list(
  recruitment = list(
    status = "matched",
    definition = "Age-1 recruitment in the same model year.",
    reason = "The source SAM and tinyAM model both use age 1 as the first modeled age."
  ),
  ssb = list(
    status = "approximate",
    definition = "Official SSB and tinyAM mature biomass for ages 1-9+.",
    reason = "SAM estimates a stock-weight process; tinyAM uses the annual reported stock weights as fixed inputs."
  ),
  abundance = list(
    status = "matched",
    definition = "Sum of abundance at ages 1-9+ in each year.",
    reason = "Official total abundance is derived from the native N-at-age surface and has no separate uncertainty estimate."
  ),
  F = list(
    status = "matched",
    definition = "Instantaneous fishing mortality by age and year.",
    reason = "SAM does not estimate age-1 F; that age remains unavailable in the official output."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Official Fbar for ages 4-8 compared with tinyAM Fbar for the same ages.",
    reason = "The source report gives the age range but not enough detail to confirm that its aggregation weights match tinyAM's abundance-weighted Fbar."
  ),
  biomass = list(
    status = "approximate",
    definition = "Official total biomass for ages 1-9+.",
    reason = "SAM estimates stock weights; tinyAM uses the observed stock-weight inputs as fixed values."
  ),
  M = list(
    status = "non_equivalent",
    definition = "Natural mortality is fixed at 0.1 per year in both models.",
    reason = "This checks consistency of a shared input and is not an estimated-M comparison."
  ),
  M_bar = list(
    status = "non_equivalent",
    definition = "Mbar is calculated from the shared fixed M input.",
    reason = "This is a fixed-input consistency check, not a comparison of estimated natural mortality."
  )
)

background <- c(
  "### Common sole in ICES Subdivisions 20–24: WGBFAS 2025",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1984–2024 at ages 1–9+, with age 9 as the plus group. | Fit the same years and ages. | None in the modeled grid. |",
  "| Recruitment and N | Random-walk recruitment; N-process variance is separate at age 1 and shared for ages 2–9. | Use random-walk recruitment, IID cohort survival deviations, and exponential initial abundance. | The N-process age-correlation structure is simplified. |",
  "| F | Random-walk F with compound-symmetry correlation across ages; ages 2–5 have separate F states and ages 6–9 share a state. | Use a stationary AR1 F process across years and ages with one mean level; report Fbar over ages 4–8. | This stable approximation does not reproduce SAM's random-walk accumulation, compound-symmetry correlation, shared older-age state, or initial-state convention. |",
  "| M and biology | M is fixed at 0.1; maturity is 0 through age 2 and 1 from age 3; annual stock and catch weights are inputs to estimated weight processes. | Use the same M, maturity, and reported annual weights, with weights held fixed. | Weight-process uncertainty is omitted. |",
  "| Catch | Landings-only catch-at-age, ages 2–9+, in thousand fish; discards are excluded. | Fit the same catch observations with a shared lognormal observation SD. | SAM's age-grouped catch variance is not reproduced. |",
  "| Indices | Q4 Sole survey (ages 1–9+, 2004–2024), gillnet CPUE (ages 2–8, 1994–2007), and trawl CPUE (ages 2–6, 1987–2008). | Retain all three at their native relative-index scales and source timings. | Index SDs are separated by series; the gillnet AR(1) residual correlation across ages is not represented. |",
  "| Catchability | Q4 ages 1, 2, 3 and 4+ share four q values; gillnet and trawl each use ages 2, 3, 4 and 5+ groups. | Use those same q-sharing groups. | None in the age-sharing pattern. |",
  "",
  "The 2026 update is retained separately as aggregate summary data. Its outputs are not substituted into the 2025 detailed fit or comparison."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "ar1", mu_form = ~1, mean_ages = 4:8),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + sd_block
    ),
    silent = silent
  )
}
