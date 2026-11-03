## Observations ----
years <- 2002:2025
ages <- 1:7
surveys <- c("Q1IBTS+BITS1-7+", "Q4IBTS+BITS1-7+")

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

q1 <- surveys[[1L]]
obs$index$q_block <- ifelse(
  obs$index$survey == q1,
  ifelse(obs$index$age <= 2, as.character(obs$index$age),
    ifelse(obs$index$age <= 5, "3-5", ifelse(obs$index$age == 6, "6", "7+"))
  ),
  ifelse(obs$index$age <= 3, as.character(obs$index$age), "4+")
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_block,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- ifelse(
  obs$index$survey == q1,
  "Q1_all_ages",
  ifelse(obs$index$age == 1, "Q4_age1",
    ifelse(obs$index$age <= 5, "Q4_ages2_5", "Q4_ages6_7")
  )
)
obs$index$sd_block <- factor(obs$index$sd_block)
obs$catch$sd_block <- factor(ifelse(obs$catch$age == 1, "catch_age1", "catch_ages2_7"))

## Model and comparison ----
N_settings <- list(process = "iid", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "ar1", mu_form = ~1, mean_ages = 3:5)
M_settings <- list(
  process = "off", mu_form = NULL,
  mu_supplied = ~M_assumption
)
catch_settings <- list(sd_form = ~0 + sd_block)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)

dat <- tinyAM::prepare_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = N_settings,
  rec_settings = rec_settings,
  F_settings = F_settings,
  M_settings = M_settings,
  catch_settings = catch_settings,
  index_settings = index_settings
)

age_plus_group <- 7L
comparison_aggregates <- c("ssb", "abundance", "recruitment")
comparison_definitions <- list(
  recruitment = list(
    status = "matched",
    definition = "Age-1 recruitment in thousand fish.",
    reason = "Both series use age 1 on the same model-year scale."
  ),
  ssb = list(
    status = "approximate",
    definition = "Published median SSB and tinyAM spawning biomass for ages 1-7+.",
    reason = "SAM estimates an annual stock-weight process; tinyAM treats source weights as fixed."
  ),
  abundance = list(
    status = "matched",
    definition = "Sum of N-at-age over ages 1-7+ in thousand fish.",
    reason = "Official total abundance is derived from the published age-specific estimates."
  ),
  F = list(
    status = "approximate",
    definition = "Fishing mortality at ages 1-7+.",
    reason = "The age-specific quantity matches, but tinyAM's AR(1) process does not reproduce the SAM age-correlated random walk."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Official Fbar over ages 3-5.",
    reason = "The report uses an unweighted age average; tinyAM reports an abundance-weighted mean over the selected ages."
  ),
  M = list(
    status = "non_equivalent",
    definition = "Fixed age- and year-varying M input.",
    reason = "This is an input-consistency check, not an estimated-M comparison."
  ),
  biomass = list(
    status = "approximate",
    definition = "Total biomass over ages 1-7+.",
    reason = "SAM estimates a stock-weight process while tinyAM uses the supplied annual weights as fixed values."
  )
)

background <- c(
  "### Baltic plaice in ICES Subdivisions 21-32: WGBFAS 2026",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM uses years 2002-2025, ages 1-7+, with age 7 as the plus group. | Use the same grid and preserve the Q1-only observation in 2025. | None in the requested observation grid. |",
  "| Catch and indices | Total catch-at-age through 2024; relative Q1 and Q4 age-specific survey indices. | Retain catch and both indices on their native scales and within-year timing. | Observation likelihoods and residual structures are simplified. |",
  "| Recruitment and N | Recruitment random walk; N variance separate at age 1 and shared for ages 2-7. | Use a random-walk recruitment process and IID cohort deviations with exponential N0. | SAM variance sharing is not reproduced exactly. |",
  "| F | Age-correlated random-walk F, ages 1-5 separate and ages 6-7 shared. Report Fbar uses ages 3-5. | Use an AR(1) F process and report over ages 3-5. | tinyAM does not reproduce the SAM random-walk accumulation or full age-correlation structure. |",
  "| M and biology | Fixed time-varying Gislason M; maturity smoothed over five-year windows; annual stock weights enter an estimated process. | Use source M, maturity and weight surfaces as fixed inputs. | Stock-weight uncertainty is omitted. |",
  "| Reported outputs | WGBFAS tables provide N/F at age and aggregate recruitment, SSB and Fbar. | Compare the same age and year metrics; derive total abundance from N-at-age. | Age-specific uncertainty is unavailable. The SAM object's recent point estimates and Fbar range differ from report tables; see the source review. |",
  "",
  "The stock was formed by merging the former 21-23 and 24-32 plaice units. Source observations are not adjusted. Reported table estimates are used for the Official reference; fitted predictions are labeled as coming from the public SAM object."
)

fit <- NULL
if (do_fit) {
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
