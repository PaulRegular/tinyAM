## Observations ----
years <- 1950:2026
ages <- 3:13
surveys <- c(
  "RU-BTr-Q4", "BS-NoRU-Q1(Aco)",
  "BS-NoRu-Q1(BTr)", "FLT007:Eco-NoRu-Q3(Btr)"
)

translation_inputs <- source$inputs[
  !(source$inputs$type == "index" & source$inputs$measure == "index_by_age_group"), ,
  drop = FALSE
]
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

obs$catch$catch_sd_group <- factor(ifelse(obs$catch$age == 3, "age3", "age4plus"))
obs$pop_covariates$F_sd_group <- factor(ifelse(
  obs$pop_covariates$age == 3, "age3", "age4plus"
))
obs$index$q_age_group <- ifelse(obs$index$age == 3, "age3", ifelse(
  obs$index$survey == "RU-BTr-Q4", "age4-8",
  ifelse(obs$index$survey == "BS-NoRU-Q1(Aco)" & obs$index$age <= 7, "age4-7",
  ifelse(obs$index$survey == "BS-NoRU-Q1(Aco)", "age8-9",
  ifelse(obs$index$survey == "BS-NoRu-Q1(BTr)" & obs$index$age <= 7, "age4-7",
  ifelse(obs$index$survey == "BS-NoRu-Q1(BTr)", "age8-10", "age4-9")))))
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age_group,
  drop = TRUE, sep = "."
)

comparison_aggregates <- "ssb"
comparison_definitions <- list(
  F_bar = list(
    status = "approximate",
    definition = "N-weighted mean fishing mortality at ages 4-7, recalculated from official F-at-age and N-at-age.",
    reason = "The published Fbar is for ages 4-7; the common comparison applies tinyAM's population-weighted definition."
  ),
  ssb = list(
    status = "matched",
    definition = "Beginning-of-year SSB from N-at-age, stock weight, maturity and zero pre-spawning mortality.",
    reason = "The official source surfaces and spawning fractions provide the same biological definition."
  )
)
age_plus_group <- 13

## Model ----
background <- c(
  "### Northeast Arctic haddock: official 2026 JRN-AFWG SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Survey observations extend through 2026; catch ends in 2025. Population ages are 3-13+, with age 13 as the plus group. | Fit 1950-2026 at ages 3-13+, retaining missing 2026 catch and 2026 survey observations. | The source population grid and the latest Q1 survey data are retained. |",
  "| N and recruitment | Age-3 recruitment follows a random walk; SAM has age-keyed abundance-state variances and its own initial-state integration. | Random-walk recruitment and IID cohort-survival deviations with exponential initial abundance. | tinyAM uses a simpler initial-state treatment and one shared SD for older abundance states. |",
  "| F | F states are separate at ages 3-7; ages 8-13+ share the age-8 state. Annual F variation has AR(1) correlation across ages; age 3 has a separate process SD. | Correlated F random walk with AR(1) cross-age increments, independent F states by age, and separate age-3 versus age-4+ increment SDs. | tinyAM does not share the official ages 8-13+ F state. |",
  "| M | Fixed at 0.2 plus reconstructed cod-predation mortality at ages 3-6; ages 7-13+ remain 0.2. | Use the final source M-at-age as fixed M. | Numerical fixed M is retained; no M process is estimated. |",
  "| Catch and indices | Residual catch-at-age is the only catch likelihood. Four lognormal survey series use quarter-specific timing, age-shared q and observation SD groups; survey residuals are AR(1)-correlated across age. | Retain the residual catch and four indices, with matching sampling times, q sharing groups, and catch/index SD groups. | Native q-power density dependence and survey residual correlation are not represented. Grouped terminal index ages are retained in the database but omitted as single-age tinyAM observations. |",
  "| Biology | Annual stock weight, maturity and fixed M are supplied; pre-spawning F and M fractions are zero. | Retain those biological surfaces and the beginning-of-year spawning convention. | The source biology and SSB calculation are preserved. |",
  "| Comparisons | N-at-age, recruitment, SSB and total biomass extend through 2026; F-at-age, M-at-age and Fbar ages 4-7 end in 2025. Only aggregate estimates have Low/High bounds, with coverage unspecified. | Compare the available outputs over 1950-2026; derive total abundance from reported N-at-age. | Fixed-M agreement is by construction; age-specific uncertainty is unavailable. |",
  "",
  "Survey observations remain on their native numerical scales because the physical index multipliers are not stated. The 69 missing SAM survey cells are not observations. This fit uses normal tinyAM starting values and supplies no explicit `start_par`."
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
    F_settings = list(
      process = "cor_rw", mu_form = NULL,
      sd_form = ~0 + F_sd_group, mean_ages = 4:7
    ),
    M_settings = list(
      process = "off", mu_form = NULL, mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~0 + catch_sd_group),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + survey
    ),
    silent = silent
  )
}
