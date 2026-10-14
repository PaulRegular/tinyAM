## Observations ----
years <- 1947:2025
ages <- 0:8
surveys <- c("HERAS", "IBTS-Q1", "IBTS0", "IBTS-Q3")

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  index_weight_source = "source",
  surveys = surveys,
  assumptions = source$assumptions
)
obs$pop_covariates$M_assumption <- obs$pop_covariates$M_assumption + 0.02

obs$catch$sd_block <- cut_ages(obs$catch$age, c(0, 2, 7, 8))
levels(obs$catch$sd_block) <- c("ages_0_1", "ages_2_6", "ages_7_8", "ages_7_8")
heras <- obs$index$survey == "HERAS"
q3 <- obs$index$survey == "IBTS-Q3"
heras_sd <- cut_ages(obs$index$age[heras], c(1:4, 7, 8))
levels(heras_sd) <- c("1", "2", "3", "4_6", "7_8", "7_8")
q3_sd <- cut_ages(obs$index$age[q3], c(0, 2, 5))
levels(q3_sd) <- c("0_1", "2_5")
obs$index$sd_block <- obs$index$survey
obs$index$sd_block[heras] <- paste0("HERAS_", heras_sd)
obs$index$sd_block[q3] <- paste0("IBTS-Q3_", q3_sd)
obs$index$sd_block <- factor(obs$index$sd_block)
heras_q <- cut_ages(obs$index$age[heras], c(1, 3, 8))
levels(heras_q) <- c("1_2", "3_8")
obs$index$q_key <- obs$index$survey
obs$index$q_key[heras] <- paste0("HERAS_", heras_q)
obs$index$q_key[q3] <- paste0("IBTS-Q3_", obs$index$age[q3])
obs$index$q_key <- factor(obs$index$q_key)


model_dat <- tinyAM::prepare_tam(
  data = obs, years = years, ages = ages, N_settings = list(
    process = "rw",
    init = "free"
  ), F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6), M_settings = list(
    process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption
  ), catch_settings = list(sd_form = ~ 0 + sd_block),
  index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + sd_block)
)
start_par <- tinyAM::make_par(model_dat)
source_N <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages, 1000)
source_F <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age", years, ages)
start_par$log_r0 <- log(source_N[1, 1])
start_par$log_n0 <- log(source_N[1, -1])
start_par$log_r <- log(source_N[-1, 1])
start_par$log_n <- log(source_N[-1, -1, drop = FALSE])
start_par$log_f <- log(source_F)


## Background and comparisons ----

age_plus_group <- 8

comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, F_bar = 1)


background <- c(
  "### North Sea autumn-spawning herring: accepted 2026 single-fleet assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",

  "| Years | The population model spans 1947–2026. Catch, weights, maturity and published M inputs end in 2025. | Fit 1947–2025 and compare over those common years. | Complete biological input tables are unavailable for 2026. |",

  "| Ages | Ages 0–8 winter rings, with age 8 as the plus group. | Retain ages 0–8+, with recruitment at age 0. | Source ages represent winter-ring classes rather than conventional chronological ages. |",

  "| N | Recruitment has a separate process variance; ages 1–8 share a survival-process variance. | Random-walk recruitment, IID cohort-survival deviations with one shared SD, and free initial older-age abundance. | The independent survival deviations and shared older-age SD broadly approximate SAM, but its initial-state treatment and full process specification are not reproduced. |",

  "| F | Random-walk F with AR(1) correlation across ages in annual increments. Ages 0–6 have separate F states, ages 7–8 share one state, and innovation variances are grouped over ages 0–1, 2–5 and 6–8. Fbar uses ages 2–6. | Correlated F random walks (cor_rw), with AR(1) correlation between age-specific annual log-F increments, one shared process SD and Fbar ages 2–6. | Cross-age correlation is now represented, but the source variance groups and shared oldest-age F state are not. |",

  "| M | Fixed annual SMS-2023 age surface with an additional 0.02 in the accepted model configuration. | Supply the translated 1947–2025 M surface as fixed numerical values. | No M process is estimated; the accepted final-M convention must be preserved without applying the offset twice. |",

  "| Catch | Direct catch-at-age in thousand fish; fishery-closure years 1978–1979 have missing observations. | Convert to individual fish and fit positive catch observations with independent lognormal errors and SD groups defined by sd_block. | Source catch-error structure may differ; zero observations cannot enter the standard lognormal likelihood directly. |",

  "| Indices | HERAS, IBTS-Q1, IBTS0 and IBTS-Q3 are age-structured. Four larval indices use partial spawning-component likelihoods. | Retain the four age-structured surveys, with source-derived q and observation-SD groups. Omit the larval indices. | Partial spawning-component indices are not abundance-at-age observations. IBTS-Q3 cross-age observation correlation is omitted, and timing remains approximate. |",

  "| Weights and maturity | Annual stock and catch weights and annual maturity are tabulated through 2025. | Retain the source annual biological surfaces for 1947–2025. | No alternative biological surfaces are substituted. |",

  "| SSB | The accepted assessment uses its spawning-time biomass convention. | Calculate SSB from translated stock weights, maturity and fitted abundance using tinyAM's timing convention. | The source spawning fractions remain unavailable, so native SSB definitions may differ. |",

  "",
  "The revised translation uses IID cohort-survival deviations and correlated F random walks in place of the previous random-walk N process and independent age-specific F random walks. This more closely approximates the accepted assessment's broad process structure while retaining a simpler parameterization.",
  "",
  "The accepted FLSAM shares its oldest F state and estimates three groups of F-process variances; tinyAM retains separate age-specific F states and uses one shared increment SD. Recruitment remains a separate random-walk process.",
  "",
  "Accepted N and F surfaces provide optimizer starting values only; tinyAM re-estimates the population and mortality states. The historical fit is restricted to 1947–2025, although accepted N and F estimates extend through 2026.",
  "",
  "Survey sampling times are approximate seasonal midpoints because the accepted fleet.txt was unavailable. Survey indices retain their published native scales, with estimated catchability absorbing unresolved units.",
  "",
  "The four larval abundance indices and their partial spawning-component likelihoods are omitted. Comparisons therefore assess whether the simplified age-structured model can approximate the accepted historical stock dynamics without reproducing the full observation likelihood."
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
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_block),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 0 + sd_block),
    silent = silent,
    start_par = start_par
  )
}
