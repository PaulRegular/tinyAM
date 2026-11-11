## Observations ----

assessment_id <- source$assessment$assessment_id[[1]]
years <- 1988:2025
ages <- 1:8
survey <- "EU Flemish Cap bottom-trawl survey"

inputs <- source$inputs
inputs$value <- suppressWarnings(as.numeric(as.character(inputs$value)))
inputs <- inputs[!(inputs$type == "index" & inputs$measure != "numbers_at_age"), , drop = FALSE]
inputs <- inputs[!(inputs$type %in% c("catch", "index") &
  (is.na(inputs$value) | inputs$value <= 0)), , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = survey,
  sampling_times = c("EU Flemish Cap bottom-trawl survey" = 0.54),
  exclude_index_years = list("EU Flemish Cap bottom-trawl survey" = 2018),
  assumptions = source$assumptions
)

obs$index$q_key <- factor(
  paste0("age", ifelse(obs$index$age >= 4, "4plus", obs$index$age)),
  levels = c("age1", "age2", "age3", "age4plus")
)
obs$index$q_power_group <- ifelse(obs$index$age == 1, "age1", NA_character_)
# Keep a valid formula row at age 1, where there is no F state.
obs$catch$F_age_group <- factor(
  ifelse(obs$catch$age >= 7, "age7plus", paste0("age", pmax(2L, obs$catch$age))),
  levels = c("age2", "age3", "age4", "age5", "age6", "age7plus")
)

translation <- attr(obs, "translation")
translation$M$status <- "source_estimated_fixed_at_prior_median"
translation$M$notes <- paste(
  "The source estimates M with age-specific priors; this approximation fixes",
  "tinyAM M at the reported prior medians. The posterior M estimates are not",
  "treated as observations or displayed as a direct M comparison."
)
attr(obs, "translation") <- translation

age_plus_group <- 8L
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_outputs <- source$outputs[!(
  source$outputs$type == "mortality" &
    source$outputs$measure == "natural_mortality_at_age"
), , drop = FALSE]
comparison_definitions <- list(
  F_bar = list(
    status = "approximate",
    definition = "Official mean F over ages 3–5 compared with tinyAM N-weighted mean F over ages 3–5.",
    reason = "Both summarize ages 3–5, but tinyAM weights age-specific F by abundance."
  ),
  abundance = list(
    status = "matched",
    definition = "Common-definition total abundance over ages 1–8+.",
    reason = "Both sum abundance over the same modeled ages."
  ),
  biomass = list(
    status = "matched",
    definition = "Common-definition stock biomass over ages 1–8+.",
    reason = "The comparison uses source N and the same translated stock weights over the fitted age range."
  ),
  ssb = list(
    status = "matched",
    definition = "Common-definition mature biomass over ages 1–8+ at start of year.",
    reason = "The comparison applies the translated stock weights and maturity to source N."
  ),
  recruitment = list(
    status = "matched",
    definition = "Beginning-of-year abundance at age 1.",
    reason = "Both models define recruitment as age-1 abundance in the same calendar year."
  )
)

## Model ----

background <- c(
  "### NAFO Division 3M cod and the 2026 assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1988–2025, ages 1–8+. | 1988–2025, ages 1–8+. | Age 8 is a single terminal group because older age-specific biological inputs are unavailable. |",
  "| Catch | Catch-at-age in thousand fish; zero values at ages 2+ are treated as missing. | Catch likelihood uses positive age-2–8+ values and F at age 1 is zero. | The source catch-at-age CV of 0.2 gives log-scale SD 0.198. The source also models annual total catch weight separately. |",
  "| Index and q | EU Flemish Cap age-specific indices; q groups are ages 1, 2, 3 and 4+. | Use all four q groups and estimate an abundance exponent only at age 1. | tinyAM constrains the exponent to be positive and omits the source prior. The translation retains the 2018 survey exclusion because one missing source age lies inside 8+. |",
  "| N and recruitment | Cohort projection, age-1 recruitment with a broad prior, and a first-year age profile. | No additional N process, free initial ages 2+, and a recruitment random walk at age 1. | The source uses Bayesian priors for recruitment and the initial age profile. |",
  "| F | Annual scalar F follows AR(1); age-specific selectivity follows random walks, age 1 is fixed at zero, and F at 8+ equals age 7. | Age-blocked AR(1) F surface with ages 7 and 8 shared and structural zero at age 1. | tinyAM does not reproduce the separate scalar and selectivity processes; it uses one state per retained age block. |",
  "| M | Age-specific M is estimated with lognormal prior medians and CV 0.15. | Hold the reported prior-median age vector fixed. | This keeps an explicit source-based baseline but does not estimate M or reproduce its prior. Posterior M is not compared to the fixed tinyAM values. |",
  "| SSB and comparisons | Annual N, F, M, recruitment, total abundance and biomass outputs are reported. | Compare N/F at ages 1–8+, age-1 recruitment, and common-definition abundance, biomass and SSB where available. | Fbar and M use different definitions or estimation treatments. |",
  "",
  "The 0.54 sampling time is the midpoint of the reported 0.50–0.58 timing interval, used as a single-time approximation. The source model integrates over a timing interval.",
  "",
  "Official N and positive F estimates supply starting values. The fit estimates those states from the translated observations."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"

  # Official N and F estimates supply starting values only.
  source_N <- .translation_start_surface(
    source$outputs, "population", "numbers_at_age", years, ages
  )
  F_block_ages <- 2:7
  source_F <- .translation_surface(
    subset(source$outputs, type == "mortality" & measure == "fishing_mortality_at_age"),
    years, F_block_ages, "Official fishing mortality"
  )
  log_F_start <- matrix(0, length(years), length(F_block_ages),
    dimnames = list(year = years, age_block = c(2:6, "7-8"))
  )
  positive_F <- is.finite(source_F$obs) & source_F$obs > 0
  log_F_start[positive_F] <- log(source_F$obs[positive_F])

  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(
      process = "ar1", mu_form = ~ 0 + F_age_group,
      age_breaks = c(F_block_ages, 9), mean_ages = 3:5
    ),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~0, sd_supplied = ~ I(0.1980422)),
    index_settings = list(
      q_form = ~ 0 + q_key + q_power(N, by = q_power_group),
      sd_form = ~0, sd_supplied = ~ I(0.2935604)
    ),
    start_par = list(
      log_r0 = log(source_N[1, 1]),
      log_n0 = log(source_N[1, -1]),
      log_r = log(source_N[-1, 1]),
      log_f = log_F_start
    ),
    silent = silent
  )
}
