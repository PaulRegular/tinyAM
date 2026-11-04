## Observations ----

assessment_id <- source$assessment$assessment_id[[1]]
years <- 1988:2025
ages <- 2:8
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

obs$index$q_key <- factor(ifelse(
  obs$index$age == 2, "age2",
  ifelse(obs$index$age == 3, "age3", "age4plus")
), levels = c("age2", "age3", "age4plus"))
obs$catch$F_age_group <- factor(ifelse(
  obs$catch$age >= 7, "age7plus", paste0("age", obs$catch$age)
), levels = c("age2", "age3", "age4", "age5", "age6", "age7plus"))

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
    status = "approximate",
    definition = "Common-definition total abundance over ages 2–8+.",
    reason = "Age 1 is omitted because the official model sets its catch and F to zero; all source and tinyAM totals use the same retained age range."
  ),
  biomass = list(
    status = "approximate",
    definition = "Common-definition stock biomass over ages 2–8+.",
    reason = "The comparison uses source N and the same translated stock weights over the fitted age range."
  ),
  ssb = list(
    status = "approximate",
    definition = "Common-definition mature biomass over ages 2–8+ at start of year.",
    reason = "Age 1 is omitted, and the official timing interval and estimated M treatment are simplified."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Beginning-of-year abundance at age 2 in the same calendar year.",
    reason = "The source defines recruitment at age 1; tinyAM starts at age 2 after omitting the structural age-1 catch group."
  )
)

## Model ----

N_settings <- list(process = "off", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(
  process = "ar1",
  mu_form = ~ 0 + F_age_group,
  age_breaks = c(2, 3, 4, 5, 6, 7, 9),
  mean_ages = 3:5
)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption)
catch_settings <- list(sd_form = ~0, sd_supplied = ~ I(0.1980422))
index_settings <- list(q_form = ~ 0 + q_key, sd_form = ~0, sd_supplied = ~ I(0.2935604))

background <- c(
  "### NAFO Division 3M cod: preliminary translation of the 2026 assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1988–2025, ages 1–8+. | 1988–2025, ages 2–8+. | Age 1 is omitted because the source fixes age-1 catch and F to zero. tinyAM age 8 is a single terminal group; no older age-specific biological inputs are available to expand it. |",
  "| Catch | Catch-at-age in thousand fish; zero values at ages 2+ are treated as missing. | Positive age-2–8+ values; source zeros and missing values are omitted. | tinyAM uses a log-observation likelihood. The source catch CV of 0.2 is represented by a fixed log-scale SD of 0.198. |",
  "| Index and q | EU Flemish Cap age-specific indices; q groups are ages 1, 2, 3 and 4+. | Retain age-2, age-3 and age-4+ q groups; omit zero observations. | Source age-1 q is size-dependent and age 1 is outside the fitted grid. The 2018 survey year is omitted because one missing source age lies inside 8+. |",
  "| N and recruitment | Cohort projection, age-1 recruitment with a broad prior, and a first-year age profile. | No additional N process, free initial ages from age 3+, and a recruitment random walk at age 2. | The age-1 recruitment and initial-age prior are not reproduced. |",
  "| F | Annual scalar F follows AR(1); age-specific selectivity follows random walks, age 1 is fixed at zero, and F at 8+ equals age 7. | Age-blocked AR(1) F surface with ages 7 and 8 shared. | tinyAM does not reproduce the separate scalar and selectivity processes; it uses one state per retained age block. |",
  "| M | Age-specific M is estimated with lognormal prior medians and CV 0.15. | Hold the reported prior-median age vector fixed. | This keeps an explicit source-based baseline but does not estimate M or reproduce its prior. Posterior M is not compared to the fixed tinyAM values. |",
  "| SSB and comparisons | Annual N, F, M, recruitment, total abundance and biomass outputs are reported. | Compare N/F at ages 2–8+, common-definition abundance, biomass and SSB where available. | Recruitment uses same-year age-2 N; source recruitment is age 1. Fbar and M are not treated as matched estimates. |",
  "",
  "The 0.54 sampling time is the midpoint of the reported 0.50–0.58 timing interval, used as a single-time approximation. The source model integrates over a timing interval."
)

fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent
  )

  retry_start <- function() {
    dat <- tinyAM::prepare_tam(
      data = obs, years = years, ages = ages,
      N_settings = N_settings, rec_settings = rec_settings,
      F_settings = F_settings, M_settings = M_settings,
      catch_settings = catch_settings, index_settings = index_settings
    )
    start <- tinyAM::make_par(dat)
    output_surface <- function(type, measure, numbers = FALSE) {
      rows <- source$outputs[
        source$outputs$type == type & source$outputs$measure == measure,
        , drop = FALSE
      ]
      age <- suppressWarnings(as.integer(as.character(rows$age)))
      grouped <- is.na(age) | !is.finite(age)
      age[grouped] <- suppressWarnings(as.integer(sub("\\+$", "", rows$age_group[grouped])))
      year <- suppressWarnings(as.integer(as.character(rows$year)))
      value <- suppressWarnings(as.numeric(as.character(rows$value)))
      if (numbers) value <- value * tinyAM:::.unit_scale(rows$unit, "numbers") / 1000
      keep <- year %in% years & age %in% ages
      rows <- data.frame(year = year[keep], age = age[keep], value = value[keep])
      if (anyDuplicated(paste(rows$year, rows$age))) {
        cli::cli_abort("Source starting values need one estimate per year and age.")
      }
      out <- matrix(NA_real_, length(years), length(ages),
        dimnames = list(as.character(years), as.character(ages)))
      out[cbind(match(rows$year, years), match(rows$age, ages))] <- rows$value
      out
    }
    N <- output_surface("population", "numbers_at_age", numbers = TRUE)
    F <- output_surface("mortality", "fishing_mortality_at_age")
    if (!is.finite(N[1L, 1L]) || N[1L, 1L] <= 0) {
      cli::cli_abort("Source N-at-age does not provide the first-year age-2 starting value.")
    }
    start$log_r0 <- log(N[1L, 1L])
    if (!is.null(start$log_n0)) {
      age <- as.integer(names(start$log_n0))
      value <- N[1L, match(age, ages)]
      keep <- is.finite(value) & value > 0
      start$log_n0[keep] <- log(value[keep])
    }
    if (length(start$log_r)) {
      year <- match(names(start$log_r), rownames(N))
      value <- N[year, 1L]
      keep <- !is.na(year) & is.finite(value) & value > 0
      start$log_r[keep] <- log(value[keep])
    }
    if (!is.null(start$log_f)) {
      year <- match(rownames(start$log_f), rownames(F))
      block_age <- dat$F_settings$age_block_start
      for (j in seq_along(block_age)) {
        value <- F[year, match(block_age[j], ages)]
        keep <- !is.na(year) & is.finite(value) & value > 0
        start$log_f[keep, j] <- log(value[keep])
      }
    }
    start
  }

  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      .assessment_cache(assessment_id, translation_tier, database = source),
      "initialization_attempts.rds"
    )
  )
}
