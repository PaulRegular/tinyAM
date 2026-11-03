## Observations ----
years <- 1957:2019
ages <- 1:9
assessment_id <- source$assessment$assessment_id[[1L]]
inputs <- source$inputs

catch_rows <- inputs$type == "catch" & inputs$measure == "numbers_at_age"
source_catch <- inputs[catch_rows, , drop = FALSE]
catch_groups <- split(
  seq_len(nrow(source_catch)),
  paste(source_catch$year, source_catch$age, sep = "\034")
)
combined_catch <- do.call(rbind, lapply(catch_groups, function(i) {
  rows <- source_catch[i, , drop = FALSE]
  row <- rows[1L, , drop = FALSE]
  values <- as.numeric(rows$value)
  row$value <- if (all(is.na(values))) NA_real_ else sum(values, na.rm = TRUE)
  row$fleet <- "Combined catch"
  row$region <- NA_character_
  row$source_type <- "analysis_derived"
  row$transformation <- "Sum source catch-at-age across Catch 6aN and Catch 6aS/7bc."
  row$notes <- paste(
    "Translation-only aggregation of the two official catch fleets.",
    "The canonical database retains each source fleet separately."
  )
  row
}))
inputs <- rbind(inputs[!catch_rows, , drop = FALSE], combined_catch)
unresolved_indices <- inputs$type == "index" & inputs$measure == "numbers_at_age"
inputs$unit[unresolved_indices] <- "native survey index"

surveys <- c("MS_HERAS", "IBTS_Q1", "IBTS_Q4")
obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0.67

q_age <- ifelse(obs$index$age >= 8, "8plus", paste0("age", obs$index$age))
obs$index$q_key <- interaction(
  obs$index$survey, q_age, drop = TRUE, lex.order = TRUE
)

sd_age <- ifelse(
  obs$index$survey == "MS_HERAS",
  ifelse(obs$index$age == 1, "age1",
    ifelse(obs$index$age == 2, "age2",
      ifelse(obs$index$age <= 6, "ages3_6", "ages7_9")
    )
  ),
  ifelse(obs$index$survey == "IBTS_Q1",
    ifelse(obs$index$age == 2, "age2",
      ifelse(obs$index$age <= 5, "ages3_5",
        ifelse(obs$index$age <= 7, "ages6_7", "ages8_9")
      )
    ),
    ifelse(obs$index$age == 2, "age2",
      ifelse(obs$index$age <= 6, "ages3_6", "ages7_9")
    )
  )
)
obs$index$sd_block <- interaction(
  obs$index$survey, sd_age, drop = TRUE, lex.order = TRUE
)

## Comparisons ----
age_plus_group <- 9L
source_F <- source$outputs[
  source$outputs$type == "mortality" &
    source$outputs$measure == "fishing_mortality_at_age",
  , drop = FALSE
]
F_groups <- split(
  seq_len(nrow(source_F)),
  paste(source_F$year, source_F$age, sep = "\034")
)
combined_F <- do.call(rbind, lapply(F_groups, function(i) {
  rows <- source_F[i, , drop = FALSE]
  row <- rows[1L, , drop = FALSE]
  row$value <- sum(as.numeric(rows$value))
  row$region <- NA_character_
  row$fleet <- NA_character_
  row$source_type <- "derived_source_output"
  row$se <- row$lwr <- row$upr <- NA_real_
  row$notes <- paste(
    "Derived total F by summing the two source area-specific F-at-age series.",
    "The original area-specific values remain in the source outputs."
  )
  row
}))
comparison_outputs <- rbind(
  source$outputs[!(
    source$outputs$type == "mortality" &
      source$outputs$measure == "fishing_mortality_at_age"
  ), , drop = FALSE],
  combined_F
)
comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Beginning-of-year abundance summed over ages 1-9+.",
    reason = "Derived from the full reported N-at-age series; no uncertainty is available for the sum."
  ),
  recruitment = list(
    status = "matched",
    definition = "Beginning-of-year recruitment at age 1 in thousand fish.",
    reason = "The age and source scale match; tinyAM uses a separate IID recruitment process."
  ),
  N = list(
    status = "matched",
    definition = "Beginning-of-year N at ages 1-9+ in thousand fish.",
    reason = "The model age grid and reported population ages match."
  ),
  F = list(
    status = "approximate",
    definition = "Total F at age formed by summing the 6.a North and 6.a South/7.b-c components.",
    reason = "tinyAM fits one pooled F process; the official model estimates separate area-specific F states."
  ),
  F_bar = list(
    status = "approximate",
    definition = "N-weighted total F over ages 3-6 in both models.",
    reason = "The comparison uses the same rate and abundance weighting, while tinyAM pools the source area-specific F processes."
  ),
  M = list(
    status = "matched",
    definition = "Fixed source M-at-age input, repeated over 1957-2019.",
    reason = "This is an input-consistency check; M is not independently estimated by either model."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass over ages 1-9+ using source N and stock weight-at-age.",
    reason = "The comparison uses the same source weights and age range."
  ),
  ssb = list(
    status = "matched",
    definition = "Spawning biomass using source N, weight, maturity, and 0.67-year survival to spawning.",
    reason = "The source spawning fractions and biological inputs are applied to both definitions."
  )
)

background <- c(
  "### Herring in ICES divisions 6.a and 7.b-c: Official 2020 FLSAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1957-2019; ages 1-9+, with age 9 as the plus group. | Retain the same years and age grid. | None in the year-age coverage. |",
  "| N and recruitment | No stock-recruit relationship; recruitment and ages 2-9 survival have separate process variances. | Use an IID recruitment process around a constant mean and IID older-cohort survival deviations. | Initial-state integration and the source's detailed cohort structure are simplified. |",
  "| F | SAM estimates separate area-specific F states, with age sharing and correlated temporal changes. Fbar is reported for ages 3-6. | Sum the two source catch streams and fit one correlated F random walk; calculate Fbar over ages 3-6. | The two fleets and their distinct selectivity/process structures are pooled. |",
  "| M and biology | Fixed age-specific M; annual stock and catch weights, maturity, and pre-spawning F/M fractions of 0.67. | Retain fixed M and source stock biology; use a 0.67 spawning time. | Catch fleets are pooled for tinyAM; the source database keeps them separate. |",
  "| Catch | Two lognormal catch fleets with age-specific SD groups. | Sum catch numbers at age across fleets and estimate one catch SD. | Separate fleet selectivity and catch residual structures are not reproduced. |",
  "| Surveys | MS_HERAS, IBTS_Q1, and IBTS_Q4 with age-specific q (ages 8-9 shared), grouped SDs, and AR(1) residual correlation across ages. | Retain each survey, q-sharing groups, sampling-time interval midpoints, and reported SD groups; leave values on the native scale because source units are unresolved. | tinyAM uses independent survey-age residuals; the source's cross-age residual correlations are omitted. |",
  "| Comparisons | Reported N and two area-specific F-at-age series, recruitment, biomass, SSB, and Fbar. | Compare N and total F at age; derive common total abundance, biomass, SSB, and N-weighted Fbar. | Summed abundance has no uncertainty; the total-F comparison is approximate because tinyAM pools area-specific F. |",
  "",
  "Catch aggregation and total F are analysis-only transformations. No source observations are adjusted using fitted values. The catch SD pooling and independent survey residuals are explicit likelihood approximations."
)

settings <- list(
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:6),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()

  fit_model <- function(start_par = NULL) {
    do.call(tinyAM::fit_tam, c(
      list(
        data = obs, years = years, ages = ages,
        silent = silent, start_par = start_par
      ),
      settings
    ))
  }

  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      model_dat <- do.call(tinyAM::prepare_tam, c(
        list(data = obs, years = years, ages = ages),
        settings
      ))
      starts <- tinyAM::make_par(model_dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      source_F <- .translation_start_surface(
        combined_F, "mortality", "fishing_mortality_at_age", years, ages
      )
      starts$log_r0 <- log(source_N[1L, 1L])
      starts$log_r <- stats::setNames(log(source_N[-1L, 1L]), as.character(years[-1L]))
      starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])
      positive_F <- is.finite(source_F) & source_F > 0
      starts$log_f[positive_F] <- log(source_F[positive_F])
      starts
    },
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      assessment_id, "initialization_attempts.rds"
    )
  )
}
