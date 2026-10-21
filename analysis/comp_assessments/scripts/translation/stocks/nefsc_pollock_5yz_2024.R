## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 1970:2023
ages <- 1:9
survey_names <- c(
  "NEFSC Albatross spring trawl survey",
  "NEFSC Bigelow spring trawl survey",
  "NEFSC Albatross fall trawl survey",
  "NEFSC Bigelow fall trawl survey"
)

inputs <- source$inputs
landings <- inputs[inputs$type == "catch" &
  inputs$measure == "landings_numbers_at_age", , drop = FALSE]
discards <- inputs[inputs$type == "catch" &
  inputs$measure == "discard_numbers_at_age", , drop = FALSE]
removals <- rbind(landings, discards)
combined_removals <- stats::aggregate(value ~ year + age, removals, sum)
catch_rows <- landings[match(
  paste(combined_removals$year, combined_removals$age),
  paste(landings$year, landings$age)
), , drop = FALSE]
catch_rows$value <- combined_removals$value
catch_rows$measure <- "numbers_at_age"
catch_rows$basis <- "numbers"
catch_rows$fleet <- "Combined removals"
catch_rows$partition <- NA_character_
catch_rows$sample_size <- NA_real_
catch_rows$unit <- "thousand fish"
catch_rows$source_type <- "reconstructed_source_input"
catch_rows$source_reference <- paste(unique(removals$source_reference), collapse = "; ")
catch_rows$transformation <- paste(
  "Analysis-only sum of official landings_numbers_at_age and",
  "discard_numbers_at_age by year and age."
)
catch_rows$notes <- paste(
  "Combined removals for the tinyAM fit; the original landings and discard",
  "streams remain separate in the canonical database."
)

age_index <- inputs[inputs$type == "index" &
  inputs$measure == "biomass_at_age" & inputs$survey %in% survey_names, , drop = FALSE]
index_compositions <- do.call(rbind, lapply(
  split(age_index, interaction(age_index$survey, age_index$year, drop = TRUE)),
  function(x) {
    x$value <- x$value / sum(x$value)
    x$measure <- "proportion_at_age"
    x$basis <- "proportion_biomass"
    x$unit <- "proportion"
    x$source_type <- "reconstructed_source_input"
    x$transformation <- paste(
      "Analysis-only normalization of native biomass_at_age values within",
      "each survey-year to create the composition input expected by the",
      "shared converter."
    )
    x$notes <- paste(
      "Normalized only for the tinyAM translation; original biomass-scale",
      "age values are retained in the canonical database."
    )
    x
  }
))
index_totals <- inputs[inputs$type == "index" &
  inputs$measure == "relative_biomass_index" & inputs$survey %in% survey_names, , drop = FALSE]
index_totals$measure <- "total_biomass"
index_totals$basis <- "biomass"
index_totals$source_type <- "reconstructed_source_input"
index_totals$transformation <- paste(
  "Analysis-only relabeling for the shared converter; values and relative",
  "kg-per-tow scale are unchanged."
)
index_totals$notes <- paste(
  "Relative survey index, not absolute biomass; retained on the source scale",
  "for reconstructing age-specific abundance indices."
)

translation_inputs <- inputs[!(inputs$type == "index" &
  inputs$measure %in% c("biomass_at_age", "relative_biomass_index")), , drop = FALSE]
translation_inputs <- rbind(translation_inputs, catch_rows,
  index_compositions, index_totals)

obs <- database_to_tam_obs(
  assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  index_weight_source = "source",
  surveys = survey_names,
  assumptions = source$assumptions
)

index_cv <- inputs[inputs$type == "index" &
  inputs$measure == "relative_standard_error" & inputs$survey %in% survey_names, , drop = FALSE]
index_key <- paste(obs$index$survey, obs$index$year)
cv_key <- paste(index_cv$survey, index_cv$year)
cv_index <- match(index_key, cv_key)
if (anyNA(cv_index)) stop("Every translated survey-age row needs its aggregate source CV.")
obs$index$relative_sd <- sqrt(log1p(index_cv$value[cv_index]^2))

## Background and comparisons ----
age_plus_group <- 9
comparison_outputs <- source$outputs
comparison_aggregates <- "abundance"
comparison_definitions <- list(
  ssb = list(
    status = "matched",
    definition = "Start-of-year spawning-stock biomass using the source maturity and SSB weight-at-age.",
    reason = "The official spawning fraction is zero and the same biological inputs are used in tinyAM."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 abundance in the same calendar year.",
    reason = "The source model represents recruitment at age 1."
  ),
  abundance = list(
    status = "matched",
    definition = "Total abundance over ages 1-9+ at January 1.",
    reason = "Both models report abundance over the same age range and annual timing."
  )
)

background <- c(
  "### Atlantic pollock 5YZ: official 2024 ASAP assessment", "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
  "| Years and ages | 1970-2023, ages 1-9+, one sex-combined stock. | Fit the same years and ages. | The complete source model grid is retained. |",
  "| N and recruitment | ASAP estimates annual abundance at age and age-1 recruitment. | Use exponential initial abundance and an IID abundance process. | tinyAM does not reproduce ASAP's initial-state integration or recruitment likelihood. |",
  "| F | One combined fishery has six age-specific selectivity blocks; directed and discard F are reported separately. | Combine F components in the reference and fit an age-specific AR1 F process around age-specific means. | The six selectivity blocks and separate removal components are not represented. |",
  "| M | Fixed M = 0.20 yr-1 at ages 1-9+. | Use the same supplied M values without an M process. | Fixed M agreement is an input check, not an estimated comparison. |",
  "| Catch | Total catch and dead discards have separate biomass likelihoods and age-composition likelihoods with different effective sample sizes. | Sum landings and discard numbers-at-age and fit one age-specific catch series with an estimated observation SD. | tinyAM does not reproduce the separate total-biomass likelihoods, composition weights, or their uncertainty. Zero age cells remain missing in the lognormal likelihood. |",
  "| Surveys and q | Four active spring/fall trawl indices use aggregate biomass indices, biomass-scale age compositions, source CVs and composition ESS. | Reconstruct age-specific abundance indices from the aggregate survey index, normalized age composition and source survey weights; estimate survey-by-age q and use the transformed aggregate CV as a shared log-scale SD. | The aggregate CV is applied to each age and age-composition sampling uncertainty is not propagated; the source annual q-deviation structure is simplified. |",
  "| Biology and SSB | Fixed maturity, three distinct annual weight matrices, and start-of-year SSB. | Use the same maturity and SSB weights; use the source survey weights to reconstruct survey abundance indices. | The source matrices retain their separate biological meanings. |",
  "| Outputs | Native point estimates include N-at-age, directed and discard F-at-age, biomass, SSB and Fbar ages 5-7. | Compare common year-age surfaces and derived total abundance. | Full time-age uncertainty is unavailable; derived totals have no uncertainty. |",
  "",
  "The 2024 model is the latest detailed, accepted official assessment found for this stock. The 2026 SASI file is a projection update only. tinyAM uses only the four active surveys; the summer index present in the ASAP input is inactive. The source observations and inputs are not adjusted using fitted q or other estimated parameters."
)

## Model ----
source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

fit_model <- function(...) tinyAM::fit_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = list(process = "iid", init = "exp"),
  F_settings = list(process = "ar1", mu_form = ~ factor(age), mean_ages = 5:7),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~ 0 + q_key, sd_form = ~0,
    sd_supplied = ~relative_sd),
  silent = silent,
  ...
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs,
      "population", "numbers_at_age", years, ages)
    f <- .translation_start_surface(source$outputs,
      "mortality", "fishing_mortality_at_age", years, ages)
    list(
      log_r0 = log(n[1L, 1L]),
      log_r = log(n[-1L, 1L]),
      log_n0 = log(n[1L, -1L]),
      log_f = log(f)
    )
  }, attempt_file = file.path(
    "analysis/comp_assessments/results/cache", assessment_id,
    "initialization_attempts.rds"
  ))
}
