## Observations ----
years <- 1999:2025
ages <- 0:7
survey <- "Combined FR-EVHOE and IE-IGFS VAST recruitment survey"
inputs <- source$inputs
density <- inputs$type == "index" & inputs$survey %in% survey
inputs$transformation[density] <- "No/km2 relabelled fish per km2 without changing numerical values."
inputs$unit[density] <- "fish per km2"
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, surveys = survey,
  sampling_times = stats::setNames(0.875, survey), assumptions = source$assumptions)
obs$index$q_group <- factor(ifelse(obs$index$age == 0, "age0", "ages1_2"))
obs$index$sd_group <- obs$index$q_group

N_settings <- list(process = "iid", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "cor_rw", mu_form = NULL,
  age_breaks = c(0:6, 8), mean_ages = 2:5)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_group, sd_form = ~0 + sd_group)
age_plus_group <- 7L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "Source N with shared January weights and production maturity over ages 0–7+.",
  reason = "The production spawning fractions were not recovered; this is common-definition beginning-year mature biomass."))

## Background ----
background <- c(
  "### Celtic Sea whiting: preliminary comparison with the Official 2026 assessment", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1999–2025; ages 0–7+. | The same historical grid; no 2026 forecasts are fitted. | No source years or ages are removed. |",
  "| N | SAM cohort dynamics; exact current configuration unavailable. | IID survival deviations, random-walk recruitment and free initial abundance. | Initial-state treatment and process variance sharing are approximate. |",
  "| F | The framework shares F for ages 6+ and correlates age increments. | A correlated random walk with one state for ages 6–7+. | The structure follows framework evidence, not a recovered 2026 parameter file. |",
  "| M | Fixed production values 1.22, 0.86, 0.65, 0.50, 0.43, 0.40, 0.38, 0.36. | Supply the same age vector. | No earlier annex vector is substituted. |",
  "| Catch | Reported total numbers; early older-age discards were estimated by SAM. | Retain those reported counts with one estimated lognormal SD. | Composition and aggregate catch likelihoods and unknown landings fractions are not reproduced. |",
  "| Index | VAST Q4 density ages 0–2; a separate commercial biomass/LPUE index. | Fit VAST at an explicit Q4 midpoint 0.875; share q and SD across ages 1–2. | The commercial biomass index has no matching age composition and is omitted. No/km² is converted to thousand fish/km², preserving its area denominator. |",
  "| Biology | Annual January stock weights and fixed production maturity. | Retain these inputs; report common-definition beginning-year mature biomass. | Spawning timing is unresolved, so native SSB equivalence is not claimed. |",
  "", "No missing observations or likelihood weights are manufactured. Package-default starts are tried first; a numerical failure permits one source N/F retry."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent)
  fit <- fit_default_first(fit_model, retry_start = function() {
    dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages,
      N_settings = N_settings, rec_settings = rec_settings,
      F_settings = F_settings, M_settings = M_settings,
      catch_settings = catch_settings, index_settings = index_settings)
    # Source F is shared at ages 6–7; use the block-start ages, not a pooled rate.
    unblocked <- dat
    unblocked$F_settings$age_breaks <- NULL
    unblocked$F_settings$age_blocks <- factor(ages, levels = ages)
    unblocked$F_settings$age_block_start <- ages
    starts <- source_start_parameters(source, unblocked)
    starts$log_f <- starts$log_f[, as.character(dat$F_settings$age_block_start), drop = FALSE]
    colnames(starts$log_f) <- levels(dat$F_settings$age_blocks)
    starts
  }, attempt_file = file.path(.assessment_cache(source$assessment$assessment_id,
    translation_tier), "initialization_attempts.rds"))
}
