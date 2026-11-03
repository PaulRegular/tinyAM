## Observations ----
years <- 1993:2025
ages <- 0:8
assessment_id <- source$assessment$assessment_id[[1]]
survey <- "FRA_IRL_IBTS"

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = survey,
  assumptions = source$assumptions
)

obs$catch$obs[obs$catch$obs == 0] <- NA_real_

catch_weight <- source$inputs[
  source$inputs$assessment_id == assessment_id &
    source$inputs$type == "catch" &
    source$inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
catch_key <- function(x) paste(x$year, x$age, sep = "\r")
catch_index <- match(catch_key(obs$catch), catch_key(catch_weight))
if (anyNA(catch_index)) stop("Catch weights do not cover the translated catch grid.")
obs$catch$relative_sd <- 1 / sqrt(catch_weight$value[catch_index])
obs$catch$sd_scale <- 1

obs$index$sd_age <- ifelse(obs$index$age <= 5, "ages0_5", "ages6_7")
obs$index$sd_group <- interaction(
  obs$index$survey, obs$index$sd_age,
  drop = TRUE, lex.order = TRUE
)

## Comparison ----
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "biomass", "ssb", "recruitment", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total start-of-year abundance over ages 0–8+ in thousand fish.",
    reason = "Official total abundance is derived from the same complete SAM N-at-age surface."
  ),
  biomass = list(
    status = "matched",
    definition = "Total stock biomass over ages 0–8+ using the native stock weights.",
    reason = "Both values use the source stock-weight surface and the same model ages."
  ),
  ssb = list(
    status = "matched",
    definition = "Spawning-stock biomass over ages 0–8+ using native stock weights, maturity and spawning fractions.",
    reason = "The source pre-spawning mortality fractions are zero, and the same biological surfaces are supplied."
  ),
  recruitment = list(
    status = "matched",
    definition = "Start-of-year abundance at age 0 in thousand fish.",
    reason = "Both values use the fitted SAM age-0 abundance."
  ),
  F_bar = list(
    status = "approximate",
    definition = "F averaged over ages 3–5.",
    reason = "The official Fbar is the arithmetic mean; tinyAM reports its population-weighted mean over the same ages."
  ),
  F = list(
    status = "approximate",
    definition = "Annual fishing mortality at ages 0–8+.",
    reason = "SAM shares the age-7/8 state and uses correlated age innovations; tinyAM estimates its own full-age process."
  ),
  M = list(
    status = "matched",
    definition = "Fixed natural mortality at ages 0–8+.",
    reason = "The same age-specific M surface is supplied to both assessments."
  )
)

background <- c(
  "### Celtic Sea haddock: 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1993–2025 at ages 0–8+, with age 8 as the plus group. | Fit the same years and ages. | The source grid is retained. |",
  "| N and recruitment | Age-0 recruitment follows a log random walk; N variance separates age 0 from ages 1–8. | Use IID N survival deviations and a separate random-walk recruitment process. | Variance sharing and initial-state integration differ. |",
  "| F | F follows a temporal random walk with AR(1) age correlation; ages 7 and 8 share a state, and age-0 F is fixed at zero. Fbar is the arithmetic mean over ages 3–5. | Use a correlated random walk over ages 0–8 and mean ages 3–5. | tinyAM estimates age-0 F and reports an N-weighted Fbar. |",
  "| M and biology | M is fixed by age. Annual stock, catch, landings and discard weights and fixed maturity are supplied. | Use the same M and native stock-weight and maturity surfaces. | Zero pre-spawning F and M fractions are retained. |",
  "| Catch | Combined landings and discards are in thousand fish; source relative precision weights are retained. | Retain catch-at-age and scale the shared observation SD by the inverse square root of the source weight. | This represents relative precision; the weighted SAM likelihood is not reproduced exactly. |",
  "| Survey and q | One VAST-standardized French–Irish Q4 index covers ages 0–7 at sampling time 0.875; q is age-specific. | Retain the native index scale, timing and age-specific q; estimate two observation SD groups for ages 0–5 and 6–7. | The index physical unit is unresolved, and tinyAM does not reproduce the full SAM likelihood. |",
  "| Outputs | SAM reports N, F, recruitment, SSB, total biomass and Fbar. | Compare N/F-at-age and the available aggregate quantities. | Total abundance is derived from SAM N-at-age and has no uncertainty. |",
  "",
  "The source fit first uses tinyAM's default initialization. If its numerical checks fail, one retry uses year- and age-aligned SAM N, F and q estimates as starting values only."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) {
    tinyAM::fit_tam(
      data = obs,
      years = years,
      ages = ages,
      N_settings = list(process = "iid", init = "exp", sd_form = ~1),
      rec_settings = list(process = "rw", mu_form = NULL),
      F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:5),
      M_settings = list(
        process = "off", mu_form = NULL,
        mu_supplied = ~M_assumption
      ),
      catch_settings = list(
        sd_form = ~0 + sd_scale,
        sd_supplied = ~relative_sd
      ),
      index_settings = list(
        q_form = ~0 + q_key,
        sd_form = ~0 + sd_group
      ),
      silent = silent,
      start_par = start_par
    )
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      model_dat <- tinyAM::prepare_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = list(process = "iid", init = "exp", sd_form = ~1),
        rec_settings = list(process = "rw", mu_form = NULL),
        F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:5),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(
          sd_form = ~0 + sd_scale,
          sd_supplied = ~relative_sd
        ),
        index_settings = list(
          q_form = ~0 + q_key,
          sd_form = ~0 + sd_group
        )
      )
      starts <- tinyAM::make_par(model_dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      starts$log_r0 <- log(source_N[1L, 1L])
      starts$log_r <- stats::setNames(log(source_N[-1L, 1L]), as.character(years[-1L]))
      starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])

      source_F <- matrix(NA_real_, length(years), length(ages),
        dimnames = list(as.character(years), as.character(ages))
      )
      source_F[, -1L] <- .translation_start_surface(
        source$outputs, "mortality", "fishing_mortality_at_age", years, ages[-1L]
      )
      positive_F <- is.finite(source_F) & source_F > 0
      starts$log_f[positive_F] <- log(source_F[positive_F])

      index_rows <- cbind(match(obs$index$year, years), match(obs$index$age, ages))
      index_numbers <- source_N[index_rows]
      index_F <- source_F[index_rows]
      pop_rows <- match(
        paste(obs$index$year, obs$index$age),
        paste(obs$pop_covariates$year, obs$pop_covariates$age)
      )
      index_M <- obs$pop_covariates$M_assumption[pop_rows]
      vulnerable_numbers <- index_numbers * exp(-(index_F + index_M) * obs$index$samp_time)
      q_design <- stats::model.matrix(~0 + q_key, data = obs$index)
      q_start <- vapply(seq_len(ncol(q_design)), function(j) {
        keep <- q_design[, j] == 1 & is.finite(obs$index$obs) &
          obs$index$obs > 0 & is.finite(vulnerable_numbers) & vulnerable_numbers > 0
        if (!any(keep)) return(NA_real_)
        stats::median(log(obs$index$obs[keep] / vulnerable_numbers[keep]))
      }, numeric(1))
      names(q_start) <- colnames(q_design)
      q_start <- q_start[is.finite(q_start)]
      starts$log_q[names(q_start)] <- q_start
      starts
    },
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      assessment_id, "initialization_attempts.rds"
    )
  )
}
