## Observations ----
years <- 1959:2025
ages <- 1:10
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("FO-GFS-Q3", "FO-GFS-Q1")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("FO-GFS-Q3" = 0.65, "FO-GFS-Q1" = 0.225),
  assumptions = source$assumptions
)

obs$catch$sd_block <- factor(ifelse(
  obs$catch$age == 10, "age10_plus", "ages2_9"
), levels = c("ages2_9", "age10_plus"))
obs$index$q_age <- factor(ifelse(
  obs$index$age <= 7, paste0("age", obs$index$age), "ages8_9"
), levels = c(paste0("age", 1:7), "ages8_9"))
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_age <- factor(ifelse(
  obs$index$age == 1, "age1", "ages2_9"
), levels = c("age1", "ages2_9"))
obs$index$sd_block <- interaction(
  obs$index$survey, obs$index$sd_age,
  drop = TRUE, lex.order = TRUE
)

## Comparison ----
age_plus_group <- 10L
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Sum of N-at-age over ages 1–10+ in the same year, in thousand fish.",
    reason = "The source total is derived from the complete native SAM N-at-age surface; uncertainty is unavailable for this derived sum."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass over ages 1–10+ using the same annual stock weights.",
    reason = "The source and tinyAM use the same age range and biological inputs."
  ),
  ssb = list(
    status = "matched",
    definition = "Mature biomass over ages 1–10+ using source stock weights and maturity; pre-spawning F and M fractions are zero.",
    reason = "The native SAM SSB definition is reproduced by the same source biology and spawning-time convention."
  ),
  recruitment = list(
    status = "matched",
    definition = "Beginning-of-year age-1 abundance in the same calendar year, in thousand fish.",
    reason = "The age and units match; the recruitment and initial-state parameterizations differ."
  ),
  F = list(
    status = "approximate",
    definition = "Annual F-at-age on ages 1–10+, with source Fbar defined as the mean over ages 3–7.",
    reason = "The SAM model has no estimated age-1 F state and shares the age-7–9 state; tinyAM estimates F over its full age grid."
  ),
  M = list(
    status = "matched",
    definition = "Fixed natural mortality of 0.2 per year at ages 1–10+.",
    reason = "The same fixed source M surface is supplied to tinyAM; neither model estimates an M process."
  )
)

background <- c(
  "### Faroe Plateau cod: 2025 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1959–2025 at ages 1–10+, with age 10 as the terminal plus group. | Fit the same 1959–2025 and 1–10+ grid. | The full fitted period and age range are retained. |",
  "| N and recruitment | Recruitment enters at age 1 with its own variance; older ages 2–10 share a survival-process variance. | Use a random-walk recruitment process and IID cohort-survival deviations with one shared SD for older ages. | The source's recruitment and older-survival processes are represented separately; initial-state treatment remains simplified. |",
  "| F | SAM estimates F at ages 2–10, shares the age-7–9 state, separates age-2 process variance from ages 3–10, and reports Fbar for ages 3–7. | Use a correlated random walk across the tinyAM age grid and define mean F over ages 3–7. | tinyAM estimates age-1 F, uses one F-process SD, and does not reproduce SAM's shared age-state or full covariance structure. |",
  "| M and biology | M is fixed at 0.2 per year; annual stock, catch, landings, and discard weights, maturity, and spawning fractions are supplied. | Use the same fixed M and biological surfaces. | Native biological inputs and zero pre-spawning mortality fractions are retained. |",
  "| Catch | Total catch-at-age in thousand fish; ages 2–9 share an observation SD and age 10+ has a separate SD. | Keep the observed catch-at-age and both source SD groups. | The remaining source likelihood structure is simplified. |",
  "| Surveys and q | FO-GFS-Q3 and FO-GFS-Q1 observe ages 1–9 at times 0.65 and 0.225. Each has separate q at ages 1–7 and a shared q for ages 8–9; summer residuals are AR(1) across ages. | Retain both native index series, timing, q groups, and age-specific observation-SD groups. | tinyAM treats survey residuals as independent and omits the summer age correlation. |",
  "| Comparisons | Native SAM reports N, F, age-1 recruitment, SSB, total biomass, and Fbar. | Compare N/F at age and the common annual metrics; derive total abundance by summing N at age. | Derived total abundance has no uncertainty; age-1 F remains an approximate comparison. |",
  "",
  "The fit first uses tinyAM's normal initialization. If that fails numerical checks, one retry initializes N, F, and q from official estimates; these values are starts only, and tinyAM estimates them freely."
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
      F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:7),
      M_settings = list(
        process = "off", mu_form = NULL,
        mu_supplied = ~M_assumption
      ),
      catch_settings = list(sd_form = ~0 + sd_block),
      index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block),
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
        F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:7),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(sd_form = ~0 + sd_block),
        index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
      )
      starts <- tinyAM::make_par(model_dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      starts$log_r0 <- log(source_N[1L, 1L])
      starts$log_r <- stats::setNames(log(source_N[-1L, 1L]), as.character(years[-1L]))
      starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])

      source_F <- source$outputs[
        source$outputs$type == "mortality" &
          source$outputs$measure == "fishing_mortality_at_age",
        , drop = FALSE
      ]
      f_matrix <- matrix(
        NA_real_, length(years), length(ages),
        dimnames = list(as.character(years), as.character(ages))
      )
      f_matrix[cbind(match(source_F$year, years), match(source_F$age, ages))] <- source_F$value
      positive_F <- is.finite(f_matrix) & f_matrix > 0
      starts$log_f[positive_F] <- log(f_matrix[positive_F])

      index <- obs$index
      q_design <- stats::model.matrix(~0 + q_key, data = index)
      index_rows <- cbind(match(index$year, years), match(index$age, ages))
      source_n <- source_N[index_rows]
      source_f <- f_matrix[index_rows]
      n_at_survey <- source_n * exp(-(source_f + 0.2) * index$samp_time)
      q_start <- vapply(seq_len(ncol(q_design)), function(j) {
        selected <- q_design[, j] == 1 & is.finite(index$obs) &
          index$obs > 0 & is.finite(n_at_survey) & n_at_survey > 0
        if (!any(selected)) return(NA_real_)
        stats::median(log(index$obs[selected] / n_at_survey[selected]))
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
