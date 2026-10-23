## Observations ----
years <- 1957:2025
ages <- 1:10
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("FO-GFS-Q3", "FO-GFS-Q1")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

catch_weight <- source$inputs[
  source$inputs$assessment_id == assessment_id &
    source$inputs$type == "catch" &
    source$inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
catch_key <- function(x) paste(x$year, x$age, sep = "\r")
catch_index <- match(catch_key(obs$catch), catch_key(catch_weight))
if (anyNA(catch_index)) stop("Catch precision weights do not cover the catch grid.")
obs$catch$relative_sd <- 1 / sqrt(catch_weight$value[catch_index])
obs$catch$sd_scale <- 1

obs$index$q_age <- ifelse(
  obs$index$survey == "FO-GFS-Q3" & obs$index$age <= 6,
  paste0("age", obs$index$age),
  ifelse(obs$index$survey == "FO-GFS-Q3", "ages7_8",
    ifelse(obs$index$age <= 5, paste0("age", obs$index$age), "ages6_7"))
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_age <- ifelse(
  obs$index$survey == "FO-GFS-Q3" & obs$index$age <= 2,
  "ages1_2",
  ifelse(obs$index$survey == "FO-GFS-Q3", "ages3_8",
    ifelse(obs$index$age == 1, "age1", "ages2_7"))
)
obs$index$sd_group <- interaction(
  obs$index$survey, obs$index$sd_age,
  drop = TRUE, lex.order = TRUE
)

## Comparison ----
age_plus_group <- 10L
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Sum of N-at-age over ages 1–10+ in the same year, in thousand fish.",
    reason = "Derived from complete native SAM N-at-age; uncertainty is unavailable for the sum."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass over ages 1–10+ using the same annual stock weights.",
    reason = "Both calculations use the native stock-weight surface and the same age range."
  ),
  ssb = list(
    status = "matched",
    definition = "Mature biomass over ages 1–10+ using native stock weights and maturity.",
    reason = "Both calculations use the source biological surfaces and zero pre-spawning mortality fractions."
  ),
  recruitment = list(
    status = "matched",
    definition = "Beginning-of-year age-1 abundance in the same year, in thousand fish.",
    reason = "The age and unit match; the process and initial-state parameterizations differ."
  ),
  F_bar = list(
    status = "approximate",
    definition = "F averaged over ages 3–7; tinyAM uses its population-weighted mean over those ages.",
    reason = "The reported SAM Fbar is an arithmetic mean, while tinyAM reports an N-weighted mean."
  ),
  F = list(
    status = "approximate",
    definition = "Annual F-at-age over ages 1–10+.",
    reason = "SAM shares its oldest F states and uses correlated age increments; tinyAM estimates its own full-age surface."
  ),
  M = list(
    status = "matched",
    definition = "Fixed natural mortality of 0.2 per year at ages 1–10+.",
    reason = "The same fixed mortality surface is supplied; neither model estimates M."
  )
)

background <- c(
  "### Faroe haddock: 2025 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1957–2025 at ages 1–10+, with age 10 as the catch plus group. | Fit the same years and ages. | The source grid is retained. |",
  "| N and recruitment | Age-1 recruitment follows a random walk; N variance differs between age 1 and ages 2–10. | Use IID N survival deviations and a separate random-walk recruitment process. | Variance sharing and initial-state integration differ. |",
  "| F | F states are fitted at ages 1–8, with ages 9–10 shared; age innovations are AR(1). Fbar is the arithmetic mean over ages 3–7. | Use a correlated random walk over ages 1–10 and mean ages 3–7. | tinyAM estimates the full age surface and reports an N-weighted Fbar. |",
  "| M and biology | M is fixed at 0.2 per year. Annual stock, catch, landings, and discard weights, maturity, and spawning fractions are supplied. | Use the same M and native stock-weight and maturity surfaces. | Zero pre-spawning F and M fractions are retained. |",
  "| Catch | Catch numbers-at-age are in thousand fish. One observation SD is shared across ages; 2025 has relative precision 0.5. | Retain catch-at-age and one estimated SD multiplier, with relative SD offsets of 1 for 1957–2024 and sqrt(2) for 2025. | This represents the relative precision by scaling the observation SD; the weighted SAM likelihood is not reproduced exactly. |",
  "| Surveys and q | FO-GFS-Q3 (1996–2025, ages 1–8) and FO-GFS-Q1 (1994–2025, ages 1–7), at times 0.65 and 0.225. Q3 ages 1–6 and Q1 ages 1–5 are separate; older sampled ages are pooled. | Retain both native series, timing, and q sharing; estimate the reported survey-age SD groups. | tinyAM does not represent the source AR(1) age correlation in survey residuals. |",
  "| Comparisons | SAM reports N, F, recruitment, SSB, total biomass, Fbar, q, and fitted observations. | Compare age-specific N/F and annual abundance, biomass, SSB, recruitment, and Fbar. | Total abundance uncertainty is unavailable; Fbar definitions differ. |",
  "",
  "The fit first uses tinyAM's normal initialization. If numerical checks fail, one retry uses aligned SAM N, F, and q estimates as starting values only."
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
        F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 3:7),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(
          sd_form = ~0 + sd_scale,
          sd_supplied = ~relative_sd
        ),
        index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_group)
      )
      starts <- tinyAM::make_par(model_dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      starts$log_r0 <- log(source_N[1L, 1L])
      starts$log_r <- stats::setNames(log(source_N[-1L, 1L]), as.character(years[-1L]))
      starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])

      source_F <- .translation_start_surface(
        source$outputs, "mortality", "fishing_mortality_at_age", years, ages
      )
      positive_F <- is.finite(source_F) & source_F > 0
      starts$log_f[positive_F] <- log(source_F[positive_F])
      index_rows <- cbind(match(obs$index$year, years), match(obs$index$age, ages))
      n_survey <- source_N[index_rows]
      f_survey <- source_F[index_rows]
      n_survey <- n_survey * exp(-(f_survey + 0.2) * obs$index$samp_time)
      q_design <- stats::model.matrix(~0 + q_key, data = obs$index)
      q_start <- vapply(seq_len(ncol(q_design)), function(j) {
        keep <- q_design[, j] == 1 & is.finite(obs$index$obs) &
          obs$index$obs > 0 & is.finite(n_survey) & n_survey > 0
        if (!any(keep)) return(NA_real_)
        stats::median(log(obs$index$obs[keep] / n_survey[keep]))
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
