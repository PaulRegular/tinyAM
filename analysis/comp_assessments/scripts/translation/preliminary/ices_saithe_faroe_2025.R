## Observations ----

years <- 1961:2025
ages <- 3:15
assessment_id <- source$assessment$assessment_id[[1L]]
surveys <- c("SummerSurvey", "SpringSurvey")

translation_inputs <- source$inputs
add_terminal_weight <- function(input_type, measure, output_measure) {
  existing <- translation_inputs[
    translation_inputs$type == input_type &
      translation_inputs$measure == measure &
      translation_inputs$year == "2025",
    , drop = FALSE
  ]
  missing_ages <- setdiff(ages, as.integer(existing$age))
  if (!length(missing_ages)) return(invisible(NULL))
  estimates <- source$outputs[
    source$outputs$type == "biology" &
      source$outputs$measure == output_measure &
      source$outputs$year == 2025 &
      source$outputs$age %in% missing_ages,
    , drop = FALSE
  ]
  if (!setequal(as.integer(estimates$age), missing_ages)) {
    stop("Native SAM terminal weight estimates do not cover the missing ages.")
  }
  added <- translation_inputs[rep(NA_integer_, length(missing_ages)), , drop = FALSE]
  added$assessment_id <- assessment_id
  added$type <- input_type
  added$measure <- measure
  added$basis <- "kg_per_fish"
  added$fleet <- added$survey <- added$sex <- added$region <- added$season <- ""
  added$year <- 2025L
  added$year_basis <- "calendar_year"
  added$age <- as.integer(estimates$age)
  added$value <- as.numeric(estimates$value)
  added$unit <- estimates$unit
  added$sampling_time <- NA_real_
  added$source_type <- "reconstructed_source_input"
  added$source_reference <- estimates$source_reference
  added$transformation <- "Used the official SAM-estimated 2025 surface because the native stock and catch weight input files have no 2025 values."
  added$notes <- "Translation-only value derived from the matching official model output; the canonical input table preserves the missing source input."
  translation_inputs <<- rbind(translation_inputs, added)
  invisible(NULL)
}
add_terminal_weight("weight", "weight_at_age", "weight_at_age")
add_terminal_weight("catch_weight", "weight_at_age", "catch_weight_at_age")

obs <- database_to_tam_obs(
  assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

catch_precision <- source$inputs[
  source$inputs$assessment_id == assessment_id &
    source$inputs$type == "catch" &
    source$inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
catch_key <- function(x) paste(x$year, x$age, sep = "\r")
catch_index <- match(catch_key(obs$catch), catch_key(catch_precision))
if (anyNA(catch_index)) stop("Catch precision weights do not cover the source catch grid.")
obs$catch$relative_sd <- 1 / sqrt(as.numeric(catch_precision$value[catch_index]))
obs$catch$sd_scale <- 1

obs$index$q_age <- factor(
  ifelse(obs$index$age <= 8, paste0("age", obs$index$age), "ages9_10"),
  levels = c(paste0("age", 3:8), "ages9_10")
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_age <- factor(obs$index$age, levels = ages)
obs$index$sd_block <- interaction(
  obs$index$survey, obs$index$sd_age,
  drop = TRUE, lex.order = TRUE
)

## Comparison ----

age_plus_group <- 15L
comparison_outputs <- source$outputs
index_predictions <- comparison_outputs$type == "index" &
  comparison_outputs$measure == "predicted_index"
comparison_outputs$unit[index_predictions] <- "thousand fish per effort-hour"
q_estimates <- comparison_outputs$type == "catchability" &
  comparison_outputs$measure == "q"
scaled_columns <- c("value", "se", "lwr", "upr")
comparison_outputs[index_predictions | q_estimates, scaled_columns] <-
  comparison_outputs[index_predictions | q_estimates, scaled_columns] * 0.001
comparison_outputs$unit[q_estimates] <- "per effort-hour"
comparison_outputs$notes[index_predictions | q_estimates] <- paste(
  comparison_outputs$notes[index_predictions | q_estimates],
  "Index and q values scaled by 0.001 for tinyAM's thousand-fish index; the database retains native values."
)
comparison_aggregates <- c("abundance", "biomass", "ssb", "recruitment", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Sum of official N-at-age over ages 3–15+ in the same calendar year.",
    reason = "Age range and units match; the source total is derived from the fitted age-specific population surface and has no separate uncertainty."
  ),
  biomass = list(
    status = "matched",
    definition = "Common-definition total biomass over ages 3–15+ using shared translated stock weights.",
    reason = "Official N and tinyAM N are evaluated with the original annual input weights and SAM-estimated 2025 weights. Native SAM biomass uses fitted weights and is retained separately."
  ),
  ssb = list(
    status = "matched",
    definition = "Common-definition mature biomass over ages 3–15+ using shared translated stock weights and maturity.",
    reason = "Both N surfaces use the same translated biology and zero pre-spawning mortality fractions. Native SAM SSB uses fitted stock weights and is retained separately."
  ),
  recruitment = list(
    status = "matched",
    definition = "Beginning-of-year abundance at age 3 in the same calendar year, in thousand fish.",
    reason = "Recruitment age and units match; the recruitment and initial-state parameterizations differ."
  ),
  F = list(
    status = "approximate",
    definition = "Annual fishing mortality at ages 3–15+.",
    reason = "SAM shares one fitted F state over ages 11–15; tinyAM estimates its age states independently."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Population-weighted F over ages 4–8, recalculated from each model's N and F.",
    reason = "The common-definition comparison uses matching ages and weighting; native SAM arithmetic Fbar is retained separately. The F-state representation remains an approximation."
  ),
  M = list(
    status = "matched",
    definition = "Fixed natural mortality of 0.2 per year over ages 3–15+.",
    reason = "This checks the shared fixed input; neither assessment estimates an M process."
  )
)

background <- c(
  "### Faroe saithe: 2025 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1961–2025 at ages 3–15+, with age 15 as the plus group. | Fit the same years and ages. | The source grid is retained. |",
  "| N and recruitment | Age-3 recruitment follows a random walk; age 3 has its own N variance and ages 4–15 share another. | Use IID N survival deviations and a separate random-walk recruitment process. | Initial-state integration and variance sharing are simplified. |",
  "| F | F states are separate at ages 3–10, with ages 11–15 shared; the age process is AR(1). Fbar is the arithmetic mean over ages 4–8. | Use a correlated random walk across ages and mean ages 4–8. | tinyAM estimates age states independently and uses an abundance-weighted Fbar. |",
  "| M and biology | M is fixed at 0.2 per year. Annual stock, catch, landings, and discard weights, maturity, landings fractions, and spawning fractions are supplied. | Use the same M, source biology, and fitted SAM 2025 stock/catch weights where the source files are missing. | The terminal weight estimates are fixed inputs in tinyAM; source values remain unchanged in the database. |",
  "| Catch | Catch numbers-at-age are in thousand fish; age 15 is the plus group. One observation variance is shared across ages. The 2025 relative precision weight is 0.5. | Retain catch-at-age and represent the relative precision by scaling the observation SD. | The source observation weighting is approximated by an SD multiplier. |",
  "| Surveys and q | SummerSurvey (1996–2025) and SpringSurvey (1994–2025) use ages 3–10 at times 0.65 and 0.225. Both have separate q at ages 3–8 and shared q at ages 9–10; survey residuals are AR(1) across ages. | Retain both age-specific indices, timing, q groups, and age-specific observation SDs. | tinyAM treats survey residuals as independent. |",
  "| Outputs | SAM reports N, F, recruitment, SSB, total biomass, arithmetic Fbar, q, and fitted observations. | Compare N/F at age and annual abundance, recruitment, biomass and SSB under shared translated biology, and population-weighted F over ages 4–8. | Native SAM aggregates are retained separately; derived comparison aggregates have no source uncertainty. |",
  "",
  "The fit uses tinyAM's normal initialization first. If numerical checks fail, one retry initializes N, F, and q from official estimates; these are starting values only, and tinyAM estimates them freely."
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
      F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 4:8),
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
        sd_form = ~0 + sd_block
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
        F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 4:8),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(
          sd_form = ~0 + sd_scale,
          sd_supplied = ~relative_sd
        ),
        index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)
      )
      starts <- tinyAM::make_par(model_dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages
      )
      starts$log_r0 <- log(source_N[1L, 1L])
      starts$log_r <- log(source_N[-1L, 1L])
      starts$log_n <- log(source_N[-1L, -1L, drop = FALSE])

      source_F <- source$outputs[
        source$outputs$type == "mortality" &
          source$outputs$measure == "fishing_mortality_at_age",
        , drop = FALSE
      ]
      f_matrix <- matrix(NA_real_, length(years), length(ages),
        dimnames = list(as.character(years), as.character(ages)))
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
