## Observations ----

years <- 1981:2021
ages <- 1:7
assessment_id <- source$assessment$assessment_id[[1L]]
surveys <- c("WCIBTS.Q1", "IRGFS.Q4", "WCIBTS.Q4", "SCO.Q1", "SCO.Q4")
obs <- database_to_tam_obs(
  assessment_id, source$inputs, years = years, ages = ages,
  surveys = surveys, assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0

precision <- source$inputs[source$inputs$measure == "relative_precision_weight", ]
relative_sd <- function(x, type) {
  weights <- precision[precision$type == type, ]
  key <- function(x) paste(x$year, x$age,
    if (type == "index") x$survey else "", sep = "/")
  at <- match(key(x), key(weights))
  value <- weights$value[at]
  value[is.na(at)] <- 1
  if (any(!is.finite(value) | value <= 0)) cli::cli_abort("Invalid source relative precision.")
  1 / sqrt(value)
}
obs$catch$relative_sd <- relative_sd(obs$catch, "catch")
obs$catch$sd_block <- factor(ifelse(obs$catch$age == 1, "age1",
  ifelse(obs$catch$age == 7, "age7plus", "ages2_6")))
obs$index$relative_sd <- relative_sd(obs$index, "index")

# Native age-to-parameter keys: equal keys share q; surveys never share q.
q_keys <- list(
  WCIBTS.Q1 = 0:5,
  IRGFS.Q4 = c(0, 1, 1),
  WCIBTS.Q4 = c(0, 1, 2, 2),
  SCO.Q1 = c(0:3, 4, 4),
  SCO.Q4 = c(0:3, 4, 4)
)
obs$index$q_key <- factor(vapply(seq_len(nrow(obs$index)), function(i) {
  survey <- obs$index$survey[i]
  paste(survey, q_keys[[survey]][obs$index$age[i]], sep = ":")
}, character(1)))

## Comparison ----

age_plus_group <- 7L
comparison_outputs <- source$outputs
index_predictions <- comparison_outputs$type == "index" &
  comparison_outputs$measure == "predicted_index"
q_estimates <- comparison_outputs$type == "catchability" & comparison_outputs$measure == "q"
columns <- c("value", "se", "lwr", "upr")
comparison_outputs[index_predictions | q_estimates, columns] <-
  comparison_outputs[index_predictions | q_estimates, columns] * 0.001
irish <- comparison_outputs$survey %in% "IRGFS.Q4"
comparison_outputs$unit[index_predictions] <- ifelse(irish[index_predictions],
  "thousand fish per minute", "thousand fish per hour")
comparison_outputs$unit[q_estimates] <- ifelse(irish[q_estimates], "per minute", "per hour")
comparison_aggregates <- c("abundance", "biomass", "ssb", "recruitment", "F_bar")
comparison_definitions <- list(
  abundance = list(status = "matched", definition = "Beginning-of-year N summed over ages 1-7+.",
    reason = "Matching calendar years, ages and units; derived official totals have no uncertainty."),
  recruitment = list(status = "matched", definition = "Age-1 recruitment in the same calendar year.",
    reason = "Recruitment age and units match."),
  N = list(status = "matched", definition = "Beginning-of-year N at ages 1-7+.",
    reason = "Matching population ages and timing; process and initial-state integration differ."),
  F = list(status = "approximate", definition = "Annual F at ages 1-7+.",
    reason = "Official F shares the ages 4-7 state; tinyAM represents each age separately."),
  F_bar = list(status = "approximate", definition = "Population-weighted F over ages 2-5 in both models.",
    reason = "Native arithmetic SAM Fbar is retained separately; the F-state representation differs."),
  M = list(status = "matched", definition = "Supplied annual Lorenzen M-at-age.",
    reason = "Input-consistency check; M is fixed, not estimated, in both models."),
  M_bar = list(status = "matched", definition = "Population-weighted mean supplied M over matching ages.",
    reason = "Differences reflect population weights, not estimated M."),
  biomass = list(status = "matched", definition = "Biomass at ages 1-7+ under shared source stock weights.",
    reason = "Source weights and beginning-of-year definition match."),
  ssb = list(status = "matched", definition = "Mature biomass at ages 1-7+ with source weights and maturity.",
    reason = "Both models use zero pre-spawning F and M fractions."))

background <- c(
  "### West of Scotland cod: final standalone 2022 assessment",
  "", print_sources(source$assessment), "",
  "This is the historical Division 6.a assessment. It was merged into the broader Northern Shelf assessment in 2023 and is not an additional independent contemporary stock.",
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1981-2021; ages 1-7+. No 2022 intermediate-year survey. | Same historical grid and plus group. | None. |",
  "| N and recruitment | Random-walk age-1 recruitment; independent older-age survival errors share one SD. Separate initial latent N at every age. | Random-walk recruitment, IID N process and free initial abundance. | Initial abundance is fitted directly rather than integrated as SAM random effects. |",
  "| F | Correlated temporal random-walk increments; ages 1, 2 and 3 distinct, ages 4-7 share one state. Age 1 has a separate process SD. | Correlated random-walk F at each age, separate age-1 SD; mean ages 2-5. | The shared older F state is approximated by separate correlated states. |",
  "| M | Fixed annual age-specific Lorenzen M, based on the source smoothed stock weights. | Supply the same complete M surface. | No estimated M process. |",
  "| Catch | Lognormal residual catch includes landings and discards; SD groups are age 1, ages 2-6 and age 7+. SD doubles from 2006. | Same observations, SD groups and relative precision. | SAM also estimates annual catch-scaling factors in 1995-2006; tinyAM omits that correction. |",
  "| Indices and q | Five surveys with independent lognormal residuals, survey-specific SDs, age-sharing q keys and relative precision for the two current Scottish surveys. | Retain all five series, q sharing, SD sharing, precision weights and timing. | Scottish rates are per hour; Irish rates remain per minute. |",
  "| Biology | Annual stock/catch weights, maturity and fixed M; zero pre-spawning fractions. | Retain source biology and calculate mature biomass at time zero. | Source biological smoothing is treated as fixed. |",
  "| Comparisons | N, F, SSB, recruitment, biomass and arithmetic Fbar. | Compare matched ages and units; recalculate population-weighted Fbar for both models. | Native outputs remain available separately; derived totals have no invented uncertainty. |",
  "",
  "Catch and index observations are unchanged. In particular, the fitted SAM catch-scaling factors are never used to adjust the tinyAM observations. The source's missing 2021 catches at ages 1 and 2 remain missing, and zeros are omitted by the lognormal observation model.",
  "",
  "The model uses normal package starting values. Numerical convergence and agreement with the official assessment are assessed separately."
)

## Model ----

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "iid", init = "free", sd_form = ~1),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mu_form = NULL,
      sd_form = ~I(age == 1), mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_block, sd_supplied = ~relative_sd),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey,
      sd_supplied = ~relative_sd),
    silent = silent
  )
}
