## Observations ----
years <- 1983:2022
ages <- 1:7
inputs <- source$inputs

collapse_biology <- function(inputs, type, measure) {
  rows <- inputs[inputs$type == type & inputs$measure == measure &
    !is.na(inputs$region), , drop = FALSE]
  rows$value <- as.numeric(rows$value)
  collapsed <- stats::aggregate(
    value ~ year + age, rows,
    function(x) mean(x[is.finite(x)])
  )
  rows <- rows[match(
    paste(collapsed$year, collapsed$age),
    paste(rows$year, rows$age)
  ), , drop = FALSE]
  rows$value <- collapsed$value
  rows$region <- NA_character_
  rows$notes <- paste(
    rows$notes,
    "Arithmetic mean of available Northwestern, Southern and Viking input values for a single-stock approximation."
  )
  rbind(inputs[!(inputs$type == type & inputs$measure == measure), ], rows)
}

inputs <- collapse_biology(inputs, "weight", "weight_at_age")
inputs <- collapse_biology(inputs, "maturity", "maturity_at_age")

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = "Survey_Q34",
  assumptions = source$assumptions
)

aggregate_rows <- function(rows, group_columns, weights = NULL) {
  key <- do.call(paste, c(lapply(group_columns, function(name) {
    value <- as.character(rows[[name]])
    value[is.na(value)] <- "<NA>"
    value
  }), sep = "\r"))
  groups <- split(seq_len(nrow(rows)), key)
  do.call(rbind, lapply(groups, function(index) {
    group <- rows[index, , drop = FALSE]
    result <- group[1, , drop = FALSE]
    values <- as.numeric(group$value)
    if (is.null(weights)) {
      result$value <- sum(values, na.rm = TRUE)
    } else {
      result$value <- stats::weighted.mean(values, weights[index], na.rm = TRUE)
    }
    result$region <- NA_character_
    for (name in intersect(c("se", "lwr", "upr"), names(result))) {
      result[[name]] <- NA_real_
    }
    result$notes <- paste(
      result$notes[1],
      "Substock estimates combined for comparison; uncertainty is not combined."
    )
    result
  }))
}

outputs <- source$outputs
population <- outputs[outputs$measure == "numbers_at_age" &
  outputs$type == "population", , drop = FALSE]
population_age <- suppressWarnings(as.integer(as.character(population$age)))
comparison_outputs <- list(aggregate_rows(population, c("year", "age")))
for (measure in c("SSB", "recruitment", "total_biomass", "predicted_catch")) {
  comparison_outputs[[length(comparison_outputs) + 1L]] <- aggregate_rows(
    outputs[outputs$measure == measure, , drop = FALSE], "year"
  )
}

for (measure in c("fishing_mortality_at_age", "natural_mortality_at_age")) {
  rows <- outputs[outputs$measure == measure, , drop = FALSE]
  n_rows <- population[match(
    paste(rows$year, rows$age, rows$region),
    paste(population$year, population$age, population$region)
  ), , drop = FALSE]
  comparison_outputs[[length(comparison_outputs) + 1L]] <- aggregate_rows(
    rows, c("year", "age"), as.numeric(n_rows$value)
  )
}

fbar <- outputs[outputs$measure == "Fbar", , drop = FALSE]
fbar_weights <- stats::aggregate(
  value ~ year + region,
  population[population_age %in% 2:4, , drop = FALSE], sum
)
fbar_weight <- fbar_weights$value[match(
  paste(fbar$year, fbar$region),
  paste(fbar_weights$year, fbar_weights$region)
)]
comparison_outputs[[length(comparison_outputs) + 1L]] <- aggregate_rows(
  fbar, "year", fbar_weight
)
comparison_outputs <- do.call(rbind, comparison_outputs)


## Background and comparisons ----

age_plus_group <- 7


comparison_scales <- c(
  N = 0.001, recruitment = 0.001, ssb = 0.001, F_bar = 1, abundance = 0.001, biomass = 0.001,
  biomass_at_age = 0.001
)


background <- c(
  "### Northern Shelf cod: accepted 2025 three-substock SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | Population model spans 1983–2025; catch ends in 2024 and supplied M observations end in 2022. | Fit 1983–2022. | Restrict the translation to years with a complete supplied M surface. |",
  "| Ages | Northwestern, Southern and Viking substocks are modelled at ages 1–7+, with recruitment at age 1. | One combined population at ages 1–7+, retaining the age-7 plus group. | Substock population structure is omitted. |",
  "| N | Recruitment follows a random walk, with separate abundance processes for the three substocks. | Random-walk recruitment, IID cohort-survival deviations with one shared SD, and freely estimated initial abundance. | Independent survival deviations replace the source's more elaborate population processes; substock-specific variation is omitted. |",
  "| F | Substock fishing mortality is linked through shared, scaled selectivity and temporally varying fishing pressure, with age-correlated F-process variation. | One age-specific correlated random-walk F process (cor_rw). | Cross-age correlation of annual F increments is approximated, but substock selectivity sharing and the complete source F structure are omitted. |",
  "| M | Natural mortality is estimated using a GMRF informed by supplied M observations. | Use the supplied year-age M values as fixed inputs. | The mortality GMRF and uncertainty in M are not estimated. |",
  "| Catch | Combined commercial catch-at-age is fitted alongside substock and quarterly composition observations. | Use combined catch-at-age with independent lognormal errors and one shared SD. | Auxiliary composition likelihoods, source catch-error correlation and variance structure are omitted. |",
  "| Survey indices | Seven fitted age-structured index streams include three substock Q1 indices, one mixed Q3+Q4 index and three substock recruitment indices. | Fit the mixed Q3+Q4 survey only, retaining its sampling time of 0.75 and supplied relative SDs. | Substock-specific abundance and recruitment indices are omitted because the translation uses one combined population. |",
  "| Survey catchability | Catchability is estimated within the accepted multistock observation model; numerical q estimates have not been independently recovered. | Use mono(q_key) to estimate non-decreasing catchability across ordered survey age groups, allowing plateaus. | A monotone age pattern is imposed to avoid estimating dome-shaped catchability from the reduced data inputs. This is a translation assumption, not a recovered SAM constraint. |",
  "| Weights and maturity | Annual weight and maturity inputs are specific to each substock. | Use arithmetic means of available substock inputs for each year and age. | Population-weighted combined biological inputs are unavailable. |",
  "| Comparison outputs | Accepted outputs are reported separately for Northwestern, Southern and Viking substocks. | Sum numbers, recruitment and biomass; calculate abundance-weighted F and M by age and abundance-weighted substock Fbar over ages 2–4. | Comparison quantities are approximate and their uncertainty is not combined. |",
  "",
  "This translation tests whether tinyAM can reproduce broad Northern Shelf cod dynamics using a single combined population and the mixed Q3+Q4 survey.",
  "",
  "The survey catchability formula mono(q_key) constrains q to be non-decreasing across the ordered age groups. It permits flat sections but prevents a dome-shaped catchability curve. This is an intentional simplifying assumption, not a constraint attributed to the accepted SAM assessment.",
  "",
  "The revised fit also uses IID cohort-survival deviations and correlated F random walks. The accepted assessment's substock structure, auxiliary composition likelihoods and biological GMRFs are not reproduced.",
  "",
  "Biological inputs are averaged across substocks without abundance weighting. Accepted outputs are aggregated or abundance-weighted for comparison, so resulting differences should not be interpreted as errors in a like-for-like reproduction of the accepted multistock model."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~ mono(q_key), sd_form = ~1, sd_supplied = ~relative_sd),
    silent = silent
  )
}
