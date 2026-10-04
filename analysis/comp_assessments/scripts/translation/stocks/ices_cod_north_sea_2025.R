translate_stock <- function(source) {
  years <- 1983:2022
  ages <- 1:7
  inputs <- source$inputs

  collapse_biology <- function(inputs, type, measure) {
    rows <- inputs[inputs$type == type & inputs$measure == measure &
                     !is.na(inputs$region), , drop = FALSE]
    rows$value <- as.numeric(rows$value)
    collapsed <- stats::aggregate(value ~ year + age, rows,
                                  function(x) mean(x[is.finite(x)]))
    rows <- rows[match(paste(collapsed$year, collapsed$age),
                       paste(rows$year, rows$age)), , drop = FALSE]
    rows$value <- collapsed$value
    rows$region <- NA_character_
    rows$notes <- paste(rows$notes,
                        "Arithmetic mean of available Northwestern, Southern and Viking input values for a single-stock approximation.")
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
      result$notes <- paste(result$notes[1],
                            "Substock estimates combined for comparison; uncertainty is not combined.")
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
    n_rows <- population[match(paste(rows$year, rows$age, rows$region),
                               paste(population$year, population$age, population$region)), , drop = FALSE]
    comparison_outputs[[length(comparison_outputs) + 1L]] <- aggregate_rows(
      rows, c("year", "age"), as.numeric(n_rows$value)
    )
  }

  fbar <- outputs[outputs$measure == "Fbar", , drop = FALSE]
  fbar_weights <- stats::aggregate(value ~ year + region,
    population[population_age %in% 2:4, , drop = FALSE], sum)
  fbar_weight <- fbar_weights$value[match(paste(fbar$year, fbar$region),
                                          paste(fbar_weights$year, fbar_weights$region))]
  comparison_outputs[[length(comparison_outputs) + 1L]] <- aggregate_rows(
    fbar, "year", fbar_weight
  )
  comparison_outputs <- do.call(rbind, comparison_outputs)

  list(
    years = years,
    ages = ages,
    age_plus_group = 7,
    obs = obs,
    comparison_outputs = comparison_outputs,
    comparison_scales = c(
      N = 1e-3, recruitment = 1e-3, ssb = 1e-3, F_bar = 1,
      abundance = 1e-3, biomass = 1e-3, biomass_at_age = 1e-3
    ),
    settings = list(
      N_settings = list(process = "rw", init = "free"),
      F_settings = list(process = "rw", mu_form = NULL),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 0 + q_key,
                            sd_form = ~ 1,
                            sd_supplied = ~ relative_sd,
                            fill_missing = FALSE)
    ),
    background = c(
      "### Northern Shelf cod: accepted 2025 three-substock SAM assessment",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | The accepted population model spans 1983–2025; catch ends in 2024 and supplied M observations end in 2022. | Fit 1983–2022. | This uses only years with a complete supplied M surface and compares outputs over that common period. |",
      "| Ages | Three substocks are modeled at ages 1–7+, with age 7 as the plus group. | Fit one combined population at ages 1–7, retaining age 7 as the plus group. | tinyAM has no substock population structure. |",
      "| N | Recruitment follows a random walk and each substock has its own abundance process. | One random-walk cohort process with free initial abundance. | Substock-specific processes and their covariance are not represented. |",
      "| F | Substock F uses scaled selectivity and age-correlated random-walk increments. | One age-specific random-walk F process. | tinyAM cannot share selectivity or correlate F increments across ages. |",
      "| M | M is estimated as a GMRF process from supplied M observations. | Treat the supplied year-age M values as fixed. | tinyAM does not fit the source biological GMRF. |",
      "| Catch | Combined commercial catch-at-age is used with additional substock and quarterly compositions. | Use combined catch-at-age with independent lognormal errors. | tinyAM does not include the auxiliary composition likelihoods or age-correlated catch errors. |",
      "| Index | Seven age-structured streams include three substock Q1 surveys, one mixed Q3+Q4 survey and three substock recruitment indices. | Use the mixed Q3+Q4 survey only, retaining its supplied relative SDs and 0.75 timing. | This is the only selected index representing the combined population; substock and recruitment indices are omitted. |",
      "| Weights and maturity | Annual weight and maturity inputs are specific to each substock. | Use the arithmetic mean of available substock inputs for each year and age. | No all-year population composition weights are available to form a combined surface. |",
      "| Comparison outputs | The database stores accepted outputs separately for all three substocks. | Sum abundance, recruitment and biomass; calculate N-weighted F and M by age and N-weighted Fbar over ages 2–4. | The single-stock comparison is approximate, and uncertainty is not combined. |",
      "",
      "This fit tests whether tinyAM can describe a simplified combined stock. Its biological inputs are unweighted across substocks, so the resulting SSB, abundance and mortality should not be interpreted as a like-for-like replication of the accepted multistock assessment."
    )
  )
}
