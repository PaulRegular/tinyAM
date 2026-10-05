translate_stock <- function(source) {
  years <- 1971:2018
  ages <- 2:12
  rv <- "DFO September RV survey"
  mobile <- "Mobile Sentinel August survey"

  rv_weights <- source$inputs[
    source$inputs$type == "weight" &
      source$inputs$measure == "weight_at_age" &
      source$inputs$survey == rv &
      source$inputs$year %in% years &
      source$inputs$age %in% 2:11,
    , drop = FALSE
  ]
  missing_years <- setdiff(years, unique(rv_weights$year))
  if (!all(missing_years %in% c(1980, 1985))) {
    cli::cli_abort("Only the documented 1980 and 1985 RV weight gaps can be interpolated.")
  }
  interpolated_weights <- do.call(rbind, lapply(missing_years, function(year) {
    rows <- rv_weights[rv_weights$year == year - 1, , drop = FALSE]
    next_year <- rv_weights[rv_weights$year == year + 1, , drop = FALSE]
    if (!setequal(rows$age, 2:11) || !setequal(next_year$age, 2:11)) {
      cli::cli_abort("Adjacent RV weight years must contain ages 2-11 before interpolation.")
    }
    rows <- rows[match(2:11, rows$age), , drop = FALSE]
    next_year <- next_year[match(2:11, next_year$age), , drop = FALSE]
    rows$year <- year
    rows$value <- (rows$value + next_year$value) / 2
    rows$source_type <- "translation_assumption"
    rows$transformation <- "Linear interpolation of RV weights for this fit only."
    rows$notes <- "No official RV weight-at-age was tabulated for this year."
    rows
  }))
  rv_weights <- rbind(rv_weights, interpolated_weights)
  plus_weight <- rv_weights[rv_weights$age == 11, , drop = FALSE]
  plus_weight$age <- 12
  plus_weight$source_type <- "translation_assumption"
  plus_weight$transformation <- "Age-11 RV weight carried forward to the tinyAM 12+ group."
  plus_weight$notes <- "Fit-only proxy; RV weights are not reported for ages 12+."

  landings_at_age <- !is.na(source$inputs$type) &
    source$inputs$type == "catch" &
    !is.na(source$inputs$measure) &
    source$inputs$measure == "landings_numbers_at_age"
  if (!any(landings_at_age)) {
    cli::cli_abort("The source-reported landings-at-age series is missing.")
  }
  catch_proxy <- source$inputs[landings_at_age, , drop = FALSE]
  catch_proxy$measure <- "numbers_at_age"
  catch_proxy$source_type <- "translation_assumption"
  catch_proxy$transformation <- paste(
    "Used source landings-at-age as a fit-only proxy for catch-at-age;",
    "no scaling to total catch was possible because the fitted age",
    "composition was not tabulated."
  )
  catch_proxy$notes <- paste(
    "This is not the accepted model's catch-at-age input. The source",
    "model fits total catch biomass and age proportions separately."
  )
  inputs <- rbind(
    source$inputs[!landings_at_age, , drop = FALSE], catch_proxy,
    interpolated_weights, plus_weight
  )

  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    inputs,
    years = years,
    ages = ages,
    weight_survey = rv,
    sampling_times = c(setNames(0.75, rv), setNames(0.625, mobile)),
    surveys = c(rv, mobile),
    exclude_index_years = setNames(list(c(1980, 1985, 2003)), rv),
    assumptions = source$assumptions
  )
  obs$weight$M_prior_mean <- ifelse(obs$weight$age <= 4, 0.65, 0.15)
  comparison_outputs <- source$outputs[
    source$outputs$measure %in% c(
      "numbers_at_age", "fishing_mortality_at_age", "SSB", "recruitment"
    ), , drop = FALSE
  ]
  terminal_m <- source$outputs[
    source$outputs$measure == "natural_mortality_at_age", , drop = FALSE
  ]
  if (nrow(terminal_m)) {
    terminal_m <- do.call(rbind, lapply(seq_len(nrow(terminal_m)), function(i) {
      row <- terminal_m[i, , drop = FALSE]
      group_ages <- switch(as.character(row$age_group),
                           "5-8" = 5:8, "9+" = 9:12, integer())
      if (!length(group_ages)) return(NULL)
      row <- row[rep(1L, length(group_ages)), , drop = FALSE]
      row$age <- group_ages
      row$age_group <- NA_character_
      row$notes <- paste(row$notes,
                         "The reported group value is expanded across its ages for comparison.")
      row
    }))
    comparison_outputs <- rbind(comparison_outputs, terminal_m)
  }

  list(
    years = years,
    ages = ages,
    age_plus_group = 12,
    obs = obs,
    comparison_outputs = comparison_outputs,
    comparison_scales = c(N = 1e-3, recruitment = 1e-3, ssb = 1e-3),
    settings = list(
      N_settings = list(process = "off", init = "exp"),
      F_settings = list(process = "ar1", mu_form = NULL),
      M_settings = list(
        process = "rw",
        mu_form = NULL,
        mu_supplied = ~ M_prior_mean,
        age_breaks = c(2, 5, 9, 12),
        first_dev_year = 1971L
      ),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(
        q_form = ~ 0 + q_key,
        sd_form = ~ 0 + survey,
        fill_missing = FALSE
      )
    ),
    warm_start_settings = list(M_settings = list(process = "iid")),
    background = c(
      "### Southern Gulf cod: accepted detailed assessment to 2018",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | The SCA model covers 1950-2018; source landings-at-age data begin in 1971. RV weights are missing for 1980 and 1985. | Fit 1971-2018. Linearly interpolate RV weights for 1980 and 1985 for this fit only; exclude RV indices in those years and the anomalous 2003 index. | This retains the full landings-at-age period without adding reconstructed values to the source database. |",
      "| Ages | The population model uses ages 2-12+, while survey age compositions cover ages 2-11. | Use ages 2-12, with age 12 as the plus group. | The reported RV age-11 weight is carried to 12+ as a fit-only weight proxy. |",
      "| N | Recruitment enters at age 2, depends on SSB two years earlier, and has autocorrelated variation; initial cohorts are reconstructed from recruitment. | Use exponential initial abundance, deterministic cohort survival, and tinyAM's recruitment process. | tinyAM does not reproduce the source stock-recruit relationship, recruitment autocorrelation, or initial-cohort estimation. |",
      "| F | The source estimates fully recruited F and period-specific logistic selectivity. | Use an age- and year-correlated AR1 F process. | This is a simpler representation of changing fishing mortality and selectivity. |",
      "| M | The source estimates log-M random walks for ages 2-4, 5-8, and 9+. Each group has one estimated initial M level through 1971, with prior means 0.65, 0.15, and 0.15 and SD 0.05; log-M increments begin in 1972 with SD fixed at 0.075. | Use the same age groups, with the first tinyAM M state in 1971 and random-walk increments from 1972. The prior means initialize those states; the final fit is initialized from a converged IID-M fit, whose estimates are starting values only. | tinyAM cannot apply the source prior to initial M or fix the random-walk SD at 0.075, so it estimates the increment SD and leaves each initial M state unpenalized. |",
      "| Catch | The source fits annual catch biomass and proportions-at-age for ages 2-12+. The database contains landed numbers-at-age for ages 3-12+ from 1971 onward, not the fitted catch composition. | Use the landings-at-age series as an explicit fit-only proxy for catch-at-age in 1971-2018; age 2 remains missing. | The original age proportions are not tabulated, so landings cannot be scaled to the source's total catch. tinyAM then uses a lognormal age-specific observation model rather than the source's total-plus-composition likelihood. |",
      "| Index | The source uses RV, mobile sentinel, and longline indices. RV 2003 is excluded by the assessment; longline combines July-October observations. | Use RV and mobile sentinel age-specific indices, with sampling times 0.75 and 0.625. Exclude RV 2003 and omit longline. | The timing values are seasonal approximations; one within-year time is not supported for the longline series. |",
      "| Weights and maturity | Survey weights and year-varying maturity are reported; the source's RV weights are available by age 2-11. | Use RV weights for population biomass and source survey weights for index reconstruction. Retain annual maturity. | Age-12+ RV weight is unavailable and uses the age-11 proxy noted above. The source does not specify maturity by sex. |",
      "| Comparison | Tables 21-23 report MLE SSB and age-specific N and F, plus recruitment; the report gives terminal M for ages 5-8 and 9+. Other estimates are posterior medians. | Compare SSB, N, F, recruitment, and the two reported terminal M groups over 1971-2018, scaling numbers and biomass to tinyAM units. | The terminal M groups are repeated across their constituent ages for the dashboard; age-specific M was not published, and the published estimates are not all on the same uncertainty basis. |",
      "",
      "This is a limited comparison of model behaviour, not a recreation of the accepted likelihood. The 2024 advice update is recorded separately; no 2024-run inputs or outputs are substituted for this 2019 detailed assessment."
    )
  )
}
