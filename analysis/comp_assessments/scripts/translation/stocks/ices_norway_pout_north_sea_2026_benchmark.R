translate_stock <- function(source) {
  years <- 1984:2024
  ages <- 0:3
  id <- source$assessment$assessment_id[[1]]
  inputs <- source$inputs
  inputs$year <- suppressWarnings(as.integer(inputs$year))
  inputs$age <- suppressWarnings(as.integer(inputs$age))
  inputs$value <- suppressWarnings(as.numeric(inputs$value))
  inputs$quarter <- suppressWarnings(as.integer(sub("^Q", "", inputs$season)))

  collapse_quarters <- function(rows, quarters_for_age, reduce, label) {
    groups <- split(seq_len(nrow(rows)), paste(rows$year, rows$age, sep = "\034"))
    do.call(rbind, lapply(groups, function(i) {
      z <- rows[i, , drop = FALSE]
      wanted <- quarters_for_age(z$age[[1]])
      z <- z[z$quarter %in% wanted, , drop = FALSE]
      if (!setequal(z$quarter, wanted)) {
        stop(label, " is missing a required quarter for year ", z$year[[1]],
             ", age ", z$age[[1]], ".")
      }
      z$value <- reduce(z$value)
      z <- z[1L, , drop = FALSE]
      z$season <- ""
      z$quarter <- NA_integer_
      z$sampling_time <- NA_real_
      z
    }))
  }

  selected_catch <- inputs[
    inputs$type == "catch" & inputs$measure == "numbers_at_age" &
      inputs$year %in% years, , drop = FALSE
  ]
  annual_catch <- collapse_quarters(
    selected_catch, function(age) 1:4, sum, "Quarterly catch"
  )
  annual_catch$transformation <- paste(
    annual_catch$transformation,
    "Quarterly catch numbers summed to calendar-year totals for tinyAM."
  )

  selected_m <- inputs[
    inputs$type == "M" & inputs$measure == "natural_mortality_at_age" &
      inputs$year %in% years, , drop = FALSE
  ]
  annual_m <- collapse_quarters(
    selected_m, function(age) if (age == 0L) 3:4 else 1:4,
    function(value) sum(value) * 0.25, "Quarterly M"
  )
  annual_m$transformation <- paste(
    annual_m$transformation,
    "Converted quarterly rates to annual-equivalent cumulative hazard; age 0 is exposed from Q3."
  )
  annual_m$notes <- paste(
    annual_m$notes,
    "Age 0 uses Q3-Q4 exposure only; older ages use all four quarters."
  )

  weight <- inputs[
    inputs$type == "weight" & inputs$measure == "weight_at_age" &
      is.na(inputs$year), , drop = FALSE
  ]
  weight <- weight[weight$quarter == ifelse(weight$age == 0L, 3L, 1L), , drop = FALSE]
  if (!setequal(weight$age, ages) || anyDuplicated(weight$age)) {
    stop("The accepted stock weights do not cover ages 0-3 once at the selected year boundary.")
  }
  weight$season <- ""
  weight$quarter <- NA_integer_
  weight$transformation <- paste(
    weight$transformation,
    "Selected Q3 weight for age-0 recruitment and Q1 weights for older ages to match the annual state boundary."
  )

  index <- inputs[
    inputs$type == "index" & inputs$measure == "numbers_at_age" &
      inputs$year %in% years, , drop = FALSE
  ]
  maturity <- inputs[
    inputs$type == "maturity" & inputs$measure == "maturity_at_age", , drop = FALSE
  ]
  annual_inputs <- rbind(
    annual_catch,
    index,
    weight,
    annual_m,
    maturity
  )
  annual_inputs$quarter <- NULL

  surveys <- c(
    "IBTS Q1 (north of 57N)", "EGFS Q3",
    "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)",
    "IBTS other countries Q3", "SGFS Q3 age 0 (2013 onward)"
  )
  sampling_times <- c(
    "IBTS Q1 (north of 57N)" = 0.125,
    "EGFS Q3" = 0.625,
    "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)" = 0.625,
    "IBTS other countries Q3" = 0.625,
    "SGFS Q3 age 0 (2013 onward)" = 0.625
  )
  obs <- database_to_tam_obs(
    id, annual_inputs, years = years, ages = ages,
    surveys = surveys, sampling_times = sampling_times,
    assumptions = source$assumptions
  )

  catch_group <- ifelse(obs$catch$age == 0, "age0",
                        ifelse(obs$catch$age <= 2, "age1_2", "age3plus"))
  obs$catch$sd_block <- factor(catch_group,
                               levels = c("age0", "age1_2", "age3plus"))
  survey_code <- c(
    "IBTS Q1 (north of 57N)" = "ibts_q1",
    "EGFS Q3" = "egfs_q3",
    "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)" = "sgfs_q3",
    "IBTS other countries Q3" = "ibts_other_q3",
    "SGFS Q3 age 0 (2013 onward)" = "sgfs_age0_q3"
  )
  index <- obs$index
  survey_name <- as.character(index$survey)
  code <- unname(survey_code[survey_name])
  age_key <- as.character(index$age)
  shared_q <- survey_name %in% c("EGFS Q3", surveys[[3]]) & index$age %in% c(2L, 3L)
  age_key[shared_q] <- "2_3"
  index$q_key <- factor(paste(code, age_key, sep = "_age"))
  index$sd_block <- paste0(code, "_age", index$age)
  sgfs_age0 <- survey_name %in% c(surveys[[3]], surveys[[5]]) & index$age == 0L
  index$sd_block[sgfs_age0] <- "sgfs_age0"
  sgfs_shared <- survey_name == surveys[[3]] & index$age %in% c(1L, 2L)
  index$sd_block[sgfs_shared] <- "sgfs_age1_2"
  obs$index <- index

  settings <- list(
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~ factor(age), mean_ages = 1:2),
    M_settings = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~ 0 + q_key,
                          sd_form = ~ 0 + sd_block,
                          fill_missing = FALSE)
  )

  quarter_cube <- function(rows, multiplier = 1) {
    out <- array(NA_real_, c(length(years), length(ages), 4L),
                 dimnames = list(as.character(years), as.character(ages),
                                 paste0("Q", 1:4)))
    rows$year <- as.integer(rows$year)
    rows$age <- as.integer(rows$age)
    rows$value <- as.numeric(rows$value)
    quarter <- as.integer(sub("^Q", "", rows$season))
    i <- cbind(match(rows$year, years), match(rows$age, ages), quarter)
    out[i] <- rows$value * multiplier
    if (any(!is.finite(out))) stop("An accepted quarterly surface is incomplete.")
    out
  }

  source_outputs <- source$outputs
  source_outputs$year <- as.integer(source_outputs$year)
  source_outputs$age <- as.integer(source_outputs$age)
  source_outputs$value <- as.numeric(source_outputs$value)
  source_inputs <- source$inputs
  n_rows <- source_outputs[
    source_outputs$measure == "numbers_at_age" & source_outputs$year %in% years, , drop = FALSE
  ]
  f_rows <- source_outputs[
    source_outputs$measure == "fishing_mortality_at_age" &
      source_outputs$year %in% years, , drop = FALSE
  ]
  m_rows <- source_inputs[
    source_inputs$type == "M" & source_inputs$measure == "natural_mortality_at_age" &
      source_inputs$year %in% years, , drop = FALSE
  ]
  n_cube <- quarter_cube(n_rows, multiplier = 1e6)
  f_cube <- quarter_cube(f_rows)
  m_cube <- quarter_cube(m_rows)

  source_n <- n_cube[, , "Q1"]
  source_n[, "0"] <- n_cube[, "0", "Q3"]
  annual_f <- annual_m_surface <- matrix(
    NA_real_, length(years), length(ages),
    dimnames = list(as.character(years), as.character(ages))
  )
  for (age in ages) {
    quarters <- if (age == 0L) 3:4 else 1:4
    annual_f[, as.character(age)] <- rowSums(f_cube[, as.character(age), quarters,
                                                       drop = FALSE]) * 0.25
    annual_m_surface[, as.character(age)] <- rowSums(
      m_cube[, as.character(age), quarters, drop = FALSE]
    ) * 0.25
  }
  if (any(!is.finite(source_n)) || any(source_n <= 0) ||
      any(!is.finite(annual_f)) || any(annual_f <= 0) ||
      any(!is.finite(annual_m_surface)) || any(annual_m_surface <= 0)) {
    stop("The annualized accepted population or mortality surface is invalid.")
  }

  dat <- do.call(tinyAM::prepare_tam, c(
    list(data = obs, years = years, ages = ages), settings
  ))
  start_par <- tinyAM::make_par(dat)
  start_par$log_r0 <- log(source_n[1L, "0"])
  start_par$log_r <- setNames(log(source_n[-1L, "0"]), as.character(years[-1L]))
  start_par$log_n <- log(source_n[-1L, as.character(ages[-1L]), drop = FALSE])
  start_par$log_f <- log(annual_f)

  index <- obs$index
  quarter <- ifelse(index$samp_time < 0.5, 1L, 3L)
  at <- cbind(match(index$year, years), match(index$age, ages), quarter)
  n_at_survey <- n_cube[at] * exp(-(f_cube[at] + m_cube[at]) * 0.125)
  q_design <- stats::model.matrix(~ 0 + q_key, data = index)
  q_start <- vapply(seq_len(ncol(q_design)), function(j) {
    selected <- q_design[, j] == 1
    log(stats::median(index$obs[selected] / n_at_survey[selected], na.rm = TRUE))
  }, numeric(1))
  names(q_start) <- colnames(q_design)
  if (!setequal(names(q_start), names(start_par$log_q)) ||
      any(!is.finite(q_start))) {
    stop("Could not create finite source-based q starting values.")
  }
  start_par$log_q[] <- q_start[names(start_par$log_q)]

  q1_n <- n_cube[, , "Q1"]
  fbar <- vapply(seq_along(years), function(i) {
    stats::weighted.mean(annual_f[i, c("1", "2")], q1_n[i, c("1", "2")])
  }, numeric(1))
  common_fbar <- source_outputs[0, , drop = FALSE]
  for (i in seq_along(years)) {
    row <- source_outputs[source_outputs$measure == "Fbar", , drop = FALSE][1L, , drop = FALSE]
    row$type <- "mortality"
    row$measure <- "Fbar"
    row$year <- years[[i]]
    row$season <- ""
    row$age <- NA_integer_
    row$age_group <- ""
    row$value <- fbar[[i]]
    row$se <- row$lwr <- row$upr <- NA_real_
    row$unit <- "per year"
    row$source_type <- "derived_common_definition"
    row$source_reference <- paste(unique(f_rows$source_reference), collapse = "; ")
    row$notes <- "N-weighted annual-equivalent F over ages 1-2, matching tinyAM's Fbar definition."
    common_fbar <- rbind(common_fbar, row)
  }

  comparison_outputs <- source_outputs[
    source_outputs$year %in% years &
      ((source_outputs$measure == "numbers_at_age" &
          source_outputs$season == "Q1" & source_outputs$age %in% 1:3) |
       (source_outputs$measure == "SSB" & source_outputs$season == "Q1") |
       source_outputs$measure == "recruitment"), , drop = FALSE
  ]
  comparison_outputs$season <- ""
  comparison_outputs$notes[comparison_outputs$measure == "numbers_at_age"] <- paste(
    comparison_outputs$notes[comparison_outputs$measure == "numbers_at_age"],
    "Q1 abundance retained for ages 1-3; age-0 abundance is omitted because recruitment enters in Q3 and is compared separately."
  )
  comparison_outputs$notes[comparison_outputs$measure == "SSB"] <- paste(
    comparison_outputs$notes[comparison_outputs$measure == "SSB"],
    "Q1 SSB uses the same start-of-year N, weight, and maturity definition as tinyAM."
  )
  comparison_f <- f_rows[0, , drop = FALSE]
  for (year in years) {
    for (age in ages) {
      row <- f_rows[f_rows$year == year & f_rows$age == age &
                      f_rows$season == "Q1", , drop = FALSE][1L, , drop = FALSE]
      row$season <- ""
      row$value <- annual_f[as.character(year), as.character(age)]
      row$se <- row$lwr <- row$upr <- NA_real_
      row$source_type <- "derived_common_definition"
      row$notes <- if (age == 0L) {
        "Annual-equivalent F from Q3-Q4, when age-0 recruits are present."
      } else {
        "Annual-equivalent F from the four quarterly rates."
      }
      comparison_f <- rbind(comparison_f, row)
    }
  }
  common_biology <- source_outputs[0, , drop = FALSE]
  for (measure in c("biomass_at_age", "mature_biomass_at_age")) {
    for (i in seq_along(years)) {
      for (j in seq_along(ages)) {
        age <- ages[[j]]
        row <- source_outputs[source_outputs$measure == "SSB" &
                                source_outputs$season == "Q1", , drop = FALSE][1L, , drop = FALSE]
        row$type <- "biomass"
        row$measure <- measure
        row$year <- years[[i]]
        row$season <- ""
        row$age <- age
        row$age_group <- ""
        weight <- obs$weight$obs[match(
          paste(years[[i]], age), paste(obs$weight$year, obs$weight$age)
        )]
        maturity <- obs$maturity$obs[match(
          paste(years[[i]], age), paste(obs$maturity$year, obs$maturity$age)
        )]
        biology <- if (measure == "mature_biomass_at_age") maturity else 1
        row$value <- source_n[i, j] * weight * biology / 1000
        row$se <- row$lwr <- row$upr <- NA_real_
        row$unit <- "tonnes"
        row$source_type <- "derived_common_definition"
        row$source_reference <- paste(unique(c(
          n_rows$source_reference, inputs$source_reference[
            inputs$type %in% c("weight", "maturity")
          ]
        )), collapse = "; ")
        row$notes <- paste(
          "Common-definition", measure,
          "from accepted age-specific abundance, stock weights, and maturity; age 0 N is Q3 and older ages are Q1."
        )
        common_biology <- rbind(common_biology, row)
      }
    }
  }
  comparison_outputs <- rbind(
    comparison_outputs, comparison_f, common_fbar, common_biology
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 3,
    obs = obs,
    start_par = start_par,
    settings = settings,
    comparison_outputs = comparison_outputs,
    comparison_scales = c(
      N = 1e-6, recruitment = 1e-6, ssb = 1e-3,
      F = 1, M = 1, F_bar = 1
    ),
    background = c(
      "### North Sea Norway pout: accepted 2026 SESAM benchmark",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Main simplification |",
      "|---|---|---|---|",
      "| Years | Quarterly inputs and states cover 1984-2025, but the fitted catch series ends in Q3 2025. | Fit 1984-2024, the last complete catch year. | The partial 2025 catch total is not treated as an annual observation. |",
      "| Ages | Ages 0-3+, with recruitment entering in Q3. | Annual ages 0-3, with age 0 represented at recruitment. | Recruitment timing is annualized; source Q3 timing is retained for survey predictions. |",
      "| N | Quarterly age-specific abundance; process variance is shared for ages 0-2 and separate for age 3+. | IID annual abundance process with exponential initial abundance. | The quarterly process and its age-specific variance sharing are not reproduced. |",
      "| F | Quarterly log-F random walk with AR(1) dependence across ages; Fbar is ages 1-2. | Annual AR(1) F process around estimated age-specific means and Fbar ages 1-2. | tinyAM applies AR(1) in both year and age and cannot combine the source random walk with its age correlation; quarterly variation and the big-jump adjustment are approximated. |",
      "| M | Fixed quarterly age-specific M from North Sea SMS; the 2022 pattern is carried through 2025. | Supply annual-equivalent M as fixed mortality. | Older ages use all four quarters; age 0 uses Q3-Q4 exposure after recruitment. |",
      "| Catch | Quarterly catch-at-age in millions of fish, with separate age-0, ages 1-2, and age-3+ error groups. | Sum quarterly catches to annual numbers and retain the same three error groups. | The joint quarterly catch likelihood is not represented. |",
      "| Index | Five fleet series sampled in Q1 or Q3, with fleet-age q and SD sharing. | Keep direct age-specific observations and use q blocks that follow the reported sharing. | The grouped age-2-3+ IBTS observation is retained in the database but omitted because tinyAM has no grouped-age index likelihood. |",
      "| Weights and maturity | Stock weights are quarter-specific; maturity is 0, 0.2, 1, 1 at ages 0-3+. | Use Q1 stock weights for ages 1-3, Q3 weight for age 0, and the accepted maturity vector. | This matches the annual state boundary; maturity at age 0 is zero. |",
      "| SSB | Quarterly SSB is calculated at the start of each quarter. | Compare Q1 SSB using the same start-of-year N, weight, and maturity. | The source and tinyAM SSB definitions then match for the selected common years. |",
      "",
      "Survey timing uses the midpoint of each quarter: 0.125 for Q1 and 0.625 for Q3. Annual-equivalent F and M sum the quarterly rate over the time each age group is present; age-0 fish enter in Q3. The assessment's Q1 age-0 N precedes current-year recruitment, so that age-specific N row is omitted from the like-for-like N comparison. Recruitment at age 0 is compared separately. The grouped IBTS age-2-3+ index and the reported 2025 IBTS other-countries Q3 observations are unresolved in the cached fitted object; neither is fabricated or inferred."
    )
  )
}
