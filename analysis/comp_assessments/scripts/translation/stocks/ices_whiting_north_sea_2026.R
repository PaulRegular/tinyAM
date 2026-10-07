translate_stock <- function(source) {
  years <- 1978:2026
  ages <- 0:8
  assessment_id <- source$assessment$assessment_id

  inputs <- source$inputs
  m_rows <- inputs[inputs$assessment_id == assessment_id &
                     inputs$type == "M" &
                     inputs$measure == "natural_mortality_at_age", , drop = FALSE]
  last_m_year <- max(as.integer(m_rows$year), na.rm = TRUE)
  extra_years <- setdiff(years, unique(as.integer(m_rows$year)))
  if (length(extra_years)) {
    m_last <- m_rows[as.integer(m_rows$year) == last_m_year, , drop = FALSE]
    extra_m <- m_last[rep(seq_len(nrow(m_last)), each = length(extra_years)), , drop = FALSE]
    extra_m$year <- rep(extra_years, times = nrow(m_last))
    extra_m$source_type <- "translation_assumption"
    extra_m$transformation <- "Held the final numerical M input through the remaining model years as a fixed tinyAM assumption."
    extra_m$notes <- paste(
      "Translation-only baseline extension. The accepted SAM run estimates",
      "later M values with its mortality GMRF."
    )
    inputs <- rbind(inputs, extra_m)
  }

  obs <- database_to_tam_obs(
    assessment_id,
    inputs,
    years = years,
    ages = ages,
    assumptions = source$assumptions
  )

  obs$index$q_key <- interaction(obs$index$survey, obs$index$q_block,
                                 drop = TRUE, lex.order = TRUE)

  source_surface <- function(type, measure, multiplier = 1) {
    rows <- source$outputs[source$outputs$type == type &
                             source$outputs$measure == measure &
                             !is.na(source$outputs$age), , drop = FALSE]
    value <- matrix(NA_real_, length(years), length(ages),
                    dimnames = list(year = as.character(years),
                                    age = as.character(ages)))
    index <- cbind(match(as.character(rows$year), as.character(years)),
                   match(as.character(rows$age), as.character(ages)))
    value[index] <- as.numeric(rows$value) * multiplier
    if (any(!is.finite(value))) stop("The accepted source surface is incomplete.")
    value
  }

  source_N <- source_surface("population", "numbers_at_age", 1000)
  source_F <- source_surface("mortality", "fishing_mortality_at_age")
  source_M <- source_surface("mortality", "natural_mortality_at_age")
  q_design <- stats::model.matrix(~ 0 + q_key, data = obs$index)
  survey_year <- match(as.character(obs$index$year), rownames(source_N))
  survey_age <- match(as.character(obs$index$age), colnames(source_N))
  survey_rows <- cbind(survey_year, survey_age)
  n_at_survey <- source_N[survey_rows] * exp(
    -(source_F[survey_rows] + source_M[survey_rows]) * obs$index$samp_time
  )
  q_start <- vapply(seq_len(ncol(q_design)), function(j) {
    selected <- q_design[, j] == 1
    stats::median(log(obs$index$obs[selected] / n_at_survey[selected]))
  }, numeric(1))
  names(q_start) <- colnames(q_design)

  start_par <- list(
    log_r0 = log(source_N[1, "0"]),
    log_r = setNames(log(source_N[-1, "0"]), as.character(years[-1])),
    log_n = log(source_N[-1, as.character(ages[-1]), drop = FALSE]),
    log_f = log(source_F),
    log_q = q_start
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 8,
    obs = obs,
    start_par = start_par,
    comparison_scales = c(
      N = 1e-3, recruitment = 1e-3, ssb = 1e-3, biomass = 1e-3,
      abundance = 1e-3, biomass_at_age = 1e-3, F_bar = 1
    ),
    settings = list(
      N_settings = list(process = "rw", init = "exp"),
      F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:5),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 0 + q_key,
                            sd_form = ~ 0 + survey,
                            sd_supplied = ~ relative_sd,
                            fill_missing = FALSE)
    ),
    background = c(
      "### North Sea whiting: accepted 2026 SAM assessment",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | Catch covers 1978–2025; Q1 covers 1983–2026; Q3 covers 1991–2025. | Fit 1978–2026. | Retains the accepted survey terminal year. |",
      "| Ages | Ages 0–8+, with recruitment at age 0. | Ages 0–8, with age 8 as the plus group. | The source age range is retained. |",
      "| N | Random-walk recruitment and abundance with age-specific variance sharing. | Random-walk recruitment, exponential initial abundance, and a cohort random walk with one SD. | tinyAM does not use SAM's age-specific variance sharing or initial-state integration. |",
      "| F | Random-walk F with AR(1) correlation between age increments; Fbar is ages 2–5. | Random-walk F with independent age increments; Fbar is ages 2–5. | tinyAM does not model age correlation in F increments. |",
      "| M | WGSAM M observations through 2022; SAM estimates later M with a GMRF. | The 2022 numerical M surface is carried forward through 2026 as a fixed assumption. | tinyAM's estimated M process did not produce a usable fit with these data; later M therefore remains a documented source of mismatch. |",
      "| Catch | Total catch at age in thousands of fish, with independent lognormal age residuals and shared variance groups. | Source catches converted to individual fish and fitted with one lognormal SD. | tinyAM does not share catch SDs by the source age groups. |",
      "| Indices | IBTS Q1 ages 1–6+ and Q3 ages 0–6+; Q1 timing is 0.125 and Q3 timing is 0.625. q varies by survey and age. | Source indices, timing, survey-age q, and relative log-scale SD factors are retained. | tinyAM uses independent age residuals instead of the source AR(1) age errors. |",
      "| Weights and maturity | Annual stock and catch weights and a smoothed maturity surface; spawning fractions are zero. | Native biological surfaces and timing are retained. | No biological inputs are substituted. |",
      "",
      "Accepted SAM N and F surfaces initialize the fit only; tinyAM re-estimates them. The native number scale is thousands of fish, so catch and population counts are translated to individual fish; biomass comparisons are scaled back to tonnes. The published 2026 table reports SSB only, while native state surfaces are available through 2026."
    )
  )
}
