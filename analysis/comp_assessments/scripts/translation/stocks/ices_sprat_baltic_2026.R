translate_stock <- function(source) {
  years <- 1974:2025
  ages <- 1:8
  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    source$inputs,
    years = years,
    ages = ages,
    surveys = c(
      "BIAS_October_SD22_29_32_recent",
      "BIAS_October_SD22_29_early",
      "BASS_May_SD24_26_28",
      "BIAS_October_age0_shifted_to_age1"
    ),
    assumptions = source$assumptions
  )

  obs$index$sd_block <- factor(obs$index$survey)
  q_age <- as.character(cut_ages(obs$index$age, c(1:6, 8)))
  obs$index$q_key <- interaction(obs$index$survey, q_age,
                                 drop = TRUE, sep = ".")

  settings <- list(
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 3:5),
    M_settings = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ 0 + q_key,
                          sd_form = ~ 0 + sd_block,
                          fill_missing = FALSE)
  )

  dat <- do.call(tinyAM::prepare_tam, c(
    list(data = obs, years = years, ages = ages), settings
  ))
  start_par <- tinyAM::make_par(dat)

  source_surface <- function(type, measure, multiplier = 1) {
    rows <- source$outputs[source$outputs$type == type &
                             source$outputs$measure == measure &
                             source$outputs$year %in% years &
                             source$outputs$age %in% ages, , drop = FALSE]
    surface <- matrix(NA_real_, length(years), length(ages),
                      dimnames = list(as.character(years), as.character(ages)))
    index <- cbind(match(as.integer(rows$year), years),
                   match(as.integer(rows$age), ages))
    surface[index] <- as.numeric(rows$value) * multiplier
    if (any(!is.finite(surface)) || any(surface <= 0)) {
      stop("The accepted ", measure, " surface is incomplete or non-positive.")
    }
    surface
  }

  source_N <- source_surface("population", "numbers_at_age", 1e6)
  source_F <- source_surface("mortality", "fishing_mortality_at_age")
  start_par$log_r0 <- log(source_N[1L, 1L])
  start_par$log_r <- log(source_N[-1L, 1L])
  start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
  start_par$log_f <- log(source_F)

  index <- obs$index
  n_at_age <- source_N[cbind(match(index$year, years), match(index$age, ages))]
  q_start <- tapply(index$obs / n_at_age, as.character(index$q_key),
                    stats::median, na.rm = TRUE)
  names(q_start) <- paste0("q_key", names(q_start))
  start_par$log_q[names(q_start)] <- log(q_start)

  list(
    years = years,
    ages = ages,
    age_plus_group = 8,
    obs = obs,
    start_par = start_par,
    settings = settings,
    comparison_scales = c(
      N = 1e-6, recruitment = 1e-6, ssb = 1e-3, F = 1, M = 1,
      F_bar = 1
    ),
    background = c(
      "### Baltic sprat: accepted 2026 SAM assessment",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | Catch and full historical inputs cover 1974–2025; the report also gives an intermediate-year age-1 estimate for 2026. | Fit 1974–2025. | The 2026 age-0 acoustic observation is shifted to age 1 in 2026 and is not a full catch-data year. |",
      "| Ages | Ages 1–8+, with age 8 as the plus group. | Ages 1–8, with age 8 as the plus group. | The reported age range is retained. |",
      "| N | Recruitment follows a random walk; process variance differs between age 1 and ages 2–8. | Random-walk recruitment, shared IID cohort residuals, and exponential initial abundance. | This approximates the separate recruitment and older-age process variation; first-year older ages follow deterministic survivorship from recruitment. |",
      "| F | SAM uses a temporal F process with AR(1) correlation across age states; Fbar is ages 3–5. | Random-walk F with independent age increments; Fbar is ages 3–5. | tinyAM does not include the source's cross-age correlation or shared oldest-age F state. |",
      "| M | Annual age-specific natural mortality varies with cod predation. | The reported M-at-age surface is supplied as fixed mortality. | The numerical source surface is retained; tinyAM does not reproduce how SMS generated M. |",
      "| Catch | One aggregate catch-at-age series in thousands of fish. | Source values are converted to fish and fitted with one catch-error SD. | The catch series is retained; the source's observation error covariance is simplified. |",
      "| Index | Four tuning fleets: three age-specific acoustic series and one age-0 series shifted to age 1 in the following year. SAM uses fleet-specific q by age and a year-class effect at age 1. | Retain four survey identities, share q over ages 6–8, and use one observation SD per survey. | The shifted age-0 series is aligned to the report's age-1 year label; q's year-class effect and source observation covariance are not represented. |",
      "| Weights and maturity | Catch and stock weights are equal and annual. Maturity is constant at 0.17, 0.93, then 1.0 for ages 3–8. | Use the reported weight and maturity inputs. | No biological surface is substituted. |",
      "| SSB | SAM reports SSB and Fbar estimates with lower and upper bounds. | tinyAM calculates SSB from the supplied biology and its spawning-time convention. | The source's 40% pre-spawning fractions for F and M are documented, but tinyAM uses its own timing convention. |",
      "",
      "SAM estimates for N and F initialize the fit only; tinyAM re-estimates them. The 2026 intermediate-year recruitment represents the 2025 year class and is excluded from this fit and common-period comparison. Survey timing is approximated as 0.8 for October and 0.375 for May."
    )
  )
}
