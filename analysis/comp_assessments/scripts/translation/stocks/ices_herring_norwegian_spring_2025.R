translate_stock <- function(source) {
  years <- 1988:2024
  ages <- 2:12
  inputs <- source$inputs
  maturity <- inputs$type == "maturity" & inputs$measure == "maturity_at_age"
  inputs$year[maturity] <- as.character(
    as.integer(inputs$year[maturity]) + as.integer(inputs$age[maturity])
  )
  inputs$year_basis[maturity] <- "calendar_year"
  natural_mortality <- inputs$type == "M" &
    inputs$measure == "natural_mortality_at_age"
  inputs <- inputs[!natural_mortality | as.integer(inputs$age) %in% ages, ]

  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    inputs,
    years = years,
    ages = ages,
    surveys = c("NASF_1988_2008", "NASF_2015_onward", "IESNS_Barents",
                "IESNS_Norwegian_Sea", "BESS"),
    assumptions = source$assumptions
  )

  settings <- list(
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 5:12),
    M_settings = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ 0 + q_key,
                          sd_form = ~ 0 + survey,
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
      stop("Source ", measure, " does not cover the tinyAM start surface.")
    }
    surface
  }
  source_N <- source_surface("population", "numbers_at_age", 1e6)
  source_F <- source_surface("mortality", "fishing_mortality_at_age")
  start_par$log_r0 <- log(source_N[1L, 1L])
  start_par$log_n0 <- log(source_N[1L, -1L])
  start_par$log_r <- log(source_N[-1L, 1L])
  start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
  start_par$log_f <- log(source_F)
  index <- obs$index
  n_at_age <- source_N[cbind(match(index$year, years), match(index$age, ages))]
  q_start <- tapply(index$obs / n_at_age, as.character(index$q_key),
                    stats::median, na.rm = TRUE)
  start_par$log_q[paste0("q_key", names(q_start))] <- log(q_start)

  list(
    years = years,
    ages = ages,
    age_plus_group = 12,
    obs = obs,
    start_par = start_par,
    settings = settings,
    comparison_scales = c(
      N = 1e-6, recruitment = 1e-6, ssb = 1e-6, F_bar = 1
    ),
    background = c(
      "### Norwegian spring-spawning herring: accepted 2025 SAM assessment",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | Catch-at-age is fitted through 2024; reported N and F estimates extend through 2025. | Fit 1988–2024. | This uses the accepted fitted period with catch data and compares the common historical years. |",
      "| Ages | Ages 2–12+, with recruitment at age 2. Source catch and biological tables retain detail through 15+. | Ages 2–12, with age 12 as the plus group. | Catch and index numbers older than 12 are summed; age-12 weight and maturity represent the plus group. |",
      "| N | The accepted assessment uses SAM's age-structured population process; the full configuration was not recovered. | IID abundance process with freely estimated initial abundance. | The exact SAM process and initial-state treatment are unavailable. |",
      "| F | SAM estimates age-specific fishing mortality; reported Fbar covers ages 5–12+. | Random-walk F with age-specific states and Fbar ages 5–12. | The native SAM transition and age-sharing settings were not recovered. |",
      "| M | Standard M is 0.9 per year at age 2 and 0.15 at ages 3+; the report notes additional annual deviations. | Supply the standard age pattern as fixed M. | The annual stock-annex deviations were not recovered. |",
      "| Catch | One aggregate catch-at-age series in thousand fish, with age-specific sampling errors. | Convert to fish and fit a lognormal likelihood with one shared observation SD. | tinyAM estimates its own error rather than applying SAM's external relative errors. |",
      "| Index | NASF, IESNS Barents, IESNS Norwegian Sea, and BESS have age-specific abundance indices. IESNS Barents is fitted at age 2 only, and its 2008 value is missing. RFID values are only available as a figure. | Retain the four numerical series; use separate q by survey and age and separate observation SDs by survey. | RFID observations are omitted because numerical values could not be recovered; SAM's q-sharing and external error weights are not reproduced. |",
      "| Weights and maturity | Annual stock and catch weights; maturity ogives vary by birth cohort. | Use annual stock weights and map maturity cohort to year using cohort = year - age. | The age-12 biological values stand in for ages 12+ because single-age abundance above age 12 is unavailable. |",
      "| SSB | Reported SSB is in thousand tonnes and Fbar is ages 5–12+. | Calculate SSB from translated weight, maturity and M; compare on the common fitted years. | The exact SAM spawning-time and annual M adjustments are not fully available. |",
      "",
      "The source reports catch, survey, and maturity inputs in different age groupings. tinyAM sums catch and index numbers above age 12 into its plus group. Maturity rows are reported by birth cohort; for a year-age observation the matching cohort is year minus age. Relative standard errors are retained in the database but are not applied in this fit. The available SAM N and F estimates and index-to-N ratios supply starting values only; tinyAM estimates them freely."
    )
  )
}
