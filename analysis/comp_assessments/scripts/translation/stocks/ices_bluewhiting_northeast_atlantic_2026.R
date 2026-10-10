## Observations ----
years <- 1981:2026
ages <- 1:10
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)

obs$catch$sd_group <- cut_ages(obs$catch$age, c(1, 2, 3, 9, 10))
levels(obs$catch$sd_group) <- c("1", "2", "3-8", "9-10", "9-10")
q_group <- as.character(cut_ages(obs$index$age, c(1:5, 8)))
obs$index$sd_group <- cut_ages(obs$index$age, c(1:4, 7, 8))
levels(obs$index$sd_group) <- c("1", "2", "3", "4-6", "7-8", "7-8")
obs$index$q_key <- interaction(obs$index$survey, q_group,
                               drop = TRUE, sep = ".")

source_surface <- function(type, measure, multiplier = 1) {
  rows <- source$outputs[source$outputs$type == type &
                           source$outputs$measure == measure &
                           !is.na(source$outputs$age), , drop = FALSE]
  value <- matrix(NA_real_, length(years), length(ages),
                  dimnames = list(year = as.character(years),
                                  age = as.character(ages)))
  row <- match(as.integer(rows$year), years)
  age <- match(as.integer(rows$age), ages)
  value[cbind(row, age)] <- as.numeric(rows$value) * multiplier
  if (any(!is.finite(value))) stop("The accepted source surface is incomplete.")
  value
}

source_N <- source_surface("population", "numbers_at_age", multiplier = 1000)
source_F <- source_surface("mortality", "fishing_mortality_at_age")
q_rows <- source$outputs[source$outputs$measure == "q" &
                           source$outputs$survey == "IBWSS", , drop = FALSE]
q_age_groups <- c("1", "2", "3", "4", "5-8")
q_start <- vapply(q_age_groups, function(group) {
  if (group == "5-8") {
    rows <- q_rows[as.character(q_rows$age_group) == group, , drop = FALSE]
  } else {
    rows <- q_rows[as.character(q_rows$age) == group, , drop = FALSE]
  }
  unique(as.numeric(rows$value)) * 1000
}, numeric(1))
names(q_start) <- paste0("q_keyIBWSS.", q_age_groups)
process_sd <- source$outputs[source$outputs$measure == "process_sd", , drop = FALSE]
catch_sd <- source$outputs[source$outputs$type == "catch" &
                             source$outputs$measure == "observation_sd", , drop = FALSE]
index_sd <- source$outputs[source$outputs$type == "index" &
                             source$outputs$measure == "observation_sd", , drop = FALSE]
sd_starts <- function(rows, groups) {
  values <- as.numeric(rows$value[match(groups, as.character(rows$age_group))])
  setNames(log(values), paste0("sd_group", groups))
}
start_par <- list(
  log_r0 = log(source_N[1, 1]),
  log_n0 = setNames(log(source_N[1, -1]), as.character(ages[-1])),
  log_r = setNames(log(source_N[-1, 1]), as.character(years[-1])),
  log_f = log(source_F),
  log_sd_r = log(as.numeric(process_sd$value[process_sd$type == "population" &
                                                  process_sd$age_group == "1"])),
  log_sd_f = log(as.numeric(process_sd$value[process_sd$type == "mortality"])),
  log_sd_catch = sd_starts(catch_sd, c("1", "2", "3-8", "9-10")),
  log_sd_index = sd_starts(index_sd, c("1", "2", "3", "4-6", "7-8")),
  log_q = log(q_start)
)


## Background and comparisons ----

age_plus_group <- 10

comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, biomass = 0.001, biomass_at_age = 0.001,
    F_bar = 1)

background <- c("### Northeast Atlantic blue whiting: accepted 2026 SAM assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Years | Catch covers 1981–2026; IBWSS covers 2004–2026 except 2010 and 2020. | Fit 1981–2026. | All native model years are retained. |",
    "| Ages | Ages 1–10+, with recruitment at age 1 and age 10 as the plus group. | Ages 1–10, with age 10 as the plus group. | The source age range is retained. |",
    "| N | Age-1 recruitment follows a random walk; age-specific abundance deviations have separate age-1 and shared ages 2–10 variances. | Random-walk recruitment and free first-year abundance by age, followed by deterministic cohort survival. | The fuller N-process fit did not converge; initial abundance is not integrated as in SAM. |",
    "| F | Annual random-walk F has age-correlated increments; age 10+ shares the age-9 state. Fbar is ages 3–7. | Random-walk F with separate age increments; Fbar is set to ages 3–7. | tinyAM does not reproduce the source correlation between age increments. |",
    "| M | Fixed natural mortality is 0.2 per year at all ages. | The same fixed M surface is supplied to tinyAM. | This assumption is retained. |",
    "| Catch | One commercial catch-at-age series in thousands of fish; lognormal errors are age-correlated and use four age groups. | Native catch converted to fish; independent lognormal errors share SDs over the same four age groups. | tinyAM does not model age-correlated observation errors. |",
    "| Index | IBWSS ages 1–8, in millions of fish, sampled at 0.245 of the year. Catchability uses age groups 1, 2, 3, 4, and 5–8; observation SD uses 1, 2, 3, 4–6, and 7–8. | Values are converted to fish, with the source catchability and observation-SD groupings preserved separately. | tinyAM omits age correlation in observation errors. |",
    "| Weights and maturity | Annual stock and catch weights; time-invariant maturity-at-age. | Native biological surfaces are retained. | No biological inputs are substituted. |",
    "", "Initial N, F, and q values are taken from the accepted SAM fit only as starting values; tinyAM re-estimates them. tinyAM converts the native catch and survey count units to fish internally. Comparison outputs are scaled back to the source units: thousands of fish and tonnes. The 2026 advice replaces the fitted recruitment estimate with a forecast assumption; the comparison reference keeps the native fitted estimate.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "free"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 3:7),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_group, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_group, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
