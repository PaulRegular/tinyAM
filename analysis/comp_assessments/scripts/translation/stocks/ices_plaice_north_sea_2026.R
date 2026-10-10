## Observations ----
years <- 1957:2025
ages <- 1:10
surveys <- c("BTS-Isis", "BTS-IBTS Q3", "SNS1", "SNS2", "IBTS Q1")
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)
obs$index$sd_block <- factor(obs$index$survey)
obs$index$q_key <- interaction(obs$index$survey, obs$index$age,
                               drop = TRUE, sep = ".")


dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages, N_settings = list(process = "iid",
    init = "exp"), F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6), M_settings = list(process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption), catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, fill_missing = FALSE))
start_par <- tinyAM::make_par(dat)


source_N <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages, 1000)
source_F <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age", years, ages)
start_par$log_r0 <- log(source_N[1L, 1L])
start_par$log_r <- log(source_N[-1L, 1L])
start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
start_par$log_f <- log(source_F)

index <- obs$index
n_at_age <- source_N[cbind(match(index$year, years),
                           match(index$age, ages))]
q_start <- tapply(index$obs / n_at_age, as.character(index$q_key),
                  stats::median, na.rm = TRUE)
if (any(!is.finite(q_start)) || any(q_start <= 0)) {
  stop("Could not obtain positive catchability starting values.")
}
names(q_start) <- paste0("q_key", names(q_start))
start_par$log_q[names(q_start)] <- log(q_start)


## Background and comparisons ----

age_plus_group <- 10

comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, F = 1, M = 1, F_bar = 1)

background <- c("### North Sea plaice: accepted 2026 SAM assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Years | Catch, biological inputs, and fitted N/F estimates cover 1957–2025; ICES also reports 2026 recruitment and SSB advice forecasts. | Fit 1957–2025. | 2026 is an advice forecast, not a full historical model year. |",
    "| Ages | Ages 1–10+, with age 10 as the plus group and recruitment at age 1. | Ages 1–10, with age 10 as the plus group. | The accepted age structure is retained. |",
    "| N | SAM estimates age-structured abundance and uses a separate process-variance key for age 1. | Exponential initial abundance with IID N process residuals; accepted N-at-age initializes the fit only. | tinyAM does not reproduce SAM's process-variance sharing; N is freely re-estimated. |",
    "| F | SAM estimates fishing mortality at age, with AR(1) correlation across ages and Fbar ages 2–6. | Random-walk F process with age-specific states and mean ages 2–6; accepted F-at-age initializes the fit only. | The correlation and parameter sharing in SAM are not reproduced exactly; F is freely re-estimated. |",
    "| M | Fixed, time-invariant age-specific values derived at the 2022 benchmark. | Supply the accepted M-at-age as fixed mortality. | The published age-specific values are retained. |",
    "| Catch | Aggregate landings and discards at age, including 50% of mature Division 7.d Q1 catch. | Use the published total catch-at-age series with one catch observation SD. | tinyAM aggregates the fishery and simplifies the catch-error structure. |",
    "| Index | Five age-specific series: BTS-Isis, BTS-IBTS Q3, SNS1, SNS2, and IBTS Q1. SAM uses fleet-specific q sharing and mixed observation-correlation structures. | Keep the five surveys separate, estimate q by survey and age, and estimate one observation SD per survey. | The report does not fully identify the printed correlation rows by survey, and the q constraints are simplified. |",
    "| Weights and maturity | Annual stock weight-at-age; maturity is constant at 0, 0.5, 0.5, then 1.0 from age 4. | Use the published stock weights and maturity ogive. | No biological surfaces are substituted. |",
    "| SSB | The accepted model reports SSB with 95% confidence intervals. | Calculate SSB using tinyAM's spawning-time convention. | A common definition is used for comparison; the ICES interval is shown only for the accepted series. |",
    "", "Survey index units are not stated in the report, so the published values remain on their native scales. Survey timing uses quarter midpoints: 0.75 for Q3 and 0.125 for Q1. These are approximations; exact SAM timing fractions were not recovered.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
