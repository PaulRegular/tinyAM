## Observations ----
years <- 1967:2025
ages <- 3:10
survey <- "NS-IBTS Q3-Q4"
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = survey,
  sampling_times = c("NS-IBTS Q3-Q4" = 0.75),
  assumptions = source$assumptions
)

obs$catch$sd_block <- cut_ages(obs$catch$age, c(3, 4, 6, 10))
levels(obs$catch$sd_block) <- c("age3", "age4_5", "age6_plus")
obs$index$q_age <- factor(obs$index$age, levels = 3:8)


dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages, N_settings = list(process = "iid",
    init = "exp"), F_settings = list(process = "ar1", mu_form = ~factor(age), mean_ages = 4:7), M_settings = list(process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption), catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_age, sd_form = ~1, fill_missing = FALSE))
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

source_N <- source_surface("population", "numbers_at_age", 1000)
source_F <- matrix(NA_real_, length(years), length(ages),
                   dimnames = dimnames(source_N))
f_rows <- source$outputs[source$outputs$type == "mortality" &
                           source$outputs$measure == "fishing_mortality_at_age" &
                           source$outputs$year %in% years, , drop = FALSE]
f_index <- cbind(match(as.integer(f_rows$year), years),
                 match(as.integer(f_rows$age), ages))
source_F[f_index] <- as.numeric(f_rows$value)
source_F[, "10"] <- source_F[, "9"]
if (any(!is.finite(source_F)) || any(source_F <= 0)) {
  stop("The accepted F surface is incomplete or non-positive.")
}

start_par$log_r0 <- log(source_N[1L, 1L])
start_par$log_r <- log(source_N[-1L, 1L])
start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
start_par$log_f <- log(source_F)

index <- obs$index
n_at_age <- source_N[cbind(match(index$year, years),
                           match(index$age, ages))]
q_start <- tapply(index$obs / n_at_age, as.character(index$q_age),
                  stats::median, na.rm = TRUE)
q_age <- sub("^q_age", "", names(start_par$log_q))
if (any(!is.finite(q_start)) || any(q_start <= 0) ||
    !setequal(q_age, names(q_start))) {
  stop("Could not obtain positive source-based q starting values.")
}
start_par$log_q[] <- log(q_start[q_age])

comparison_outputs <- source$outputs
f_plus <- comparison_outputs[
  comparison_outputs$type == "mortality" &
    comparison_outputs$measure == "fishing_mortality_at_age" &
    !is.na(comparison_outputs$age_group) &
    comparison_outputs$age_group == "9+", , drop = FALSE
]
f_plus$age <- 10L
f_plus$age_group <- "10+"
f_plus$notes <- paste(f_plus$notes,
                      "The source model constrains ages 9 and 10+ to one F state; duplicated here only to align common comparison ages.")
comparison_outputs <- rbind(comparison_outputs, f_plus)


## Background and comparisons ----

age_plus_group <- 10


comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, F = 1, M = 1, F_bar = 1)

background <- c("### North Sea saithe: accepted 2026 SAM assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Years | The fitted assessment covers 1967–2025. A separate 2026 short-term forecast is reported. | Fit 1967–2025. | The 2026 forecast is not a full historical assessment year. |",
    "| Ages | Ages 3–10+, with recruitment at age 3. SAM shares the F state for ages 9 and 10+. | Ages 3–10, with age 10 as the plus group. | The comparison repeats the single accepted 9+ F state across ages 9 and 10+ as specified by SAM. |",
    "| N | SAM estimates abundance at ages 3–10+, with separate process variance for recruitment and shared variance for older ages. | Random-walk recruitment with its own SD, exponential initial abundance and shared IID survival SD; accepted N-at-age initializes the fit only. | Recruitment-versus-survival variance sharing is retained; initial-state integration differs. |",
    "| F | SAM estimates age-specific F with AR(1) age correlation, four innovation-SD groups and shared states at ages 9–10+. Fbar is ages 4–7. | AR(1) F process around estimated age-specific means, one process SD and Fbar ages 4–7; accepted F-at-age initializes the fit only. | The shared 9+ state, age-specific innovation SDs and exact RW covariance are not reproduced. |",
    "| M | Fixed, age-specific natural mortality derived from mean stock weights under Lorenzen's relationship. | Use the published fixed M-at-age values. | The accepted numerical values are retained. |",
    "| Catch | Combined catch numbers-at-age and catch weight-at-age, ages 3–10+. | Use both reported surfaces with three catch-SD age groups. | The fleet is combined and SAM's cross-age observation correlation is not reproduced. |",
    "| Index | Q3–Q4 research-vessel index at ages 3–8, plus annual commercial CPUE tuned to exploitable biomass. | Fit the age-specific survey index with q estimated by age. | tinyAM does not represent the aggregate biomass-targeted CPUE observation; it is retained in the database but omitted from the fit. |",
    "| Weights and maturity | Annual stock weights, annual maturity-at-age, and annual catch weights are supplied to SAM as known. | Use the published annual surfaces without replacing them. | Stock weights are kept distinct from catch weights. |",
    "| SSB | SAM reports SSB with 95% confidence intervals; its exact spawning-time detail is unresolved. | Compare start-year SSB calculated from accepted N and tinyAM N using the same translated weights and maturity. | This common-definition calculation is separate from native reported SAM SSB and its intervals. |",
    "", "The age-specific survey timing is approximated at 0.75, the midpoint of Q3–Q4. The report does not give the exact 2026 sampling fraction. The assessment configuration reproduced in the 2026 report is timestamped 2024.")

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
    F_settings = list(process = "ar1", mu_form = ~factor(age), mean_ages = 4:7),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_age, sd_form = ~1, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
