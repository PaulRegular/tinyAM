## Observations ----
years <- 1967:2025
ages <- 3:8
source_ages <- 3:10
survey <- "NS-IBTS Q3-Q4"
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = source_ages,
  surveys = survey,
  sampling_times = c("NS-IBTS Q3-Q4" = 0.75),
  assumptions = source$assumptions
)

obs$catch$sd_block <- cut_ages(obs$catch$age, c(3, 4, 6, 10))
levels(obs$catch$sd_block) <- c("age3", "age4_5", "age6_plus")
obs$index$q_age <- cut_ages(obs$index$age, c(3, 4, 5, 6, 8))

dat <- tinyAM::prepare_tam(
  data = obs, years = years, ages = ages, N_settings = list(
    process = "iid",
    init = "exp"
  ), F_settings = list(process = "ar1", mu_form = ~ factor(age), mean_ages = 4:7), M_settings = list(
    process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption
  ), catch_settings = list(sd_form = ~ 0 + sd_block),
  index_settings = list(q_form = ~ 0 + q_age, sd_form = ~1)
)
start_par <- tinyAM::make_par(dat)


source_N <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, source_ages, 1000)
source_F <- matrix(NA_real_, length(years), length(source_ages),
  dimnames = dimnames(source_N)
)
f_rows <- source$outputs[source$outputs$type == "mortality" &
  source$outputs$measure == "fishing_mortality_at_age" &
  source$outputs$year %in% years, , drop = FALSE]
f_index <- cbind(
  match(as.integer(f_rows$year), years),
  match(as.integer(f_rows$age), source_ages)
)
source_F[f_index] <- as.numeric(f_rows$value)
source_F[, "10"] <- source_F[, "9"]
if (any(!is.finite(source_F)) || any(source_F <= 0)) {
  stop("The accepted F surface is incomplete or non-positive.")
}

plus <- source_ages >= max(ages)
source_F[, as.character(max(ages))] <-
  rowSums(source_F[, plus] * source_N[, plus]) / rowSums(source_N[, plus])
source_N[, as.character(max(ages))] <- rowSums(source_N[, plus])
source_F <- source_F[, as.character(ages), drop = FALSE]
source_N <- source_N[, as.character(ages), drop = FALSE]

start_par$log_r0 <- log(source_N[1L, 1L])
start_par$log_r <- log(source_N[-1L, 1L])
start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
start_par$log_f <- log(source_F)

index <- obs$index
n_at_age <- source_N[cbind(
  match(index$year, years),
  match(index$age, ages)
)]
q_start <- tapply(index$obs / n_at_age, as.character(index$q_age),
  stats::median,
  na.rm = TRUE
)
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
    comparison_outputs$age_group == "9+", ,
  drop = FALSE
]
f_plus$age <- 10L
f_plus$age_group <- "10+"
f_plus$notes <- paste(
  f_plus$notes,
  "The source model constrains ages 9 and 10+ to one F state; duplicated here only to align common comparison ages."
)
comparison_outputs <- rbind(comparison_outputs, f_plus)
shared_f <- comparison_outputs$measure == "fishing_mortality_at_age" &
  comparison_outputs$age_group %in% "9+"
comparison_outputs$age_group[shared_f] <- NA_character_
comparison_outputs$notes[shared_f] <- paste(
  comparison_outputs$notes[shared_f],
  "The shared F state is displayed at age 9 and repeated at age 10+; it is not a pooled age-group rate."
)


## Background and comparisons ----

age_plus_group <- 8


comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, F = 1, M = 1, F_bar = 1)

background <- c(
  "### North Sea saithe: accepted 2026 SAM assessment", "", print_sources(source$assessment),
  "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
  "| Years | The fitted assessment covers 1967–2025. A separate 2026 short-term forecast is reported. | Fit 1967–2025. | The 2026 forecast is not a full historical assessment year. |",
  "| Ages | Ages 3–10+, with recruitment at age 3. SAM shares the F state for ages 9 and 10+. | Ages 3–8+, with age 8 as the plus group. | Catch and accepted abundance at ages 8–10+ are summed; accepted F is abundance-weighted over those ages. |",
  "| N | SAM estimates abundance at ages 3–10+, with separate process variance for recruitment and shared variance for older ages. | Random-walk recruitment with its own SD, exponential initial abundance and shared IID survival SD; accepted N-at-age initializes the fit only. | Recruitment-versus-survival variance sharing is retained; initial-state integration differs. |",
  "| F | SAM estimates age-specific F with AR(1) age correlation, four innovation-SD groups and shared states at ages 9–10+. Fbar is ages 4–7. | Age-correlated F random walks, one innovation SD and Fbar ages 4–7; accepted F-at-age initializes the fit only. | The age-8 plus group replaces separate older-age F states; age-specific innovation SDs are not reproduced. |",
  "| M | Fixed, age-specific natural mortality derived from mean stock weights under Lorenzen's relationship. | Use published M at modeled ages; age-8 M applies to the whole 8+ group. | Older fish share age-8 mortality in the reduced age structure. |",
  "| Catch | Combined catch numbers-at-age and catch weight-at-age, ages 3–10+. | Use both reported surfaces with three catch-SD age groups. | The fleet is combined and SAM's cross-age observation correlation is not reproduced. |",
  "| Index | Q3–Q4 research-vessel index at ages 3–8, plus annual commercial CPUE tuned to exploitable biomass. | Use `mono(q_age)` for non-decreasing survey catchability across age groups 3, 4, 5 and 6–8, with one shared q for ages 6–8. The age-8 index is linked to the modeled 8+ group. | The monotone constraint limits the trade-off between declining catchability and inflated abundance at older ages. The terminal survey age represents only age 8 in SAM. tinyAM does not represent the aggregate biomass-targeted CPUE observation; it is retained in the database but omitted from the fit. |",
  "| Weights and maturity | Annual stock weights, annual maturity-at-age, and annual catch weights are supplied to SAM as known. | Retain biological ages 3–10 for within-plus-group weight and maturity calculations. | Stock weights are kept distinct from catch weights. |",
  "| SSB | SAM reports SSB with 95% confidence intervals; its exact spawning-time detail is unresolved. | Compare start-year SSB calculated from accepted N and tinyAM N using the same translated weights and maturity. | This common-definition calculation is separate from native reported SAM SSB and its intervals. |",
  "", "The age-8 plus group and non-decreasing survey catchability improved agreement with the Accepted assessment and reduced excess biomass accumulation at older ages. These are tinyAM translation choices rather than assumptions imposed by SAM.",
  "", "The age-specific survey timing is approximated at 0.75, the midpoint of Q3–Q4. The report does not give the exact 2026 sampling fraction. The assessment configuration reproduced in the 2026 report is timestamped 2024."
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
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "cor_rw", mean_ages = 4:7),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_block),
    index_settings = list(q_form = ~ mono(q_age), sd_form = ~1),
    silent = silent,
    start_par = start_par
  )
}
