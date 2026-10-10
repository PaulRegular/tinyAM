## Observations ----
years <- 1991:2025
ages <- 0:8
assessment_id <- source$assessment$assessment_id
inputs <- source$inputs

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = c("HERAS", "GERAS", "N20", "IBTS/BITSQ1"),
  assumptions = source$assumptions
)

catch_block <- cut_ages(obs$catch$age, c(0, 1, 2, 8))
levels(catch_block) <- c("sam_sd_5", "sam_sd_6", "sam_sd_0")
obs$catch$sd_block <- factor(as.character(catch_block))
geras <- obs$index$survey == "GERAS"
geras_q <- cut_ages(obs$index$age[geras], 1:3)
levels(geras_q) <- c("sam_q_2", "sam_q_2", "sam_q_7")
obs$index$q_key <- with(obs$index, ifelse(
  survey == "HERAS", ifelse(age == 2, "sam_q_0", "sam_q_1"),
  ifelse(survey == "N20", "sam_q_3", paste0("sam_q_", age + 1L))
))
obs$index$q_key[geras] <- as.character(geras_q)
obs$index$q_key <- factor(obs$index$q_key)
obs$index$sd_block <- factor(ifelse(
  obs$index$survey == "HERAS", "sam_sd_1",
  ifelse(obs$index$survey == "GERAS", "sam_sd_2",
         ifelse(obs$index$survey == "N20", "sam_sd_3", "sam_sd_4"))
))
obs$index$relative_sd[!is.finite(obs$index$relative_sd)] <- 1


dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages, N_settings = list(process = "rw",
    init = "exp"), F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:5), M_settings = list(process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption), catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, sd_supplied = ~relative_sd, fill_missing = FALSE))
start_par <- tinyAM::make_par(dat)

source_surface <- function(type, measure, multiplier = 1) {
  rows <- source$outputs[
    source$outputs$type == type & source$outputs$measure == measure &
      !is.na(source$outputs$year) & !is.na(source$outputs$age) &
      as.integer(source$outputs$year) %in% years &
      as.integer(source$outputs$age) %in% ages,
    , drop = FALSE
  ]
  surface <- matrix(
    NA_real_, length(years), length(ages),
    dimnames = list(as.character(years), as.character(ages))
  )
  index <- cbind(match(as.integer(rows$year), years),
                 match(as.integer(rows$age), ages))
  if (anyNA(index) || anyDuplicated(paste(index[, 1], index[, 2]))) {
    stop("Source ", measure, " rows do not map uniquely to the model grid.")
  }
  surface[index] <- as.numeric(rows$value) * multiplier
  if (any(!is.finite(surface)) || any(surface <= 0)) {
    stop("The accepted ", measure, " surface is incomplete or non-positive.")
  }
  surface
}

source_N <- source_surface("population", "numbers_at_age", 1000)
source_F <- source_surface("mortality", "fishing_mortality_at_age")
source_M <- source_surface("mortality", "natural_mortality_at_age")
start_par$log_r0 <- log(source_N[1L, "0"])
start_par$log_r <- setNames(log(source_N[-1L, "0"]),
                            as.character(years[-1L]))
start_par$log_n <- log(source_N[-1L, as.character(ages[-1L]), drop = FALSE])
start_par$log_f <- log(source_F)

index <- obs$index
q_design <- stats::model.matrix(~ 0 + q_key, data = index)
row_index <- cbind(match(as.character(index$year), rownames(source_N)),
                   match(as.character(index$age), colnames(source_N)))
n_at_survey <- source_N[row_index] *
  exp(-(source_F[row_index] + source_M[row_index]) * index$samp_time)
q_start <- vapply(seq_len(ncol(q_design)), function(j) {
  selected <- q_design[, j] == 1
  log(stats::median(index$obs[selected] / n_at_survey[selected],
                    na.rm = TRUE))
}, numeric(1))
names(q_start) <- colnames(q_design)
if (!setequal(names(q_start), names(start_par$log_q)) ||
    any(!is.finite(q_start))) {
  stop("Could not create finite source-based starting values for q.")
}
start_par$log_q[] <- q_start[names(start_par$log_q)]


## Background and comparisons ----

age_plus_group <- 8

comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, biomass = 0.001, abundance = 0.001,
    biomass_at_age = 0.001, F = 1, M = 1, F_bar = 1)

background <- c("### Western Baltic spring-spawning herring: accepted 2026 SAM assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Main simplification |", "|---|---|---|---|",
    "| Years and ages | 1991–2025; ages 0–8+, with recruitment at age 0. | Fit ages 0–8, with age 8 as the plus group. | The source age range is retained. |",
    "| N | SAM uses a segmented recruitment relationship, a separate recruitment variance, and shared survival variance for ages 1–8. | Random-walk recruitment with its own SD, a cohort random walk with shared SD, and exponential initial abundance. | The hockey-stick relationship, survival-error dynamics and initial-state integration differ; separate recruitment and survival SDs are available. |",
    "| F | SAM has age-specific F, AR(1) dependence across ages, and one shared state for ages 7–8; Fbar is ages 2–5. | Random-walk F with independent age processes. | tinyAM does not represent the source cross-age covariance or the shared oldest-age state. |",
    "| M | Fixed, time-invariant age-specific M derived from NSAS herring and profiled at the 2025 benchmark. | Use the same fixed M-at-age values in every year. | M is not estimated in either fit. |",
    "| Catch | Total catch at age in thousands of fish; lognormal errors have separate variances for ages 0, 1, and 2–8. | Convert to individual fish and retain those three observation-SD groups. | The source's other model details are simplified. |",
    "| Indices | HERAS ages 2–6+, GERAS ages 1–3, N20 age 0, and IBTS/BITS Q1 ages 3–5+; timings are 0.625, 0.8, 0.4, and 0.136365. | Retain all available observations, exact timing, q-sharing groups, and IBTS/BITS relative SD factors. | tinyAM uses independent age residuals rather than the source's AR(1) age errors. Five 1999 HERAS values are missing in the source. |",
    "| Weights and maturity | Annual stock and catch weights; constant maturity ogive; spawning fractions are 0.168 for F and 0.25 for M. | Use source weights and maturity, without the pre-spawning survival adjustment. | Biological inputs are retained; spawning timing differs. |",
    "| SSB | Based on source N, weights and maturity after survival to spawning. | Calculate start-year SSB; use accepted N with the same translated biology for common-definition comparisons. | The common comparison excludes spawning-time survival and is distinct from native SAM SSB. |",
    "", "Source N, F, and q values initialize the tinyAM fit but do not constrain it. The observation precision weight provided for IBTS/BITS Q1 is retained with its equivalent SD multiplier; absent source weights on the other indices use a neutral multiplier of one.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "rw", init = "exp"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, sd_supplied = ~relative_sd, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
