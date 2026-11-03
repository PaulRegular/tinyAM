source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1982:2025
ages <- 1:11
assessment_id <- source$assessment$assessment_id[[1L]]
surveys <- c("UK(E&W)-BTS-Q3", "UK(E&W)-YFS", "FR-YFS")

inputs <- source$inputs
index_measure <- inputs$type == "index" & inputs$measure == "numbers_at_age"
index_rows <- inputs[index_measure, , drop = FALSE]
index_rows$year <- as.integer(index_rows$year)
index_rows$age <- as.integer(index_rows$age)
index_rows$value <- as.numeric(index_rows$value)
keep_index <- with(index_rows,
  (survey == "UK(E&W)-BTS-Q3" & year != 2022 &
    (year < 2010 | age >= 4)) |
    (survey == "UK(E&W)-YFS" & year %in% 1987:2006 & age == 1) |
    (survey == "FR-YFS" & year %in% 1987:2023 & age == 1)
)
index_rows <- index_rows[keep_index & is.finite(index_rows$value) & index_rows$value > 0, , drop = FALSE]
inputs <- rbind(inputs[!index_measure, , drop = FALSE], index_rows)

obs <- database_to_tam_obs(
  assessment_id, inputs,
  years = years, ages = ages,
  surveys = surveys,
  sampling_times = c(
    "UK(E&W)-BTS-Q3" = 0.625,
    "UK(E&W)-YFS" = 0.625,
    "FR-YFS" = 0.625
  ),
  assumptions = source$assumptions
)

obs$catch$sd_block <- factor(ifelse(
  obs$catch$age == 1, "age_1",
  ifelse(obs$catch$age == 2, "age_2", "ages_3_plus")
))
obs$index$q_block <- ifelse(
  obs$index$survey == "UK(E&W)-BTS-Q3" & obs$index$age >= 5,
  "ages_5_6", as.character(obs$index$age)
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_block,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- ifelse(
  obs$index$survey == "UK(E&W)-BTS-Q3",
  ifelse(obs$index$age == 1, "BTS_age_1",
    ifelse(obs$index$age == 2, "BTS_age_2", "BTS_ages_3_6")),
  obs$index$survey
)
obs$maturity$spawn_time <- 0

## Model and comparison ----
age_plus_group <- 11L
N_settings <- list(process = "iid", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 3:7)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~0 + sd_block)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)

dat <- tinyAM::prepare_tam(
  data = obs, years = years, ages = ages,
  N_settings = N_settings, rec_settings = rec_settings,
  F_settings = F_settings, M_settings = M_settings,
  catch_settings = catch_settings, index_settings = index_settings
)

comparison_outputs <- source$outputs
source_n <- comparison_outputs[
  comparison_outputs$type == "population" &
    comparison_outputs$measure == "numbers_at_age", , drop = FALSE
]
source_weight <- inputs[
  inputs$type == "weight" & inputs$measure == "weight_at_age", , drop = FALSE
]
source_maturity <- inputs[
  inputs$type == "maturity" & inputs$measure == "maturity_at_age", , drop = FALSE
]
ssb_at_age <- comparison_outputs[
  comparison_outputs$type == "biomass" &
    comparison_outputs$measure == "mature_biomass_at_age", , drop = FALSE
]
key <- function(x) paste(x$year, x$age, sep = "\r")
n_index <- match(key(ssb_at_age), key(source_n))
w_index <- match(key(ssb_at_age), key(source_weight))
m_index <- match(as.character(ssb_at_age$age), as.character(source_maturity$age))
if (anyNA(c(n_index, w_index, m_index))) {
  cli::cli_abort("Official N, weight, and maturity do not cover the SSB-at-age comparison grid.")
}
ssb_at_age$value <- source_n$value[n_index] *
  source_weight$value[w_index] * source_maturity$value[m_index]
ssb_at_age$se <- ssb_at_age$lwr <- ssb_at_age$upr <- NA_real_
ssb_at_age$source_type <- "derived_source_output"
ssb_at_age$source_reference <- paste(unique(c(
  source_n$source_reference, source_weight$source_reference,
  source_maturity$source_reference
)), collapse = "; ")
ssb_at_age$notes <- paste(
  "Analysis-only SSB-at-age derived as official N-at-age × stock weight-at-age × maturity.",
  "It replaces Table 17.19 only in this comparison because its 2021–2025 age values do not sum to reported annual SSB.",
  "No age-specific uncertainty is available."
)
comparison_outputs <- rbind(
  comparison_outputs[!(
    comparison_outputs$type == "biomass" &
      comparison_outputs$measure == "mature_biomass_at_age"
  ), , drop = FALSE],
  ssb_at_age
)

comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total abundance summed over ages 1–11+.",
    reason = "Derived from the reported N-at-age over the full modeled age range; uncertainty is unavailable."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass from N-at-age and annual stock weight-at-age over ages 1–11+.",
    reason = "The database total is derived from those same inputs and has no uncertainty estimate."
  ),
  ssb = list(
    status = "matched",
    definition = "Start-of-year mature biomass from N-at-age, stock weight-at-age, and maturity over ages 1–11+.",
    reason = "These source inputs reproduce reported annual SSB within 0.02%. The age-specific table is replaced in the comparison by the same calculation because its 2021–2025 values do not sum to annual SSB."
  ),
  recruitment = list(
    status = "matched",
    definition = "Same-year abundance at age 1, in thousand fish.",
    reason = "The source recruitment age and units match; population process assumptions differ."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Common population-weighted F over ages 3–7.",
    reason = "ICES reports the arithmetic mean over ages 3–7; the dashboard recalculates the source series using N weights to match tinyAM's definition."
  ),
  M = list(
    status = "approximate",
    definition = "Fixed natural mortality input of 0.1 per year at ages 1–11+.",
    reason = "This is an input consistency check, not an independently estimated M comparison."
  )
)

background <- c(
  "### Eastern English Channel common sole: 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1982–2025; ages 1–11+, with age 11 as the plus group. | Fit the same years and ages. | The source model grid is retained. |",
  "| N and recruitment | SAM uses an IID log-N process and random-walk recruitment. | Use an IID N process, random-walk recruitment, and exponential initial abundance. | Process variances and initial-state integration are simplified. |",
  "| F and M | SAM uses a random-walk F process with some adjacent ages sharing states; M is fixed at 0.1 per year. | Use age-specific temporal random-walk F and fixed M=0.1. | SAM's age-state sharing and other process details are not reproduced. |",
  "| Catch | Combined catch-at-age, with catch SDs shared for age 1, age 2, and ages 3–11. | Retain the combined series and the same three SD groups. | The one zero catch observation is omitted by tinyAM's lognormal likelihood. |",
  "| Indices and q | UK BTS ages 1–6, UK YFS age 1, and French YFS age 1. BTS ages 1–3 from 2010 onward and all 2022 BTS values are excluded; BTS q is shared for ages 5–6. | Fit the same retained observations and q sharing. Use survey-specific SD groups, with BTS SDs for age 1, age 2, and ages 3–6. | A single 0.625 timing represents the reported 0.5–0.75 window. TinyAM does not reproduce BTS's cross-age residual correlation. |",
  "| Commercial indices | Three commercial LPUE series are fitted as aggregate fishable-stock-biomass indices. | Retain these source inputs but do not fit them. | tinyAM's age-structured index likelihood cannot represent this aggregate likelihood. |",
  "| Biology and outputs | Annual stock weights, maturity by age, fixed M, and start-of-year mature biomass. | Use source biology and compare SSB, N-at-age, F-at-age, recruitment, abundance, and biomass. | Age-specific N, F, and derived SSB have no uncertainty estimates. The SSB-at-age comparison is calculated from N × weight × maturity because the published 2021–2025 values appear misaligned. |",
  "",
  "The survey-index scale is unresolved, so the reported index values are not rescaled. Stock weights for 1982–2003 are reconstructed in the source report. The comparison retains published aggregate SSB and uses a common biological calculation for SSB-at-age; canonical Table 17.19 values remain unchanged."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) {
    arguments <- c(
      list(data = obs, years = years, ages = ages),
      list(
        N_settings = N_settings, rec_settings = rec_settings,
        F_settings = F_settings, M_settings = M_settings,
        catch_settings = catch_settings, index_settings = index_settings,
        silent = silent
      )
    )
    if (!is.null(start_par)) arguments$start_par <- start_par
    do.call(tinyAM::fit_tam, arguments)
  }
  retry_start <- function() {
    source_N <- .translation_start_surface(
      source$outputs, "population", "numbers_at_age", years, ages
    )
    source_F <- .translation_start_surface(
      source$outputs, "mortality", "fishing_mortality_at_age", years, ages
    )
    start_par <- tinyAM::make_par(dat)
    start_par$log_r0 <- log(source_N[1L, 1L])
    start_par$log_r <- log(source_N[-1L, 1L])
    start_par$log_n <- log(source_N[-1L, -1L, drop = FALSE])
    start_par$log_f <- log(source_F)
    start_par
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      assessment_id, "initialization_attempts.rds"
    )
  )
}
