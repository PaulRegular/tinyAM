## Observations ----
years <- 1977:2025
ages <- 0:8

obs <- database_to_tam_obs(
  source$assessment$assessment_id[[1]],
  source$inputs,
  years = years,
  ages = ages,
  surveys = "GRAHS",
  sampling_times = c(GRAHS = 0.675),
  assumptions = source$assumptions
)

obs$catch$obs_sd_group <- factor(ifelse(obs$catch$age == 0, "age0", "age1plus"))
obs$index$q_key <- factor(ifelse(
  obs$index$age == 1, "age1",
  ifelse(obs$index$age == 2, "age2", "age3plus")
), levels = c("age1", "age2", "age3plus"))
obs$index$obs_sd_group <- factor(ifelse(
  obs$index$age == 1, "age1",
  ifelse(obs$index$age == 8, "age8plus", "age2to7")
), levels = c("age1", "age2to7", "age8plus"))

fraction_surface <- function(measure) {
  x <- source$inputs[
    source$inputs$type == "biology" & source$inputs$measure == measure,
    c("year", "age", "value")
  ]
  x$year <- as.integer(x$year)
  x$age <- as.integer(x$age)
  x
}

fraction_f <- fraction_surface("fraction_F_before_spawning")
fraction_m <- fraction_surface("fraction_M_before_spawning")
key <- function(x) paste(x$year, x$age, sep = "\034")
maturity_key <- key(obs$maturity)
f_index <- match(maturity_key, key(fraction_f))
m_index <- match(maturity_key, key(fraction_m))
if (anyNA(f_index) || anyNA(m_index)) {
  cli::cli_abort("Spawning fractions do not cover the translated maturity surface.")
}
obs$maturity$spawn_time <- (fraction_f$value[f_index] + fraction_m$value[m_index]) / 2

comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Official N-at-age with common source weights and maturity, applying an effective within-year spawning fraction.",
    reason = "SAM applies separate F and M fractions before spawning; tinyAM has one spawning-time fraction for total mortality. The arithmetic mean of the source fractions is used."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-0 recruitment in thousand fish over 1977-2025.",
    reason = "The age and units match, but tinyAM's recruitment and initial-state process treatment simplifies the SAM random-walk model."
  )
)
age_plus_group <- 8

## Model ----
background <- c(
  "### Gulf of Riga herring: official 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1977-2025, ages 0-8+; age 8 is the plus group. | Fit the full 1977-2025 and 0-8+ grid. | The historical model grid is retained. |",
  "| N and recruitment | Age-0 recruitment random walk; separate recruitment and shared age-1-8 survival process variances. | IID N process with exponential initial abundance and random-walk recruitment. | tinyAM uses a different initial-state and variance-sharing structure. |",
  "| F | Correlated temporal random walk; ages 7 and 8+ share a state; Fbar ages 2-6. | Correlated random-walk F with independent age states. | tinyAM does not reproduce the exact SAM age-state sharing and age-specific process-variance groups. |",
  "| M | Fixed at 0.20 per year, except 0.25 in 1979-1983. | Use the source M-at-age surface. | The numerical M assumptions are retained. |",
  "| Catch and index | Lognormal catch-at-age and August GRAHS acoustic index; q at age 1, age 2, and ages 3-8; observation SD groups follow the SAM keys. | Retain the same observations, GRAHS timing, q groups, and catch/index SD groups. | Acoustic residual correlation across ages is not represented. |",
  "| Biology and SSB | Annual stock weight equals catch weight; maturity is fixed through 1994 and GAM-smoothed thereafter; SAM uses separate F/M pre-spawning fractions. | Retain weight and maturity, with a single effective spawning fraction equal to the mean of source F and M fractions. | tinyAM uses one fraction of total mortality and cannot represent separate F and M timing exactly. |",
  "| Comparisons | SAM reports N-at-age, F-at-age, SSB, recruitment and Fbar; 95% intervals are available for SSB, recruitment and Fbar. | Compare reported trajectories using the source units and age grid. | tinyAM uncertainty and F process structure differ from SAM; source output uncertainty remains native. |",
  "",
  "The fit uses package-default starting values. The acoustic index is kept on its native scale because the report gives no physical multiplier."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + obs_sd_group),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + obs_sd_group
    ),
    silent = silent
  )
}
