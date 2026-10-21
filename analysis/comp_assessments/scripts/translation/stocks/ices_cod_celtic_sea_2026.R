source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1980:2025
ages <- 1:7
assessment_id <- source$assessment$assessment_id[[1]]

inputs <- source$inputs
inputs$value <- suppressWarnings(as.numeric(as.character(inputs$value)))
inputs <- inputs[!(inputs$type == "catch" &
  inputs$measure == "numbers_at_age" & is.na(inputs$value)), , drop = FALSE]
inputs <- inputs[!(inputs$type == "index" &
  inputs$survey == "VAST_IGFS_EVHOE" &
  (is.na(inputs$value) | inputs$value <= 0)), , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = "VAST_IGFS_EVHOE",
  sampling_times = c(VAST_IGFS_EVHOE = 0.855),
  assumptions = source$assumptions
)

obs$catch$obs_sd_group <- factor(ifelse(
  obs$catch$age == 1, "age1",
  ifelse(obs$catch$age == 2, "age2", "age3plus")
), levels = c("age1", "age2", "age3plus"))
obs$index$q_key <- factor(paste0("age", obs$index$age), levels = paste0("age", 1:3))
obs$index$obs_sd_group <- factor("VAST", levels = "VAST")

fraction_surface <- function(measure) {
  x <- source$inputs[
    source$inputs$type == "biology" & source$inputs$measure == measure,
    c("year", "age", "value")
  ]
  x$year <- as.integer(x$year)
  x$age <- as.integer(x$age)
  x$value <- as.numeric(x$value)
  x
}

key <- function(x) paste(x$year, x$age, sep = "\034")
fraction_f <- fraction_surface("fraction_F_before_spawning")
fraction_m <- fraction_surface("fraction_M_before_spawning")
maturity_key <- key(obs$maturity)
f_index <- match(maturity_key, key(fraction_f))
m_index <- match(maturity_key, key(fraction_m))
if (anyNA(f_index) || anyNA(m_index)) {
  cli::cli_abort("Spawning fractions do not cover the translated maturity surface.")
}
obs$maturity$spawn_time <- (fraction_f$value[f_index] + fraction_m$value[m_index]) / 2

age_plus_group <- 7L
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  ssb = list(
    status = "matched",
    definition = "Native SAM SSB and tinyAM mature biomass use the same ages, stock weights, maturity and zero pre-spawning mortality fractions.",
    reason = "The fitted models differ, but their reported SSB calculation uses the same source biological surfaces and timing."
  ),
  biomass = list(
    status = "matched",
    definition = "Total stock biomass across ages 1-7+ using the source stock weights.",
    reason = "The SAM and tinyAM population estimates differ; the biomass definition and age range match."
  ),
  abundance = list(
    status = "matched",
    definition = "Total abundance is the sum of N-at-age over ages 1-7+.",
    reason = "The total is derived from the official SAM N-at-age table and has no propagated uncertainty."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 recruitment in the same calendar year, in thousand fish.",
    reason = "The age and units match; the recruitment-process and initial-state parameterizations differ."
  )
)

## Model ----
background <- c(
  "### Celtic Sea cod: official 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1980-2025, ages 1-7+; age 7 is the plus group. | Fit the full 1980-2025 and 1-7+ grid. | The full fitted historical grid is retained. |",
  "| Catch | Total catch-at-age, including landings and discards; ages 1 and 2 are missing through 2003. | Retain the same available catch observations and ages. | Unavailable source cells are omitted, not filled or set to zero. |",
  "| Indices | The final SAM run uses French OTDEF Q2+3+4 aggregate LPUE and the combined Irish IGFS/French EVHOE Q4 VAST index at ages 1-3. | Use the age-specific VAST index only. | The aggregate commercial LPUE likelihood is not supported by tinyAM's age-specific index model. The report's VAST age-4 series is not used by the fitted SAM object. |",
  "| N and recruitment | SAM uses age-1 recruitment and a shared age-2+ survival process. | Use an exponential initial state, IID N deviations and a random-walk recruitment process. | tinyAM does not reproduce the source variance groups or complete initial-state treatment. |",
  "| F | The source has six age-specific F states, AR(1)-correlated process increments across ages, and a shared age-6/7+ state; Fbar is ages 2-5. | Use a correlated-random-walk F process across ages 1-7, with mean level aligned to ages 2-5. | The exact source state sharing, selectivity and covariance structure are simplified. |",
  "| M | Fixed age-specific M: 0.501, 0.330, 0.264, 0.233, 0.211, 0.211 and 0.211 per year for ages 1-7+. | Use the same fixed M input. | The source M surface is retained; M is not estimated. |",
  "| Catchability and observation error | VAST q is separate at ages 1, 2 and 3; its residuals use AR correlation across ages. Catch SD groups are age 1, age 2 and ages 3-7. | Use separate q values by VAST age, one VAST SD, and the three catch SD groups. | tinyAM does not model the VAST cross-age residual correlation. |",
  "| Biology and biomass | The SAM run uses annual stock weights, annual catch/landings/discard weights and constant maturity; F and M fractions before spawning are zero. | Retain source stock weights and maturity, and set spawning time to zero. | The same SSB timing and biological inputs are used. |",
  "| Comparisons | SAM reports N, F, SSB, total biomass, age-1 recruitment and Fbar. | Compare those metrics over the full fitted period; total abundance is derived from N-at-age. | Total abundance has no uncertainty because it is calculated from N-at-age without covariance information. |",
  "",
  "The source report and native SAM fit object are cached locally. The stock-annex PDF was attempted but returned HTTP 401; the final SAM object supplies the fitted settings and numerical inputs. The package-default fit was attempted first. Its optimizer failed with an NA/NaN gradient; one retry initialized N, recruitment and F states from aligned official estimates. That retry returned false convergence (code 1), so no dashboard is available yet."
)

fit_model <- function(...) tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + obs_sd_group),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~1),
    silent = silent,
    ...
  )

fit_start <- function() {
  model_dat <- tinyAM::prepare_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + obs_sd_group),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~1)
  )
  source_N <- .translation_start_surface(
    source$outputs, "population", "numbers_at_age", years, ages
  )
  source_F <- .translation_start_surface(
    source$outputs, "mortality", "fishing_mortality_at_age", years, ages
  )
  start <- tinyAM::make_par(model_dat)
  start$log_r0 <- log(source_N[1L, 1L])
  start$log_r <- log(source_N[match(names(start$log_r), rownames(source_N)), 1L])
  start$log_n <- log(source_N[rownames(start$log_n), colnames(start$log_n), drop = FALSE])
  start$log_f <- log(source_F[rownames(start$log_f), colnames(start$log_f), drop = FALSE])
  start
}

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = fit_start,
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache",
      source$assessment$assessment_id[[1]],
      "initialization_attempts.rds"
    )
  )
}
