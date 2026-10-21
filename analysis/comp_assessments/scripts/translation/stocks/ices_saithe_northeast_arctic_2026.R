## Observations ----
years <- 1960:2025
ages <- 3:12
surveys <- c(
  "FLT13 Norwegian acoustic survey",
  "FLT14 Norwegian acoustic survey"
)

obs <- database_to_tam_obs(
  source$assessment$assessment_id[[1]],
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)

obs$index$q_age <- factor(
  ifelse(obs$index$age <= 5, paste0("age", obs$index$age), "age6_7"),
  levels = c("age3", "age4", "age5", "age6_7")
)
obs$index$q_key <- interaction(obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$relative_precision_weight[!is.finite(
  obs$index$relative_precision_weight
)] <- 1
obs$index$relative_sd <- 1 / sqrt(obs$index$relative_precision_weight)

age_plus_group <- 12L
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total abundance is the sum of N-at-age over ages 3-12+.",
    reason = "The total is derived from the official age-specific abundance surface and has no propagated uncertainty."
  ),
  biomass = list(
    status = "matched",
    definition = "Total stock biomass across ages 3-12+ using official stock weights.",
    reason = "The assessment and tinyAM use the same age range and stock-weight surface."
  ),
  ssb = list(
    status = "matched",
    definition = "Spawning-stock biomass at the beginning of the year using official abundance, stock weights and maturity.",
    reason = "The source model has zero fishing and natural mortality before spawning, so its SSB uses the same annual state and biological surfaces."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-3 recruitment in the same calendar year, in thousand fish.",
    reason = "Age and units match; tinyAM uses an estimated Beverton-Holt relationship with IID deviations."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Fbar over ages 4-7.",
    reason = "The source reports an arithmetic mean while tinyAM reports an abundance-weighted mean."
  )
)

## Background and comparisons ----
background <- c(
  "### Northeast Arctic saithe: official 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Fitted years 1960-2025; ages 3-12+, with recruitment at age 3. | Fit the same years and ages, with age 12 as the plus group. | The separate 2026 forecast is not part of the historical fit. |",
  "| N and recruitment | Beverton-Holt recruitment with a three-year SSB lag; separate recruitment and survival variances; no initial process density. | IID recruitment deviations around a three-year-lag Beverton-Holt curve, IID survival deviations and free initial abundance. | tinyAM retains the broad process structure but not SAM's variance grouping and initial-state treatment. |",
  "| F | Log-F random walk with AR(1) correlation across age increments; ages 3-7 have separate states and ages 8-12 share one state. Fbar covers ages 4-7. | Correlated log-F random walks across ages, with age-specific states and one shared process SD; Fbar covers ages 4-7. | The shared state for ages 8-12 and other SAM-specific process details are not represented. |",
  "| M | Fixed M=0.2 per year at ages 3-12. | Use the same fixed annual M surface. | M is an input, not an estimated comparison. |",
  "| Catch | One catch-at-age series in thousand fish, with a shared lognormal SD and a prediction-variance link. | Retain the catch-at-age observations and estimate one lognormal SD. | The prediction-variance link is not represented. |",
  "| Acoustic indices | Two Norwegian acoustic surveys, ages 3-7, at sampling time 0.8; FLT13 errors are independent and FLT14 errors are AR(1) across age. Relative precision weights are available for part of FLT14. | Retain both surveys, source timing, source q-sharing by age, and relative precision weights; estimate one SD per survey. | The FLT14 cross-age error correlation is not represented. Unweighted rows use the source model's unit multiplier. |",
  "| Biology and SSB | Annual stock and catch weights, maturity, and spawning fractions are supplied; no mortality occurs before spawning. | Retain the source annual biological surfaces and compare SSB at the same beginning-of-year timing. | Native SSB uncertainty is retained only for the official reference. |",
  "| Other outputs | SAM reports N-at-age, F-at-age, fitted observations, q, recruitment, SSB, total biomass and arithmetic Fbar. | Compare available age-specific and aggregate outputs; derive total abundance from N-at-age. | Derived abundance has no propagated uncertainty; Fbar definitions differ. |",
  "",
  "The source age-12 row represents the 12+ group. The acoustic index scale is retained as reported because its physical unit is unresolved. The raw StoX CV table is preserved in the database, but its year-age pairing with native model precision weights is unresolved; the fit uses only the weights aligned to the fitted observations.",
  "The fit first uses tinyAM's normal starting values. If numerical checks fail, one retry uses official N-at-age and F-at-age as starting values only; the observations and model settings stay the same."
)

## Model ----
source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
fit_model <- function(...) tinyAM::fit_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = list(process = "iid", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ bh(ssb, lag = 3)),
  F_settings = list(process = "cor_rw", mean_ages = 4:7),
  M_settings = list(
    process = "off", mu_form = NULL, mu_supplied = ~M_assumption
  ),
  catch_settings = list(sd_form = ~1),
  index_settings = list(
    q_form = ~0 + q_key,
    sd_form = ~0 + survey,
    sd_supplied = ~relative_sd
  ),
  silent = silent,
  ...
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      dat <- tinyAM::prepare_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = list(process = "iid", init = "free"),
        rec_settings = list(process = "iid", mu_form = ~ bh(ssb, lag = 3)),
        F_settings = list(process = "cor_rw", mean_ages = 4:7),
        M_settings = list(
          process = "off", mu_form = NULL, mu_supplied = ~M_assumption
        ),
        catch_settings = list(sd_form = ~1),
        index_settings = list(
          q_form = ~0 + q_key,
          sd_form = ~0 + survey,
          sd_supplied = ~relative_sd
        )
      )
      start <- tinyAM::make_par(dat)
      source_N <- .translation_start_surface(
        source$outputs, "population", "numbers_at_age", years, ages, 1000
      )
      source_F <- .translation_start_surface(
        source$outputs, "mortality", "fishing_mortality_at_age", years, ages
      )
      start$log_r0 <- log(source_N[1L, 1L])
      start$log_n0 <- log(source_N[1L, -1L])
      start$log_r <- log(source_N[-1L, 1L])
      start$log_n <- log(source_N[-1L, -1L, drop = FALSE])
      start$log_f <- log(source_F)
      start
    },
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache",
      source$assessment$assessment_id[[1]], "initialization_attempts.rds"
    )
  )
}
