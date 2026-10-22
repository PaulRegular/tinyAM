## Observations ----
years <- 1955:2026
ages <- 3:14
surveys <- c("IS-SMB", "IS-SMH")
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("IS-SMB" = 0.2, "IS-SMH" = 0.8),
  exclude_index_years = list("IS-SMH" = 2011L),
  assumptions = source$assumptions
)
obs$index <- obs$index[!(obs$index$survey == "IS-SMH" & obs$index$age > 13), , drop = FALSE]
obs$index$q_key <- interaction(obs$index$survey, obs$index$age, drop = TRUE, lex.order = TRUE)
obs$index$sd_block <- factor(obs$index$survey)

## Comparison ----
age_plus_group <- NULL
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total abundance across ages 3-14 in thousand fish.",
    reason = "The reference series is derived from official N-at-age over the same ages as the tinyAM fit."
  ),
  ssb = list(
    status = "approximate",
    definition = "Mature biomass from fitted N, stock weights and maturity over ages 3-14.",
    reason = "tinyAM cannot use the separate spawning-weight surface or the source pre-spawning mortality fractions."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-3 recruitment in thousand fish.",
    reason = "The assessment reports age-3 recruitment, but its cohort-migration and recruitment parameterization is not reproduced."
  ),
  F = list(
    status = "approximate",
    definition = "Annual F-at-age for ages 3-14.",
    reason = "A temporal random walk does not reproduce the source selectivity periods or exact covariance structure."
  ),
  M = list(
    status = "matched",
    definition = "Fixed M = 0.2 per year.",
    reason = "The same fixed age-specific M is supplied to tinyAM."
  )
)

background <- c(
  "### Icelandic cod: 2026 MFRI assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | Catch-at-age covers 1955-2025 at ages 3-14; the March index and population estimates extend to 2026. | Fit 1955-2026 at ages 3-14. | The 2026 year has no catch-at-age; the age-14 terminal group is not treated as a plus group. |",
  "| Catch | Numbers-at-age are aggregated across 20 métiers. | Retain the published aggregated catch-at-age. | The component métiers and their separate likelihood weights are unavailable. |",
  "| Surveys | March (IS-SMB) covers ages 1-14; October (IS-SMH) covers ages 3-13 and has no 2011 survey. | Use ages 3-14 for March and 3-13 for October; use 0.2 and 0.8 sampling times from the historical MUPPET configuration. | Current production timing and correlated residual settings are not fully verified; zero March-index cells do not contribute to tinyAM's lognormal likelihood. |",
  "| N and recruitment | The source estimates age-3 recruitment and age-specific population numbers. | Use deterministic cohort survival and a random-walk recruitment process. | Source migration events and exact process configuration are not reproduced. |",
  "| F | MUPPET describes a temporal random walk with selectivity periods. | Use a random walk through time at each age and define mean F over ages 5-10. | The selectivity blocks and source covariance are not represented. |",
  "| M and biology | M is fixed at 0.2; the source reports stock, spawning-stock and catch weights, and maturity. | Use fixed M, stock weights and maturity. | tinyAM uses one weight surface for population and SSB and does not apply source pre-spawning F/M fractions. |",
  "| Outputs | Official tables report SSB, B4+, recruitment, N-at-age and F-at-age. | Compare SSB, age-3 recruitment, N/F at age and derived total abundance over ages 3-14. | Age-specific uncertainty is unavailable; B4+ is retained as a native output but is not equivalent to tinyAM total biomass. |",
  "",
  "The 2026 report narrative gives SSB of 371,000 tonnes, while its machine-readable summary table gives 377,277 tonnes. The database retains the table value and flags the difference. The report also describes pre-1985 stock weights and maturity as held at 1985 values, while the published tables show varying values; those table values are preserved unchanged. The exact fitted 2026 MUPPET object and full observation covariance were not recovered."
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
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 5:10),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block),
    silent = silent
  )
}
