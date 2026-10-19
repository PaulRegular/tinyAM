## Observations ----
years <- 1969:2025
ages <- 1:10

obs <- database_to_tam_obs(
  source$assessment$assessment_id[[1L]],
  source$inputs,
  years = years,
  ages = ages,
  surveys = character(),
  assumptions = source$assumptions
)

## Comparisons and model ----
age_plus_group <- 10
comparison_outputs <- source$outputs
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  N = list(
    status = "approximate",
    definition = "Numbers at age, with age 10 as 10+",
    reason = "The source and tinyAM use the same age groups, but differ in population and observation processes."
  ),
  F = list(
    status = "approximate",
    definition = "Annual fishing mortality at age, with age 10 as 10+",
    reason = "The source selectivity and fishing-mortality structure are simplified in tinyAM."
  ),
  M = list(
    status = "matched",
    definition = "Fixed source natural mortality supplied to tinyAM",
    reason = "M is fixed at 0.28 per year in the source assessment; agreement is by construction, not an estimated comparison."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-1 abundance in the same calendar year",
    reason = "The reported age definition matches, while recruitment and population processes differ."
  ),
  abundance = list(
    status = "approximate",
    definition = "Sum of numbers at ages 1-10+",
    reason = "The source total is derived from N-at-age; population processes differ."
  ),
  ssb = list(
    status = "approximate",
    definition = "Beginning-of-year spawning biomass using shared source weights and maturity",
    reason = "The comparison uses a common definition derived from source N-at-age because tinyAM does not separately apply the source F- and M-before-spawning fractions."
  ),
  biomass = list(
    status = "approximate",
    definition = "Total biomass from N-at-age and shared source weights",
    reason = "The common-definition source value is derived from source N-at-age and the same weight surface used by tinyAM."
  )
)

background <- c(
  "### Northern contingent Atlantic mackerel: DFO 2025 assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1969-2025, ages 1-10+, with age 10 as the plus group. | 1969-2025, ages 1-10, with age 10 as the plus group. | The reported model grid is retained. |",
  "| N and recruitment | CCAM uses Beverton-Holt recruitment, IID recruitment deviations, and IID cohort-survival residuals; ages 2-10+ share one N-process variance. | Beverton-Holt recruitment with IID deviations, IID survival residuals, and exponential initial abundance. | The recruitment and survival process forms are retained, but tinyAM does not reproduce the source variance grouping or initial-state integration. |",
  "| F | One annual log-F level follows a random walk, with age-specific selectivity held constant; ages 5-10+ share the terminal selectivity key. | Age-specific log-F states follow a random walk with correlated age increments. | Correlated increments approximate the shared annual F change, but do not impose the source's exact selectivity pattern. |",
  "| M and biology | M is fixed at 0.28 per year. Annual weight and maturity are supplied; the source also provides separate F- and M-before-spawning fractions of 0.47. | Use the same fixed M, annual stock weights and maturity. | tinyAM has one spawning-time adjustment using combined Z, so the source's separate F and M fractions are not applied. |",
  "| Catch | The source fits bounded total catch and catch-at-age composition with source-specific likelihoods. | Fit the unchanged catch-at-age observations directly with lognormal error. | The direct catch likelihood is an approximation; zero cells are treated as missing by the log observation model. |",
  "| Egg-production index | An aggregate total egg production index is fitted using maturity, female proportion, fecundity and mortality to sampling time. | Omit the index and fit the catch-only approximation. | tinyAM does not calculate this aggregate fecundity-weighted prediction; assigning it ages or treating eggs as fish abundance would change the observation. |",
  "| SSB and biomass | The official run reports native SSB and total biomass with uncertainty. | Compare common definitions derived from source N-at-age using the same annual weight and maturity surfaces as tinyAM. | Derived comparisons do not have uncertainty; native official estimates remain available separately. |",
  "",
  "The source's reported M is fixed, so any M agreement is by construction. The fitted model is a structural approximation, not a close likelihood replication. Because the aggregate egg index is omitted, catch-at-age is the only observation series; population scale and F may therefore be weakly identified. The fit call uses package-default starting values."
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
    rec_settings = list(process = "iid", mu_form = ~ bh(ssb)),
    F_settings = list(
      process = "cor_rw", mu_form = NULL,
      mean_ages = 5:10
    ),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0, sd_form = ~0),
    silent = silent
  )
}
