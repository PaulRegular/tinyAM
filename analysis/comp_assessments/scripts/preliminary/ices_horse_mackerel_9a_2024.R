## Observations ----

years <- 1992:2023
ages <- 0:11
ibts_survey <- "Portuguese-Spanish IBTS"

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = ibts_survey,
  sampling_times = c("Portuguese-Spanish IBTS" = 0.8),
  assumptions = source$assumptions
)

obs$index <- obs$index[obs$index$age > 0 & obs$index$obs > 0, , drop = FALSE]
obs$index$q_block <- factor(obs$index$age, levels = 1:11)
obs$index$sd_block <- 1
obs$index$relative_sd <- 0.30

## Background and comparisons ----

age_plus_group <- 11L
comparison_outputs <- source$outputs
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  recruitment = list(
    status = "approximate",
    definition = "Age-0 abundance in the same calendar year, reported in thousands of fish",
    reason = "The reported age matches; tinyAM uses an estimated Beverton-Holt curve and IID deviations instead of the source equilibrium initialization and stock-recruitment parameterization."
  ),
  ssb = list(
    status = "non_equivalent",
    definition = "Native Stock Synthesis SSB versus tinyAM start-year SSB",
    reason = "The source spawning-time calculation is not fully documented in the recovered files, and age-specific source abundance is unavailable for a common-definition calculation."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Source F over ages 1-10 versus tinyAM population-weighted F over ages 1-10",
    reason = "The annual age-specific source N and F values needed to recalculate a matching weighted mean are unavailable."
  )
)

background <- c(
  "### Southern horse mackerel: 2024 ICES Stock Synthesis assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | Annual 1992-2023, ages 0-10 and 11+, sexes combined. | Fit the same years and ages, with age 11 as the plus group. | The source age and year grid is retained. |",
  "| N and recruitment | Beverton-Holt recruitment, sigmaR 0.74, steepness 0.8, and equilibrium initial age structure based on a 30,000-tonne equilibrium catch. | Beverton-Holt mean with IID recruitment deviations and exponential initial abundance. | The curve parameters, process bias treatment, and equilibrium initialization are not reproduced exactly. |",
  "| F | Double-normal fishery selectivity with two time-varying parameters. | Age-specific F follows a random walk. | Selectivity and its temporal parameterization are simplified. |",
  "| M and biology | Fixed Gislason M by age, time-invariant female maturity, and stock weights equal to catch weights. | Retain the reported M, maturity, and weight-at-age surfaces. | M is fixed in both representations; tinyAM uses the female maturity ogive for the combined-sex population. |",
  "| Catch | Combined catch biomass and numbers-at-age; the source uses catch CV 0.05 and a Dirichlet-multinomial age likelihood. | Fit the reported direct catch numbers-at-age with a lognormal error model; retain source catch biomass in the database only. | The source biomass-plus-composition likelihood and age-reading error are not represented. |",
  "| Survey index | Portuguese-Spanish October IBTS biomass and age composition; ages 1-11+, with a nominal CV of 0.30. | Fit the recoverable stratified mean abundance-at-age values at ages 1-11+, with q estimated by age and the source 0.30 weight supplied. | TinyAM fits the published age-specific index directly, rather than aggregate biomass and age composition. Reported zero cells are omitted because the fitted log-index likelihood cannot use zero observations. |",
  "| Other indices | Triennial DEPM SSB and commercial biomass CPUE are also in the official model. | Keep the DEPM observations in the database, but omit them and CPUE from the fit. | tinyAM cannot translate aggregate SSB to age-specific predictions without age composition; the numeric CPUE series was not recovered. |",
  "| SSB and Fbar | The report gives aggregate SSB and F(1-10), with CVs, but no numerical age-specific N or F. | Compare recruitment approximately; show native SSB and Fbar as non-equivalent context. | No age-specific source states are invented to derive total abundance, F-at-age, or common-definition SSB/Fbar. |",
  "",
  "The fit uses package-default initialization and no official values as parameter starts. It is a structural approximation based on the recoverable age-specific series, not a replication of the source likelihood."
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
    rec_settings = list(process = "iid", mu_form = ~ bh(ssb)),
    F_settings = list(process = "rw"),
    M_settings = list(
      process = "off",
      mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~0 + q_block,
      sd_form = ~0 + sd_block,
      sd_supplied = ~relative_sd
    ),
    silent = silent
  )
}
