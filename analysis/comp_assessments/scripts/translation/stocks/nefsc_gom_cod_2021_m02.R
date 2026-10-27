## Observations ----
assessment_id <- source$assessment$assessment_id[[1L]]
years <- 1982:2019
ages <- 1:9

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  sampling_times = c("NEFSC spring" = 0.25, "NEFSC fall" = 0.75),
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0
obs$index$q_key <- interaction(obs$index$survey, obs$index$age, drop = TRUE)
obs$catch$relative_sd <- sqrt(log1p(0.05^2))

## Comparison and background ----
age_plus_group <- 9
comparison_outputs <- source$outputs[
  source$outputs$measure %in% c("SSB", "recruitment"), , drop = FALSE
]
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Reported annual SSB using January 1 stock weights and the digitized female maturity curve.",
    reason = "The source does not publish the numeric M02 maturity vector or an explicit spawning-time fraction. tinyAM uses the source January 1 weight convention and sets spawning time to 0; its population estimates also come from a different likelihood and model."
  ),
  recruitment = list(
    status = "matched",
    definition = "Annual age-1 recruitment in thousand fish.",
    reason = "Both assessments report age-1 recruitment in the same calendar year; tinyAM uses a different process likelihood."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "The source reports fully selected F, while tinyAM summarizes age-specific F over selected ages.",
    reason = "The source report does not identify the fully selected age range and no numeric official F-at-age surface was recovered."
  ),
  abundance = list(
    status = "unavailable",
    definition = "Total beginning-of-year abundance.",
    reason = "The source report does not publish numeric fitted N-at-age, so total abundance cannot be derived."
  )
)

background <- c(
  "### Gulf of Maine Atlantic cod: official 2021 M02 management-track model",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1982–2019, ages 1–9+, with recruitment at age 1. | Fit 1982–2019, ages 1–9, with age 9 as the plus group. | The reported model grid is retained. |",
  "| N and recruitment | ASAP age-structured model with annual age-1 recruitment deviations. | Deterministic cohort survival with IID recruitment around a constant mean and exponential initial abundance. | tinyAM does not reproduce the source initial-state integration or recruitment likelihood. |",
  "| F and M | One combined fishery, three selectivity periods, fixed M=0.2 per year. | AR1 F process with logistic age means; fixed source M. | Selectivity blocks and the source catch likelihood are simplified; M is an input check, not an estimated comparison. |",
  "| Catch | Total fishery catch-at-age after gear-specific commercial discard mortality and 15% recreational release mortality; source CV 0.05. | One total catch-at-age stream with the source CV converted to log-residual SD. | tinyAM uses one lognormal age-specific likelihood, not the source's full likelihood structure. |",
  "| Surveys and q | NEFSC spring and fall, plus MADMF spring; survey catchability is time-constant. The source also applies survey-specific CV reweighting. | Use NEFSC spring and fall abundance-at-age with survey-specific age-logistic q and estimated survey-level observation SD. | Published abundance-at-age estimates stand in for the source's index and composition likelihoods. MADMF age values and annual survey errors were not recovered. Sampling times 0.25 and 0.75 are seasonal approximations. |",
  "| Biology and SSB | Female time-series-average maturity and January 1 / stock weights; exact numeric M02 maturity values are not published. | Use weights from Table 8.3 and a rounded digitization of the female maturity curve; set spawn time to 0. | The digitized maturity values and January 1 timing choice are approximate. |",
  "| Outputs | SSB, age-1 recruitment, and fully selected F are tabulated for 2012–2019; fitted N-at-age and F-at-age are shown only in figures. | Compare SSB and recruitment over 2012–2019. | F-at-age, total abundance, and uncertainty cannot be recovered. The fully selected F age range is not specified, so it is not compared to tinyAM F summaries. |",
  "",
  "The 2021 combined Gulf of Maine unit is historical. The 2023 research track replaced it with four stock units; none of the old-unit inputs or outputs are assigned to those successors."
)

## Model ----
source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
model_settings <- list(
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~ 1),
  F_settings = list(process = "ar1", mu_form = ~ logistic(age)),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
  catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd),
  index_settings = list(
    q_form = ~ 0 + survey + logistic(age, by = survey),
    sd_form = ~ 0 + survey
  )
)
fit_model <- function(...) do.call(tinyAM::fit_tam, c(
  list(data = obs, years = years, ages = ages, silent = silent),
  model_settings,
  list(...)
))

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(
    fit_model,
    retry_start = function() cli::cli_abort(
      "No official age-specific N or F estimates are available for a source-based retry."
    ),
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache", assessment_id,
      "initialization_attempts.rds"
    )
  )
}
