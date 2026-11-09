## Observations ----

years <- 1986:2024
ages <- 1:7
survey <- "INIDEP global winter survey"

model_inputs <- source$inputs
is_weight <- !is.na(model_inputs$type) & model_inputs$type == "weight" &
  !is.na(model_inputs$measure) & model_inputs$measure == "weight_at_age"
keep <- !is_weight | (!is.na(model_inputs$season) & model_inputs$season == "1 July")
model_inputs <- model_inputs[keep, , drop = FALSE]

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  model_inputs,
  years = years,
  ages = ages,
  surveys = survey,
  sampling_times = stats::setNames(0.5, survey),
  assumptions = source$assumptions
)

comparison_outputs <- source$outputs
ssb_rows <- comparison_outputs$type == "biomass" &
  comparison_outputs$measure == "spawning_biomass"
comparison_outputs$measure[ssb_rows] <- "SSB"
age_plus_group <- 7
comparison_aggregates <- c("ssb", "recruitment", "biomass")
comparison_definitions <- list(
  N = list(
    status = "approximate",
    definition = "Official beginning-of-year numbers at ages 1-6 and 7+",
    reason = "Population dynamics and the source age-structured likelihood are simplified."
  ),
  ssb = list(
    status = "approximate",
    definition = "Official spawning biomass, reported as spawning_biomass",
    reason = "tinyAM uses the source July weight, maturity and fixed M schedules, but simplifies the population and spawning calculations."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Official recruitment at age 1",
    reason = "tinyAM uses an IID recruitment process rather than the source model's recruitment structure."
  ),
  biomass = list(
    status = "approximate",
    definition = "Official total biomass",
    reason = "The population dynamics and biomass calculation are simplified in tinyAM."
  )
)

## Model ----

background <- c(
  "### Southern Argentine hake, south of 41 degrees S: INIDEP 2025 assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1986-2024; ages 1-6 and 7+. | Same years and age groups. | 2025 values are preliminary projections and are outside the fitted historical period. |",
  "| Catch | Estimated catch numbers-at-age combine reported catch, discards and bycatch across sources. | One aggregate catch series using the reported numbers-at-age. | Fleet-specific age catches and the source composition likelihood are not reproduced. |",
  "| Survey | Winter global-survey numbers-at-age in selected years; other inputs include juvenile indices, reproductive biomass indices and CPUE. | Winter global-survey numbers-at-age only. | Aggregate or age-2 indices are omitted because they do not provide the same age-specific observation series. |",
  "| Biology and timing | Fixed M = 0.3 per year; population weights at 1 January and 1 July; source maturity schedule. | Fixed source M and maturity; the 1 July weight vector is used for the annual model. Survey timing is set to 0.5. | tinyAM uses one annual weight surface; the report describes the survey as winter but gives no precise sampling fraction. |",
  "| Processes | Statistical catch-at-age model with source recruitment and fleet-specific fishing mortality. | IID recruitment, an RW process on log F with independent increments across ages, and no additional N process. | tinyAM does not reproduce the source model's recruitment, selectivity and fleet structure. |",
  "| Outputs | Annual N-at-age, biomass, spawning biomass, recruitment and a scalar F summary are reported. | Compare N-at-age, total biomass, SSB and recruitment. | Age-specific F is unavailable. The source spawning_biomass field is mapped to SSB for the dashboard. |",
  "",
  "The default fit uses package starting values. If it fails numerical checks, one retry uses aligned official N-at-age values only as starting values; unavailable parameters retain package defaults. Starting values do not constrain the estimates."
)

N_settings <- list(process = "off", init = "exp")
rec_settings <- list(process = "iid", mu_form = ~1)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 3:6)
M_settings <- attr(obs, "translation")$M$settings
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~1, sd_form = ~1)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)

  fit_model <- function(start_par = NULL) {
    tinyAM::fit_tam(
      data = obs,
      years = years,
      ages = ages,
      N_settings = N_settings,
      rec_settings = rec_settings,
      F_settings = F_settings,
      M_settings = M_settings,
      catch_settings = catch_settings,
      index_settings = index_settings,
      silent = silent,
      start_par = start_par
    )
  }

  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      model_dat <- tinyAM::prepare_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = N_settings,
        rec_settings = rec_settings,
        F_settings = F_settings,
        M_settings = M_settings,
        catch_settings = catch_settings,
        index_settings = index_settings
      )
      source_start_parameters(source, model_dat)
    },
    attempt_file = file.path(
      .assessment_cache_root(translation_tier, database = source),
      source$assessment$assessment_id[[1]],
      "initialization_attempts.rds"
    )
  )
}
