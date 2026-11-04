## Observations ----

assessment_id <- source$assessment$assessment_id[[1]]
years <- 1960:2020
ages <- 5:15
surveys <- c(
  "Canadian spring RV bottom-trawl survey",
  "Canadian autumn RV bottom-trawl survey",
  "EU-Spain Divisions 3NO spring bottom-trawl survey"
)
sampling_times <- c(
  "Canadian spring RV bottom-trawl survey" = 0.38,
  "Canadian autumn RV bottom-trawl survey" = 0.84,
  "EU-Spain Divisions 3NO spring bottom-trawl survey" = 0.38
)

inputs <- source$inputs
inputs$value <- suppressWarnings(as.numeric(as.character(inputs$value)))
inputs <- inputs[!(inputs$type %in% c("catch", "index") &
  (is.na(inputs$value) | inputs$value <= 0)), , drop = FALSE]

obs <- database_to_tam_obs(
  assessment_id, inputs,
  years = years, ages = ages,
  surveys = surveys, sampling_times = sampling_times,
  assumptions = source$assumptions
)
obs$index$q_key <- factor(paste(obs$index$survey, obs$index$age, sep = "_"))

age_plus_group <- 15L
comparison_aggregates <- character()
comparison_outputs <- source$outputs[0, , drop = FALSE]

## Model ----

N_settings <- list(process = "off", init = "free")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "rw", mu_form = NULL, mean_ages = 6:9)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption)
catch_settings <- list(sd_form = ~0, sd_supplied = ~ I(0.01))
index_settings <- list(
  q_form = ~ 0 + q_key,
  sd_form = ~0,
  sd_supplied = ~ I(1)
)

background <- c(
  "### NAFO Divisions 3LNO American plaice: preliminary translation of the historical 2021 ADAPT assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Assessment period | The 2026 record is interim monitoring, not a full assessment. The separate 2021 ADAPT run supplies detailed age-specific inputs through 2020. | Fit the historical 2021 ADAPT input series, 1960–2020, ages 5–15+. | The current monitoring record is not used as if it contained fitted population states. |",
  "| Catch | Commercial numbers-at-age, treated as exact. | Direct positive catch-at-age with a fixed log-scale SD of 0.01. | tinyAM requires a finite observation SD; this is a tight approximation to exact catch, not a source error estimate. |",
  "| Indices and q | Canadian spring, Canadian autumn and EU-Spain spring age-specific indices; q is survey- and age-specific. | Retain all three series, with survey-age q and fixed log-residual SD of 1. | Timing values 0.38 and 0.84 are seasonal midpoints from the Canadian survey date ranges; 0.38 approximates the EU spring series. Zero indices are omitted. |",
  "| N, recruitment and F | ADAPT population estimates through 2020; terminal survivor parameters for 2021. Historical absolute estimates were considered unreliable. | Deterministic cohort survival, age-5 recruitment random walk, and age-specific F random walk. | tinyAM does not reproduce ADAPT. No historical N/F surface is available to compare or use for a source-based retry. |",
  "| M and biology | M=0.2 per year, except 0.53 in 1989–1996; stock weights and maturity are annual. | Use the recorded year-age M surface and original stock weights/maturity. | Fixed M agreement is by construction, not evidence of estimated-M performance. |",
  "| Outputs | The accessible numerical output is a small set of terminal survivor parameter estimates with SEs; annual bias-adjusted N, F and SSB are not tabulated. | Fit for structural exploration; no annual Official–tinyAM comparison is claimed. | Survivor parameters are not treated as bias-adjusted population abundance, and no annual quantities are reconstructed from figures. |",
  "",
  "The 2021 source warns that its absolute population estimates are unreliable. A fit or trend comparison would not validate absolute stock size. The 2026 interim-monitoring record remains a separate current record and is not a model-fitting target."
)

fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent
  )
  fit <- fit_default_first(
    fit_model,
    retry_start = function() {
      cli::cli_abort("The 2021 report does not provide aligned annual N-at-age or F-at-age surfaces for a source-based retry.")
    },
    attempt_file = file.path(
      .assessment_cache(assessment_id, translation_tier, database = source),
      "initialization_attempts.rds"
    )
  )
}
