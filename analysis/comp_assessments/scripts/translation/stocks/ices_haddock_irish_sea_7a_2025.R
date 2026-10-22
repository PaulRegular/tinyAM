source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

years <- 1993:2024
ages <- 0:5

analysis_inputs <- source$inputs
missing_index <- analysis_inputs$type == "index" &
  is.na(analysis_inputs$value)
analysis_inputs <- analysis_inputs[!missing_index, , drop = FALSE]

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  analysis_inputs,
  years = years,
  ages = ages,
  surveys = c("NIGFS-WIBTS-Q1", "NIGFS-WIBTS-Q4", "NI-MIK", "UK-FSPW"),
  assumptions = source$assumptions
)
obs$index$q_key <- interaction(obs$index$survey, obs$index$age,
  drop = TRUE, lex.order = TRUE
)

age_plus_group <- 5L
comparison_outputs <- source$outputs[
  source$outputs$assessment_id == source$assessment$assessment_id[[1]] &
    source$outputs$year >= min(years) & source$outputs$year <= max(years), , drop = FALSE
]
comparison_aggregates <- c("ssb", "recruitment", "F_bar")
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Annual SSB as reported by ICES, compared with tinyAM mature biomass.",
    reason = "The report does not provide all spawning-time details needed to confirm that the tinyAM and ASAP biomass calculations use identical timing."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-0 recruitment in the same year, in thousands of fish.",
    reason = "The age and unit match; tinyAM's process and initial-state treatment are simplified."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "Official Fbar over ages 2-4; tinyAM Fbar is N-weighted.",
    reason = "Official N-at-age is unavailable, so the official Fbar cannot be re-expressed using tinyAM's population weighting."
  ),
  M = list(
    status = "approximate",
    definition = "Fixed M-at-age input from the 2016 Lorenzen table, ages 0-5.",
    reason = "The exact M vector in the 2025 ASAP fit was not recovered; equality in the dashboard checks input mapping only."
  )
)

background <- c(
  "### Irish Sea haddock: 2025 ICES ASAP assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1993-2024, ages 0-5. | Fit 1993-2024 and ages 0-5, with age 5 as the plus group. | The source grids are retained. |",
  "| N and recruitment | Single-fleet ASAP; recruitment is at age 0. | IID N process with exponential initialization and a separate random-walk recruitment process. | Process variances and initial-state integration differ from ASAP. |",
  "| F and M | Fbar ages 2-4; M is based on Lorenzen (1996). | Correlated random-walk F, with mean ages 2-4; use the 2016 Lorenzen M vector as fixed input. | Exact final-fit M values and the full ASAP process settings were not recovered. |",
  "| Catch | Catch numbers-at-age; landings and discards are combined in the model catch. | Retain catch-at-age with one estimated lognormal observation SD. | Catch scale is inferred as thousands of fish from table values. |",
  "| Indices and q | NIGFS-WIBTS-Q1, NIGFS-WIBTS-Q4, NI-MIK and UK-FSPW with age ranges and survey-specific likelihood settings. | Keep all four series and estimate survey-age q and survey-specific observation SDs. | Native index multipliers are unstated; CV/ESS settings and the complete ASAP likelihood are simplified. |",
  "| Biology and SSB | Annual stock weights and female maturity are reported. | Use the reported age-specific tables to calculate tinyAM mature biomass. | Stock-weight units are inferred and spawning-time details are not fully documented. |",
  "| Reported comparisons | Aggregate recruitment, SSB and Fbar are available with 95% intervals. | Compare those aggregate series; N-at-age and F-at-age panels are unavailable. | The published sources do not provide those age-specific population outputs. |",
  "",
  "The report table includes missing survey cells and zeros. Missing cells are excluded from the analysis observation table. tinyAM treats zero observations as missing under its lognormal likelihood. Source values remain unchanged in the database.",
  "",
  "The default-start fit converged (optimizer code 0, maximum absolute gradient 0.00087, and successful positive-definite sdreport). Agreement remains poor: mean absolute percent differences are 178.9% for recruitment and 523.5% for approximate SSB; terminal differences are +258.5% and +842.3%, respectively. Recruitment trend correlation is 0.695. This is a runnable structural approximation, not a close replication. Official N-at-age and F-at-age are unavailable, so no source-aligned population-start retry was possible."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) {
    arguments <- list(
      data = obs,
      years = years,
      ages = ages,
      N_settings = list(process = "iid", init = "exp"),
      rec_settings = list(process = "rw", mu_form = NULL),
      F_settings = list(process = "cor_rw", mean_ages = 2:4),
      M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
      catch_settings = list(sd_form = ~1),
      index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey),
      silent = silent
    )
    if (!is.null(start_par)) arguments$start_par <- start_par
    do.call(tinyAM::fit_tam, arguments)
  }
  retry_start <- function() {
    stop("The official outputs have no N-at-age or F-at-age surfaces for aligned source starts.")
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      source$assessment$assessment_id[[1]], "initialization_attempts.rds"
    )
  )
}
