## Observations ----
years <- 1991:2025
ages <- 1:8

weight_proxy <- source$inputs[
  source$inputs$type == "catch_weight" &
    source$inputs$measure == "weight_at_age", , drop = FALSE
]
weight_proxy$type <- "weight"
weight_proxy$basis <- "kg_per_fish"
weight_proxy$source_type <- "analysis_approximation"
weight_proxy$source_reference <- paste(
  weight_proxy$source_reference,
  "used only as an analysis-side stock-weight proxy"
)
weight_proxy$transformation <- paste(
  "Catch weights copied to the analysis data as a proxy because the full",
  "SAM stock-weight surface is not published."
)
weight_proxy$notes <- paste(
  "Analysis-only proxy from catch weights; this is not recorded as official",
  "stock weight in the canonical database."
)
translation_inputs <- rbind(source$inputs, weight_proxy)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  surveys = c("Rockall pt1", "Modelled Q3"),
  sampling_times = c("Rockall pt1" = 0.725, "Modelled Q3" = 0.75),
  assumptions = source$assumptions
)

obs$index$q_age_group <- ifelse(obs$index$age >= 5, "5-8", as.character(obs$index$age))
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_age_group,
  drop = TRUE, sep = "."
)
obs$index$relative_sd <- 1
q3_cv <- source$inputs[
  source$inputs$type == "index" &
    source$inputs$survey == "Modelled Q3" &
    source$inputs$measure == "relative_standard_error", , drop = FALSE
]
q3_key <- paste(q3_cv$year, q3_cv$age, sep = "\034")
obs_q3 <- obs$index$survey == "Modelled Q3"
obs_key <- paste(obs$index$year[obs_q3], obs$index$age[obs_q3], sep = "\034")
cv_index <- match(obs_key, q3_key)
if (anyNA(cv_index)) stop("Modelled Q3 indices need matching published CVs.")
obs$index$relative_sd[obs_q3] <- sqrt(log1p(q3_cv$value[cv_index]^2))

comparison_outputs <- source$outputs[
  source$outputs$year >= min(years) & source$outputs$year <= max(years), , drop = FALSE
]
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Common-definition SSB from official N-at-age, translated maturity and the catch-weight proxy.",
    reason = "The source stock-weight surface is estimated by SAM but is not published as a complete numerical table; catch weights are used as an analysis-only proxy."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-1 recruitment in the same calendar year.",
    reason = "Recruitment age matches, but the tinyAM random-walk and initial-state treatment simplify the SAM model."
  )
)
age_plus_group <- 8

## Model ----
background <- c(
  "### Rockall haddock: official 2026 ICES SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1991-2025, ages 1-8+; age 8 is the plus group. | Fit 1991-2025 and ages 1-8+. | The reported historical assessment grid is retained. |",
  "| Recruitment and N | Age-1 recruitment random walk; N process SD is separate for age 1 and shared over ages 2-8. | Random-walk recruitment and IID cohort-survival deviations for ages 2-8 with one shared process SD; exponential initial abundance. | tinyAM treats recruitment separately and cannot retain SAM's age-1 N variance group or full initial-state integration. |",
  "| F | Annual F process with AR(1) correlation across ages; ages 7-8 share a state and process variance is shared. | Correlated random-walk F, with Fbar summarized over ages 2-5. | The temporal process and age-correlation structure are represented, but other SAM parameter and likelihood details are simplified. |",
  "| M | Fixed, time-invariant scaled-Lorenzen M-at-age. | Use the published fixed M-at-age values. | The numerical M profile is retained. |",
  "| Catch | Total catch numbers-at-age, landings fraction, and catch weights; catch observation variance is shared by age with additional uncertainty from 2011. | Use total catch numbers-at-age and one shared observation SD. | tinyAM does not reproduce the additional year-age catch uncertainty structure. |",
  "| Indices and q | Rockall pt1 (1991-1997) and Modelled Q3 (1999-2025); q is shared at ages 5-8, and Modelled Q3 uses age-specific CV weights. | Preserve both relative-scale indices and source CV precision for Modelled Q3; estimate survey-age-group q and survey-specific SD. | Density-dependent q powers, exact SAM observation settings and the age-specific CV likelihood treatment are simplified. |",
  "| Weights and maturity | Fixed maturity; SAM estimates a stock-weight surface using catch-weight observations and a correlated process. | Use source maturity and an analysis-only catch-weight proxy for stock weights. | The complete official stock-weight surface is not published; the proxy is identified in the dashboard and is not added to the source database. |",
  "| Comparisons | N-at-age and F-at-age are reported for 1991-2025; SSB, recruitment, Fbar ages 2-5 and estimated catch are reported in summary tables. | Compare age-specific N/F, fixed M, common-definition SSB, recruitment, total abundance and Fbar over 1991-2025. | Summary Low/High bounds are retained without calling them 95% intervals; age-specific N/F uncertainty is not available. |",
  "",
  "The exact 2026 survey fractions are not stated in the report. The Q3 timing values 0.725 and 0.75 were checked in the overlapping 2025 SAM object and carried forward as documented assumptions. Modelled Q3 CVs are transformed to relative log-SD factors using sqrt(log(1 + CV^2)); Rockall pt1 has no published CV and uses a neutral relative factor of 1. The fit uses package-default starting values."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp", sd_form = ~1),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + survey,
      sd_supplied = ~relative_sd
    ),
    silent = silent
  )
}
