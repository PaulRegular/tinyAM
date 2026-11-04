## Observations ----
years <- 1969:2025
ages <- 2:12
assessment_id <- source$assessment$assessment_id[[1L]]
surveys <- c("UK-CBT-late", "UK-COT", "Q1SWBeam-nonoffset", "FSP-UK")

index_rows <- source$inputs[
  source$inputs$type == "index" & source$inputs$measure == "numbers_at_age",
  , drop = FALSE
]
index_rows$year <- as.integer(index_rows$year)
index_rows$age <- as.integer(index_rows$age)
active <- with(index_rows,
  (survey == "UK-CBT-late" & year %in% 2003:2025 & age %in% 3:11) |
    (survey == "UK-COT" & year %in% 1988:2015 & age %in% 3:11) |
    (survey == "Q1SWBeam-nonoffset" & year %in% 2006:2025 & year != 2022 & age %in% 2:11) |
    (survey == "FSP-UK" & year %in% 2004:2025 & age %in% 2:11)
)
index_rows <- index_rows[active & is.finite(index_rows$value) & index_rows$value > 0, , drop = FALSE]
index_measure <- source$inputs$type == "index" & source$inputs$measure == "numbers_at_age"
inputs <- rbind(source$inputs[!index_measure, , drop = FALSE], index_rows)

sampling_times <- c(
  "UK-CBT-late" = 0.5,
  "UK-COT" = 0.5,
  "Q1SWBeam-nonoffset" = 0.175,
  "FSP-UK" = 0.725
)
obs <- database_to_tam_obs(
  assessment_id, inputs, years = years, ages = ages,
  surveys = surveys, sampling_times = sampling_times,
  assumptions = source$assumptions
)
obs$index$sd_block <- factor(obs$index$survey, levels = surveys)
obs$index$q_key <- interaction(
  obs$index$survey, pmin(obs$index$age, 7L), drop = TRUE, sep = "."
)

## Model ----
age_plus_group <- 12L
N_settings <- list(process = "iid", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "ar1", mu_form = ~factor(age), mean_ages = 3:9)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)

dat <- tinyAM::prepare_tam(
  data = obs, years = years, ages = ages,
  N_settings = N_settings, rec_settings = rec_settings,
  F_settings = F_settings, M_settings = M_settings,
  catch_settings = catch_settings, index_settings = index_settings
)

comparison_aggregates <- c("abundance", "biomass", "ssb", "F_bar")
comparison_definitions <- list(
  abundance = list(
    status = "matched",
    definition = "Total start-of-year abundance summed over ages 2-12+.",
    reason = "Official N-at-age covers the full fitted age range; total abundance is derived without uncertainty."
  ),
  biomass = list(
    status = "matched",
    definition = "Total biomass over ages 2-12+ using the source stock weights.",
    reason = "The same-age N and stock-weight surfaces reproduce the source total biomass to table rounding."
  ),
  ssb = list(
    status = "matched",
    definition = "SSB from start-of-year N, annual stock weight, and the source age-specific maturity.",
    reason = "The source N, stock-weight, and maturity tables reproduce reported SSB within 0.04%; small differences reflect report rounding."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-2 abundance in the same calendar year, in thousands of fish.",
    reason = "The recruitment age and units match; population process assumptions differ."
  ),
  F_bar = list(
    status = "matched",
    definition = "Abundance-weighted F over ages 3-9, reconstructed from the source N-at-age and F-at-age.",
    reason = "The database retains native unweighted ICES Fbar separately; the comparison uses a common population-weighted definition."
  ),
  M = list(
    status = "approximate",
    definition = "Natural mortality is compared with the fixed source input of 0.1 per year.",
    reason = "Both sides use the supplied fixed value; equality is an input check, not evidence that M was estimated by either model."
  ),
  M_bar = list(
    status = "approximate",
    definition = "Abundance-weighted Mbar over the model ages, compared with the fixed source input.",
    reason = "Both sides use the supplied fixed value; equality is an input check, not evidence that M was estimated by either model."
  )
)

background <- c(
  "### Western English Channel common sole: 2026 ICES XSA assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | XSA fitted through 2025, ages 2-12+, with age 12 as the plus group. | Fit the same years and ages. | 2026 aggregate recruitment and SSB are forecasts and are excluded from model fitting. |",
  "| Catch and population | XSA treats landings numbers-at-age as input to the survivor equations. | Fit the same catch-at-age values in tinyAM's catch observation model. | tinyAM estimates population states with a likelihood on catch, rather than treating catch as the XSA accounting input. |",
  "| N, F and M | XSA estimates N/F through its survivors analysis and fixes M at 0.1 per year. | Use IID cohort deviations, a random walk for age-2 recruitment, AR(1) F, and fixed M=0.1. | The process model and catch likelihood are structural approximations to XSA. |",
  "| Tuning series and q | Four lognormal series; ages 3-11, with q plateauing at age 7. XSA applies a 15-year tricubic time taper. | Retain all four series, separate q and observation SD by series, and share q from age 7 onward. | tinyAM does not apply XSA's temporal taper; zero index cells are omitted from its lognormal fit. |",
  "| Survey timing | Source timing windows are 0-1 for UK-CBT/UK-COT, 0.1-0.25 for Q1SWBeam, and 0.7-0.75 for FSP. | Use the window midpoints: 0.5, 0.175, and 0.725. | These are point-time approximations to source windows. |",
  "| Biology and outputs | Annual stock weights, age-specific maturity, fixed M, N/F-at-age, and aggregate recruitment, biomass, SSB and Fbar. | Use source biology and compare available N/F and aggregate series. | Age-specific output uncertainty is not published; native Fbar is unweighted, while the dashboard's common-definition Fbar is N-weighted. |",
  "",
  "The source tuning file stores effort followed by catch numbers-at-age. The canonical database retains both source fields and reconstructs each index as catch-at-age divided by effort. Index values outside the official fleet year/age ranges, the Q1 2022 series, zero-effort rows, and zero-valued index cells are not passed to tinyAM. These choices are documented approximations; source rows remain unchanged. The fit uses package-default initialization without source-based starts."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings,
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    silent = silent
  )
}
