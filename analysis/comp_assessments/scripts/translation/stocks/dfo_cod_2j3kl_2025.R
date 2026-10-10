## Observations ----
if (!exists("smith_sound", inherits = FALSE)) smith_sound <- TRUE
if (!exists("juveniles", inherits = FALSE)) juveniles <- FALSE
years <- 1968:2024
ages <- if (juveniles) 0:14 else 2:14
surveys <- "DFO fall RV survey"
inputs <- source$inputs
if (smith_sound) {
  surveys <- c(surveys, "Smith Sound acoustic survey")
  # Ages 0, 14, 15 and 16 have zero sample counts in every reported year.
  # Remove them before biomass conversion: no weights or q can be identified.
  smith_sample <- inputs$survey %in% "Smith Sound acoustic survey" &
    inputs$unit %in% "fish_sampled"
  stopifnot(all(inputs$value[smith_sample & inputs$age %in% c(0, 14:16)] == 0))
  inputs <- inputs[!(smith_sample & inputs$age %in% c(0, 14:16)), ]
}
if (juveniles) {
  surveys <- c(surveys, "Fleming juvenile survey", "Newman Sound juvenile survey")
}

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("Fleming juvenile survey" = 9.5 / 12,
                     "Newman Sound juvenile survey" = 9 / 12),
  assumptions = source$assumptions
)

# Use the accepted assessment's median reported Mbar as the fixed baseline.
# This is a replication-oriented plug-in choice; tinyAM estimates deviations
# around this level because estimating both the M level and process was unstable.
obs$weight$M_assumption <- median(
  source$outputs$value[source$outputs$measure == "Mbar"],
  na.rm = TRUE
)

rv <- obs$index$survey == "DFO fall RV survey"
juvenile <- grepl("juvenile", obs$index$survey)
obs$index$q_key <- ifelse(
  obs$index$age <= 5,
  paste("RV age", obs$index$age),
  "RV age 6+"
)
obs$index$q_key[!rv] <- paste("Smith age", obs$index$age[!rv])
obs$index$q_key[juvenile] <- paste("juvenile age", obs$index$age[juvenile])
obs$index$q_key <- factor(obs$index$q_key)
obs$index$sd_key <- factor(ifelse(juvenile, "juvenile", obs$index$survey))
obs$index$smith_sound_q_key <- cut_ages(
  pmax(obs$index$age, 2), c(seq(2, 12, 2), 14)
)
levels(obs$index$smith_sound_q_key) <- paste("age", levels(obs$index$smith_sound_q_key))
# The interaction must be zero, rather than NA, on non-RV observations.
obs$index$smith_sound_q_key[!rv] <- levels(obs$index$smith_sound_q_key)[1L]
obs$index$smith_sound_year <- as.integer(
  rv & obs$index$year %in% 1995:2007
) # Approximate offshore-availability effect during years with >10 kt in Smith Sound


## Background and comparisons ----

comparison_scales <- c(N = 1e-06, recruitment = 1e-06, ssb = 1e-06, abundance = 1e-06, biomass = 1e-06,
    F_bar = 1, M_bar = 1, Z_bar = 1)

background <- c("### Northern cod (2J3KL): accepted 2025 xteNCAM assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Years | The accepted xteNCAM model covers 1954–2024; commercial catch-at-age is reported from 1962. | Fit 1968–2024, retaining the historical tinyAM Northern cod analysis window. | The shorter period is an analysis simplification rather than a feature of the accepted assessment. |",
    "| Ages | The accepted population model includes ages 0–14; commercial catch and the fall RV survey cover ages 2–14. | Model ages 2–14, with recruitment entering at age 2. | Comparisons use accepted N at age 2 and totals over ages 2–14; native age-0 recruitment and age-0+ totals have different definitions. |",
    "| Recruitment | Recruitment enters the accepted model at age 0 and is informed by juvenile indices and a Beverton–Holt relationship with same-year SSB. | Age-2 recruitment follows a random walk. | BH is available in tinyAM, but the retained model does not equate age-2 survivors with the source age-0 relationship. |",
    "| N | Population abundance in xteNCAM follows cohort survival within a broader state-space model. | Cohort abundance after recruitment follows deterministic survival, with first-year abundance initialized using exponential survivorship. | No additional N-process deviations are estimated in this simplified translation. The exponential initialization is closer to the accepted model's initial-age structure than freely estimating all initial abundances. |",
    "| F | Fishing mortality varies across ages and years with strongly correlated process variation; the accepted assessment estimates very high temporal correlation and substantial age correlation. | Fit independent age-specific temporal random walks in F and summarize mean F over ages 5–14. | The random walk approximates the strong temporal persistence in xteNCAM, but does not reproduce correlation among ages. A closer two-dimensional AR1 approximation did not converge reliably. |",
    "| M | Natural mortality varies through time around a baseline level and includes correlated age-year process variation plus a Capelin-to-cod biomass effect. | Use the median reported accepted-assessment Mbar as a fixed baseline and estimate an AR1 M process from 1984 onward, with neighbouring ages coupled into blocks; summarize mean M over ages 5–14. | Estimating both the overall M level and its process was unstable. The accepted Capelin effect and full age-year M structure cannot be reproduced, so the reported Mbar level is used as a replication-oriented plug-in calibration. |",
    "| Catch | xteNCAM uses commercial catch-age composition together with reported landings treated as bounded information on total removals. | Fit the reported catch-at-age numbers directly using tinyAM's lognormal observation model and omit the bounded landings component. | The catch likelihood and treatment of total removals differ substantially between models. Zero catch cells are treated as missing by tinyAM's log-scale observation model. |",
    "| Index | RV catchability is separate at ages 2–5 and shared at ages 6–14. Smith Sound biomass and sampled counts inform an inshore population component and offshore availability. | Retain RV q sharing and its 1995–2007 age-block availability adjustment. Include reconstructed Smith Sound numbers in years with both biomass and samples, with independent q by age and separate SD. RV timing is approximated as 0.75; Smith timing uses reported months. | The availability adjustment and reconstructed index approximate rather than reproduce the local-population likelihood. Sentinel and juvenile indices are omitted from the main age-2+ fit; juvenile integration is tested separately. |",
    "| Weights | xteNCAM uses annual age-specific biological weights, including beginning-of-year stock weights and separate catch weights. | Use the reported annual stock weight-at-age values for population biomass calculations. | tinyAM does not reproduce all source-specific uses of separate weight series within the accepted likelihood. |",
    "| Maturity | Female maturity-at-age varies by calendar year; the underlying maturity model includes a cohort effect. | Use the reported calendar-year-by-age female maturity values directly. | The reported table already contains the resulting annual maturity-at-age values, so no cohort-year shift is applied. tinyAM has no explicit sex structure. |",
    "| SSB | The accepted assessment calculates spawning biomass from the age-0+ population using its full mortality, maturity, and population structure. | Calculate SSB from modeled ages 2–14 using the reported weights and maturity-at-age. | Ages 0–1 and several source-model processes are absent, so SSB is an approximate rather than exact reproduction of the accepted quantity. |",
    "", "Tables 19–24 provide numerical age-specific N, biomass, mature biomass, Z, M and F. Native age-0 recruitment is not compared with age-2 recruitment. Common-age abundance and biomass comparisons use only the ages represented in both models; rounded age-specific tables have no reported uncertainty.",
    "", if (smith_sound) "Smith Sound numbers-at-age are reconstructed as B * p / sum(p * w), using sample proportions p and beginning-of-year stock weights w (Table 9). The framework biomass equation uses stock weights (2025/034, equation 2.14). These weights approximate survey-time body mass; survey-specific growth/condition is not recoverable. Only the seven years with both samples and biomass are used. Ages 0 and 14–16 contain only zero samples. Separate Smith Sound q is fitted at each informative age; it includes local availability and does not reproduce xteNCAM's separately modeled inshore component or its q normalization. These reconstructed values are translation products, not published abundance estimates." else "Smith Sound reconstruction is omitted in this sensitivity.",
    if (juveniles) "This sensitivity extends ages to 0–14 and includes Fleming and Newman indices with shared q by age, as in the accepted model. Timing is approximated as 9.5/12 for Fleming (September–October) and 9/12 for Newman (July–November); these are season midpoints, not recovered model timing parameters. tinyAM estimates F at ages 0–1, whereas xteNCAM fixes it to zero; that unsupported constraint is a material limitation of this trial." else "Juvenile surveys are retained in the database but excluded from this age-2+ fit. An age-0+ integration trial is recorded separately.")

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
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 5:14),
    M_settings = list(process = "ar1", mu_form = NULL, mu_supplied = ~M_assumption, first_dev_year = 1984,
        age_breaks = c(3, 5, 7, 9, 11, 14), mean_ages = 5:14),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~q_key + smith_sound_year:smith_sound_q_key, sd_form = if (length(surveys) >
        1L) ~0 + sd_key else ~1, fill_missing = FALSE),
    silent = silent
  )
}
