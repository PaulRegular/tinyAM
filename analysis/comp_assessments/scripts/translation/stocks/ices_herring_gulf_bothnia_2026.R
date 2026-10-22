## Observations ----
years <- 1980:2025
ages <- 1:15
survey <- "BIAS acoustic survey"

model_inputs <- source$inputs[
  !(source$inputs$type == "catch" & source$inputs$measure == "numbers_at_age"),
]
stock_weights <- model_inputs[
  model_inputs$type == "catch_weight" & model_inputs$measure == "weight_at_age",
]
stock_weights$type <- "weight"
stock_weights$notes <- paste(
  stock_weights$notes,
  "Used as a proxy stock-weight surface because no separate stock-weight series was recovered."
)
model_inputs <- rbind(model_inputs, stock_weights)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  model_inputs,
  years = years,
  ages = ages,
  surveys = survey,
  sampling_times = c("BIAS acoustic survey" = 0.75),
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0

settings <- list(
  N_settings = list(process = "iid", init = "exp"),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "ar1", mu_form = ~logistic(age), mean_ages = 3:7),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~1),
  index_settings = list(q_form = ~0 + q_block, sd_form = ~1)
)
prepared <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = years, ages = ages), settings
))

## Model and comparison ----
age_plus_group <- 15L
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Published ICES aggregate SSB compared with tinyAM SSB over 1980-2025.",
    reason = "tinyAM begins at age 1, uses catch weight as a proxy for stock weight, and applies start-of-year spawning biomass because the exact source spawning-time settings and age-specific official biomass outputs were not recovered."
  ),
  recruitment = list(
    status = "non_equivalent",
    definition = "ICES recruitment is age 0; the tinyAM model begins at age 1.",
    reason = "Catch compositions begin at age 1 and official N-at-age is unavailable to align recruitment definitions."
  ),
  M = list(
    status = "non_equivalent",
    definition = "The fixed 2021 benchmark M input is reused in tinyAM; this is not an estimated-output comparison.",
    reason = "Any numerical match is by construction because tinyAM uses the same fixed source assumption."
  ),
  F_bar = list(
    status = "non_equivalent",
    definition = "ICES Fbar is the unweighted mean over ages 3-7; tinyAM Fbar is abundance-weighted over ages 3-7.",
    reason = "The report does not provide age-specific official F and N needed for a like-for-like comparison."
  ),
  M_bar = list(
    status = "non_equivalent",
    definition = "The ICES assessment uses fixed M from the 2021 benchmark; tinyAM applies the same fixed input.",
    reason = "Input consistency is not an independent comparison of estimated M."
  )
)

background <- c(
  "### Gulf of Bothnia herring: Official 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | Total catch 1963-2025; age compositions 1980-2025, ages 1-15+. | Fit 1980-2025 and ages 1-15+, with age 15 as the plus group. | Older catch-only years are not fitted because the age compositions are unavailable. |",
  "| Abundance and recruitment | SS3 hybrid-F model; age-0 recruitment and an exploited initial population. | Age-1 recruitment with an exponential initial abundance and IID cohort deviations. | The recruitment age and population-state structure differ. |",
  "| Fishing mortality | SS3 hybrid F; Fbar is the unweighted mean over ages 3-7. | AR1 F process with a logistic age mean; tinyAM reports abundance-weighted Fbar over ages 3-7. | The F process and Fbar definition differ. |",
  "| Natural mortality | Fixed age-varying M from the 2021 benchmark, held constant through time. | Use the benchmark M vector after linear interpolation to integer ages. | Age 15 is the tinyAM 15+ group and uses the source age-15 value. |",
  "| Catch and acoustic index | Total catch in tonnes with age compositions; BIAS acoustic index at ages 1-15+ from 2007. | Reconstruct catch numbers at age from total biomass, number proportions, and catch weights; fit the area-corrected BIAS age index with age-specific q. | tinyAM uses lognormal catch-at-age and index likelihoods; the SS3 composition weights and exact index precision were not recovered. |",
  "| Biology and SSB | Annual catch weights and maturity ogives are reported; source spawning settings are not available in the recovered outputs. | Use catch weights as the stock-weight proxy, annual maturity, and start-of-year SSB. | SSB is an approximate aggregate comparison. |",
  "| Outputs | Annual SSB, age-0 recruitment, and unweighted Fbar are available; age-specific N/F/M outputs were not recovered. | Compare aggregate SSB; show recruitment and Fbar as non-equivalent. | Total abundance and age-specific state comparisons remain unavailable. |",
  "",
  "The BIAS survey occurs in late September to early October; 0.75 is an approximate within-year time. Table 4.3.6 values are already area-corrected and are not multiplied by AREA_CORR_FACTOR again. The assessment also uses a 1990-2006 trapnet index, but its numeric input series and likelihood details were not recoverable. The fit uses package-default starting values."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = settings$N_settings,
    rec_settings = settings$rec_settings,
    F_settings = settings$F_settings,
    M_settings = settings$M_settings,
    catch_settings = settings$catch_settings,
    index_settings = settings$index_settings,
    silent = silent
  )
}
