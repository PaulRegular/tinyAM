## Observations ----
years <- 1995:2024
ages <- 2:11

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = "Sentinel mobile",
  sampling_times = c("Sentinel mobile" = 0.54),
  assumptions = source$assumptions
)

obs$pop_covariates$F_age_block <- cut(
  obs$pop_covariates$age,
  breaks = c(1, 3, 5, 7, 9, 11),
  labels = c("2-3", "4-5", "6-7", "8-9", "10-11"),
  include.lowest = TRUE
)
obs$pop_covariates$F_period <- cut(
  obs$pop_covariates$year,
  breaks = c(1994, 1999, 2005, 2013, 2021, 2024),
  labels = c("1995-1999", "2000-2005", "2006-2013", "2014-2021", "2022-2024")
)

comparison_outputs <- source$outputs[
  source$outputs$year %in% years & source$outputs$measure != "Fbar" &
    (source$outputs$measure != "total_numbers" | source$outputs$age_group == "2+"),
  , drop = FALSE
]
comparison_aggregates <- "ssb"
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Annual SSB, ages 2-11+",
    reason = "The population fit uses only the Sentinel mobile index and simplifies the source likelihoods and state structure."
  ),
  recruitment = list(
    status = "approximate",
    definition = "Age-2 abundance in the same calendar year",
    reason = "The source and tinyAM both report age-2 abundance, but use different recruitment and population-state formulations."
  )
)
age_plus_group <- 11

## Background and model ----
background <- c(
  "### Northern Gulf cod (3Pn4RS): 2025 DFO assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | 1973-2024; ages 2-11+. | 1995-2024; ages 2-11+. | Restrict to the years with the translated Sentinel mobile age index. |",
  "| N and recruitment | First-year abundance is estimated freely; age-2 recruitment has two mean periods and stationary AR(1) variation. | Free initial abundance in 1995; AR(1) recruitment with one mean for the post-1990 period. | The initial-state integration and pre-1995 history are omitted. |",
  "| F | Age/year-correlated deviations around 16 age- and year-block means. | AR(1) age/year process around additive age-block and period means. | TinyAM uses a simpler mean surface and correlation structure. |",
  "| M | Fixed age means plus AR(1) age/year deviations, grouped at ages 4-5, 6-7, 8-9 and 10-11+. | The documented baseline means with an AR(1) process; age 4 and age 5 have separate process blocks, then 6-7, 8-9 and 10-11+. | TinyAM requires supplied means to be constant within each process block; the source groups age 4 and 5 despite different means. The process also starts in 1995 rather than 1984. |",
  "| Catch | Catch-at-age and censored catch observations are modeled with the Baranov catch equation. | Published catch-at-age fitted with lognormal observation error. | Censoring and the source catch likelihood are not reproduced. |",
  "| Index | Six age-specific survey series with survey-specific timing and observation models. | Sentinel mobile only, with a separate q by age and approximate mid-July timing (0.54). | Five official index series and their observation weights are omitted; the retained index is relative mean fish per standardized tow. |",
  "| Biology | Annual stock weights and female maturity are from separate upstream models. | Digitized annual fitted curves from Figures 90 and 70. | Values are approximate readings from published plots; empirical points are excluded. |",
  "| Comparison | Published annual abundance, biomass, F, M, recruitment and SSB. | Compare annual SSB, recruitment, N, F and M over 1995-2024. | Source estimates are rounded and the model definitions differ. |",
  "",
  "The fit uses package-default starting values. No official population estimates are used as inputs or starting values."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "ar1", mu_form = ~1),
    F_settings = list(
      process = "ar1",
      mu_form = ~0 + F_age_block + F_period
    ),
    M_settings = list(
      process = "ar1",
      mu_form = NULL,
      mu_supplied = ~M_assumption,
      age_breaks = c(4, 5, 6, 8, 10, 11),
      first_dev_year = 1995
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~1),
    silent = silent
  )
}
