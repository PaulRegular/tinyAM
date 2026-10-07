translate_stock <- function(source) {
  years <- 1968:2016
  ages <- 1:10
  surveys <- c("NEFSC Albatross spring trawl", "NEFSC Bigelow spring trawl")

  inputs <- source$inputs
  keep_index <- inputs$type == "index" &
    inputs$measure == "numbers_at_age" &
    inputs$survey %in% surveys &
    ((inputs$survey == surveys[[1L]] & inputs$age >= 3 & inputs$age <= 10) |
       (inputs$survey == surveys[[2L]] & inputs$age >= 3 & inputs$age <= 7))
  inputs <- inputs[inputs$type != "index" | keep_index, , drop = FALSE]

  obs <- database_to_tam_obs(
    source$assessment$assessment_id[[1L]],
    inputs,
    years = years,
    ages = ages,
    surveys = surveys,
    sampling_times = c(
      "NEFSC Albatross spring trawl" = 95.7 / 365,
      "NEFSC Bigelow spring trawl" = 95.7 / 365
    ),
    assumptions = source$assumptions
  )

  obs$catch$age_factor <- factor(obs$catch$age, levels = ages)

  list(
    years = years,
    ages = ages,
    age_plus_group = 10,
    obs = obs,
    settings = list(
      N_settings = list(process = "iid", init = "exp"),
      F_settings = list(process = "ar1", mu_form = ~ 0 + age_factor, mean_ages = 6:10),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 1,
                            fill_missing = FALSE)
    ),
    comparison_scales = c(
      N = 1e-6, recruitment = 1e-6, ssb = 1e-3, biomass = 1e-3,
      F = 1, M = 1, F_bar = 1
    ),
    background = c(
      "### Atlantic mackerel: accepted 2018 ASAP assessment",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | Final ASAP Run 118 covers 1968-2016. | Fit 1968-2016. | The full detailed time series is represented. |",
      "| Ages | Combined sexes, ages 1-10+, with recruitment at age 1. | Ages 1-10, with age 10 as the plus group. | The reported age structure is preserved. |",
      "| N | ASAP estimates numbers at age and annual age-1 recruitment. | An IID abundance process with exponential initial abundance. | This is a simpler state process than the accepted ASAP model. |",
      "| F | One fishery with time-constant selectivity, fixed at full selection for ages 6-10+. | Age-specific mean F with AR1 deviations; mean F is summarized over ages 6-10. | tinyAM does not represent the accepted fixed selectivity curve directly. |",
      "| M | Fixed at 0.2 per year for all ages and years. | Use the fixed 0.2 schedule. | This part of the accepted model is retained directly. |",
      "| Catch | The accepted model uses combined U.S.-Canadian catch-at-age, plus its own catch and age-composition likelihoods. | Use published combined catch-at-age once with one lognormal observation error. | tinyAM does not reproduce ASAP's catch likelihood or catch-composition formulation; zero cells are treated as missing by the log-scale observation model. |",
      "| Index | The accepted model uses a range-wide egg SSB index and NEFSC spring trawl indices at ages 3+: Albatross ages 3-10 and Bigelow ages 3-7. | Use the two published age-specific spring trawl series and separate survey-age catchability; a constant timing of 95.7/365 is used. | The intermittent aggregate egg SSB index cannot be represented as an age-specific index without inventing age composition, so it is omitted from the fit. Timing is based on a reported time-series mean, not year-specific dates. |",
      "| Weights and maturity | Annual combined U.S.-Canadian catch/SSB weights and annual Canadian maturity ogives are used. | Apply both annual source surfaces directly. | The maturity observations represent the northern spawning contingent in the combined-stock assessment. |",
      "| SSB | The accepted ASAP model reports annual SSB from its full population and biological structure. | Calculate SSB from tinyAM's fitted abundance, translated weights, and maturity. | Differences in process equations and the omitted aggregate egg index remain. |",
      "",
      "The database also retains the intermittent combined egg-survey SSB estimates and the published January-1 and exploitable biomass series. The 2025 management-track assessment informs current advice but does not provide the full age-specific tables available in the accepted 2018 detailed assessment."
    )
  )
}