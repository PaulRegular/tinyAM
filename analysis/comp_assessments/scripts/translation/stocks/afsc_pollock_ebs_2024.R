translate_stock <- function(source) {
  years <- 1964:2024
  ages <- 1:15
  surveys <- c(
    "NMFS bottom-trawl VAST",
    "NMFS acoustic-trawl",
    "NMFS acoustic-trawl age-1 index"
  )
  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    source$inputs,
    years = years,
    ages = ages,
    weight_survey = "",
    index_weight_source = "source",
    maturity_reference_year = 1964,
    maturity_multiplier = 0.5,
    surveys = surveys,
    assumptions = source$assumptions
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 15,
    obs = obs,
    comparison_scales = c(
      N = 1e-9, recruitment = 1e-6, ssb = 1e-6, biomass_at_age = 1e-6
    ),
    settings = list(
      N_settings = list(process = "off", init = "exp"),
      F_settings = list(process = "rw", mu_form = NULL),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 0 + q_key,
                            sd_form = ~ 0 + survey,
                            fill_missing = FALSE)
    ),
    background = c(
      "## Eastern Bering Sea pollock: accepted 2024 Model 23.0",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | Population model years 1964-2024. | Fit 1964-2024. | The accepted model period is retained. |",
      "| Ages | Ages 1-15, with age 15 as the model plus group. | Ages 1-15, with age 15 as the plus group. | The model age range is retained; report tables may group ages 10-15 as 10+. |",
      "| N | Age-1 recruitment varies by year; initial ages 2-15 have a shared log mean and regularized age deviations. | Exponential initial abundance, deterministic older-age survival, and tinyAM's recruitment process. | The source initial-age penalty and recruitment variance are not reproduced. |",
      "| F | Fishing mortality varies annually and uses age-selectivity curves with time-varying parameters. | An age-year random-walk process for F. | This is a broad process approximation, not the source selectivity model. |",
      "| M | Fixed age-specific values: 0.9 at age 1, 0.45 at age 2, and 0.3 at ages 3-15. | The same supplied age-specific M; no M process is estimated. | The accepted numerical vector is retained. |",
      "| Catch | Total fishery biomass and annual age compositions; the 2024 fishery composition is unavailable. | Reconstruct numbers-at-age using source catch weights; fit only observed age-year cells. | tinyAM fits age-specific observations rather than separate total-biomass and composition likelihoods. |",
      "| Index | VAST and acoustic-trawl biomass indices with age compositions; ATS age 1 is a separate index. | Reconstruct age-specific indices with source survey weights; use time-invariant q by survey and age, and fit the age-1 series separately. | tinyAM uses independent lognormal errors and does not reproduce the VAST covariance or composition likelihoods. |",
      "| Weights and maturity | Annual stock, catch, and survey weights; a fixed maturity vector and female fraction 0.5. | Stock weights for biology, source weights for observation conversion, and the source maturity vector multiplied by 0.5. | The accepted biological inputs and female spawning convention are retained. |",
      "",
      "Historical fishery CPUE and acoustic-vessel-only indices are omitted because the database has no age compositions for translating them into age-specific indices. The source also applies observation covariance and composition weights that tinyAM does not represent.",
      "",
      "The extracted accepted outputs do not include F-at-age values, so the accepted F surface is unavailable in the comparison dashboard. A free-initial-age fit failed under both random-walk and AR1 F processes; the exponential N0 fit converged with the random-walk F process."
    )
  )
}
