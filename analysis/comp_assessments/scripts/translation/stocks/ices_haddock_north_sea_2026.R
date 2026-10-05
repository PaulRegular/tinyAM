translate_stock <- function(source) {
  years <- 1972:2026
  ages <- 0:8
  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    source$inputs,
    years = years,
    ages = ages,
    assumptions = source$assumptions
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 8,
    obs = obs,
    comparison_scales = c(
      ssb = 1e-3, recruitment = 1e-3, N = 1e-3, abundance = 1e-3,
      biomass = 1e-3, biomass_at_age = 1e-3
    ),
    settings = list(
      N_settings = list(process = "iid", init = "free"),
      F_settings = list(process = "ar1", mu_form = NULL),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = TRUE),
      index_settings = list(q_form = ~ 0 + q_key,
                            sd_form = ~ 0 + survey,
                            sd_supplied = ~ relative_sd,
                            fill_missing = TRUE)
    ),
    background = c(
      "### Northern Shelf haddock: accepted 2026 SAM assessment",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | Population and survey series span 1972–2026; catch ends in 2025. | Fit years 1972–2026. | Retains the final survey year; catch is unavailable in 2026. |",
      "| Ages | Ages 0–8+, with age 8 as a plus group. | Ages 0–8, with age 8 as the plus group. | The age range and plus group are retained. |",
      "| N | Random-walk recruitment; N-process variance differs for recruits, ages 1–7 and the plus group; the first state has no process density. | Random-walk recruitment, IID abundance process and free initial abundance. | tinyAM cannot match the source's exact variance sharing. |",
      "| F | Time-varying F with age-correlated process deviations. | AR1 F process over model ages and years. | This approximates the source correlation structure; its temporal process is not identical. |",
      "| M | Fixed annual age-specific M from the accepted fit. | Fixed supplied M from the accepted run. | The source numerical surface is retained. |",
      "| Catch | One total-catch fleet, recorded as numbers-at-age. | Source numbers-at-age and one estimated catch-error scale. | tinyAM uses independent lognormal errors. |",
      "| Index | Q1 ages 1–8+ and Q3+Q4 ages 0–8+; q varies by survey and age; selected ages use a density-dependent q power. | Native age-specific survey indices, survey-by-age q, and supplied relative SD factors. | tinyAM retains the observation weights but has no q-power term or independent-error correlation structure. |",
      "| Weights and maturity | Annual stock weights and maturity; spawning fractions are zero. | Annual source values retained. | No biological surface is substituted. |",
      "",
      "Survey indices retain their native scale because the source documents do not provide a calibrated physical unit.",
      "",
      "Comparison units: tinyAM numbers, recruitment and biomass are converted to the accepted assessment's reported units in the dashboard and comparison tables."
    )
  )
}
