translate_stock <- function(source) {
  years <- 1946:2026
  ages <- 3:12
  source_ages <- 3:15
  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    source$inputs,
    years = years,
    ages = source_ages,
    assumptions = source$assumptions
  )
  obs$index$q_key <- interaction(
    obs$index$survey,
    pmin(obs$index$age, 11),
    drop = TRUE,
    lex.order = TRUE
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 12,
    obs = obs,
    comparison_scales = c(
      ssb = 1e-3, recruitment = 1e-3, N = 1e-3, abundance = 1e-3,
      biomass = 1e-3, biomass_at_age = 1e-3
    ),
    settings = list(
      N_settings = list(process = "iid", init = "free"),
      F_settings = list(process = "rw", mu_form = NULL),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = TRUE),
      index_settings = list(q_form = ~ 0 + q_key,
                            sd_form = ~ 0 + survey,
                            fill_missing = TRUE)
    ),
    background = c(
      "### Northeast Arctic cod: accepted 2026 SAM assessment",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|---|---|---|",
      "| Years | Population and survey series span 1946–2026; catch ends in 2025. | Fit years 1946–2026. | Retains the final survey year; the 2026 catch row is missing. |",
      "| Ages | Population model ages 3–15, with survey observations reported as ages 3–12+. | Fit ages 3–12 with source biology retained through age 15; tinyAM's age-12 plus group uses hidden ages 12–15. | The grouped survey values remain 12+ observations and are compared with the corresponding grouped predictions. |",
      "| N | Stochastic recruitment, age-structured process variation, and no first-state process density. | Random-walk recruitment, IID abundance process, and free initial abundance. | These are the closest available settings; process variance sharing differs from SAM. |",
      "| F | Log-F random-walk innovations, with separate age-3 behavior and shared older states. | Age-specific log-F random walk. | The process family is retained, but SAM's exact state and variance sharing are not represented. |",
      "| M | Fixed annual age surface includes background mortality and externally iterated cannibalism. | Fixed supplied M retained through age 15. | The source age-specific mortality values contribute to tinyAM's hidden age-12+ group. |",
      "| Catch | Residual-catch numbers-at-age with missing observations retained. | Source values are kept; missing rows are estimated by tinyAM. | tinyAM uses independent lognormal catch errors. |",
      "| Index | Five age-specific series; ages 11–12 share catchability within each series; errors are correlated across age in SAM. | Time-invariant survey-by-age q, with ages 11–12 shared, and independent survey-specific lognormal errors. | tinyAM has no matching cross-age observation correlation. |",
      "| Weights and maturity | Original annual age-specific inputs through age 15; spawning fractions are zero. | Original annual values retained through age 15. | tinyAM uses hidden-age biology to form the age-12+ group. |",
      "",
      "Source: accepted 2026 JRN-AFWG assessment and its saved SAM fit. SAM's short-term recruitment forecast is not treated as a fitted recruitment estimate.",
      "",
      "Comparison units: tinyAM numbers, recruitment and biomass are converted to the accepted assessment's reported units in the dashboard and comparison tables."
    )
  )
}
