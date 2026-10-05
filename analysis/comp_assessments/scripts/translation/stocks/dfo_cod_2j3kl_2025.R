translate_stock <- function(source) {
  years <- 1968:2024
  ages <- 2:14

  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    source$inputs,
    years = years,
    ages = ages,
    surveys = "DFO fall RV survey",
    assumptions = source$assumptions
  )
  obs$weight$M_assumption <- 0.312
  obs$index$q_key <- factor(ifelse(
    obs$index$age <= 5,
    paste("age", obs$index$age),
    "age 6+"
  ))
  smith_breaks <- c(
    seq(min(obs$index$age), 13, 2),
    max(obs$index$age) + 1
  )
  obs$index$smith_sound_q_key <- cut(
    obs$index$age,
    breaks = smith_breaks,
    labels = paste0(
      "age ",
      head(smith_breaks, -1),
      "-",
      tail(smith_breaks, -1) - 1
    ),
    right = FALSE
  )
  obs$index$smith_sound_year <- as.integer(obs$index$year %in% 1995:2005)

  list(
    years = years,
    ages = ages,
    obs = obs,
    comparison_scales = c(
      recruitment = 1e-6, ssb = 1e-6, abundance = 1e-6,
      biomass = 1e-6, F_bar = 1, M_bar = 1, Z_bar = 1
    ),
    settings = list(
      N_settings = list(
        process = "off",
        init = "exp"
      ),
      F_settings = list(
        process = "rw", # because xteNCAM forces the 2DAR1 process to be close to a random walk
        mu_form = NULL,
        mean_ages = 5:14
      ),
      M_settings = list(
        process = "ar1",
        mu_form = NULL,
        mu_supplied = ~ M_assumption,
        first_dev_year = 1984,
        age_breaks = c(3, 5, 7, 9, 11, 14), # Coupling ages to avoid convergence issues
        mean_ages = 5:14
      ),
      catch_settings = list(
        sd_form = ~ 1,
        fill_missing = TRUE
      ),
      index_settings = list(
        q_form = ~ mono(q_key) + smith_sound_year:smith_sound_q_key, # Allow adjustment of RV q to account for Smith Sound years
        sd_form = ~ 1,
        fill_missing = TRUE
      )
    ),
    background = c(
      "### Northern cod (2J3KL): accepted 2025 xteNCAM assessment",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | The accepted population model covers 1954–2024; catch-at-age is reported from 1962. | Retain the existing 1968–2024 tinyAM comparison window and use the reported year-age values directly. | This keeps the simplified fit period unchanged while correcting the maturity-year indexing. |",
      "| Ages | The population model begins at age 0; catch and fall RV observations cover ages 2–14. | Model ages 2–14, with recruitment entering at age 2. | Juvenile index streams at ages 0–1 have no documented sampling time for tinyAM. Age-2 recruitment is not directly comparable with the source's age-0 recruitment. |",
      "| N | The source estimates annual abundance with cohort and age/year process structure. | IID cohort-process deviations and freely estimated initial abundance. | tinyAM cannot reproduce the source covariance or initial-state integration. |",
      "| F | Annual F has correlated variation across ages and years, with correlation forced to be close to 1. | A random walk across years was applied with random deviations across age. | Process is similar, but correlation across ages is unacounted for. |",
      "| M | M has a baseline of 0.312 per year, a Capelin-to-cod biomass effect, and correlated process errors. | Fix M at the reported baseline of 0.312 per year. | The numerical Capelin covariate is not tabulated, so the annual source M pattern cannot be reconstructed. |",
      "| Catch | Catch-at-age counts cover ages 2–14; reported landings are a separate bounded input. | Fit the direct catch-at-age values with tinyAM's lognormal observation model; omit reported landings. | The likelihoods differ. Source zero cells are treated as missing by tinyAM's log-scale observation model. |",
      "| Index | The fall RV age-specific index covers ages 2–14 in surveyed years; other inputs include sentinel, Smith Sound, and juvenile surveys. | Use the fall RV series, with separate q for ages 2–5 and shared q for ages 6–14. Also impose an effect to allow the model to reduce RV survey catchability to account for the Smith Sound event. Set its seasonal timing to 0.75. | The remaining series have unresolved timing, units, or age mappings. The reported 0.75 is a seasonal approximation. |",
      "| Weights and maturity | Annual stock weights and female maturity-at-age are reported by calendar year through 2024. | Use the reported year-by-age values directly. | The cohort effect describes how maturity was estimated; Table 8 reports calendar-year-by-age values, so no year shift is needed. tinyAM has no sex-structured population. |",
      "| Maturity and SSB | The source reports age-0+ population and spawning biomass. | Calculate SSB over modeled ages 2–14. | Age-0 and age-1 abundance is outside this fit, so SSB and total-abundance comparisons are approximate. |",
      "",
      "The source report shows age-specific population and mortality outputs graphically but does not provide the numeric surfaces needed for an age-by-age comparison. The model comparison is therefore limited to the numerical aggregate outputs stored in the database and should be interpreted alongside these differences."
    )
  )
}
