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

  # Use the accepted assessment's median reported Mbar as the fixed baseline.
  # This is a replication-oriented plug-in choice; tinyAM estimates deviations
  # around this level because estimating both the M level and process was unstable.
  obs$weight$M_assumption <- median(
    source$outputs$value[source$outputs$measure == "Mbar"],
    na.rm = TRUE
  )

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
  obs$index$smith_sound_year <- as.integer(
    obs$index$year %in% 1995:2007
  ) # Approximate offshore-availability effect during years with >10 kt in Smith Sound

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
        process = "rw", # Approximates xteNCAM's near-unit temporal F correlation; the closer 2D AR1 approximation did not converge reliably.
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
        fill_missing = FALSE
      ),
      index_settings = list(
        q_form = ~ q_key + smith_sound_year:smith_sound_q_key, # Allow adjustment of RV q to account for Smith Sound years
        sd_form = ~ 1,
        fill_missing = FALSE
      )
    ),
    background = c(
      "### Northern cod (2J3KL): accepted 2025 xteNCAM assessment",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | The accepted xteNCAM model covers 1954–2024; commercial catch-at-age is reported from 1962. | Fit 1968–2024, retaining the historical tinyAM Northern cod analysis window. | The shorter period is an analysis simplification rather than a feature of the accepted assessment. |",
      "| Ages | The accepted population model includes ages 0–14; commercial catch and the fall RV survey cover ages 2–14. | Model ages 2–14, with recruitment entering at age 2. | Ages 0–1 and their juvenile survey information are omitted, so recruitment and total abundance are not directly comparable with the accepted assessment. |",
      "| Recruitment | Recruitment enters the accepted model at age 0 and is informed by juvenile indices and a Beverton–Holt stock-recruitment relationship. | Recruitment enters at age 2 and follows tinyAM's temporal recruitment process. | tinyAM does not reproduce the age-0 recruitment process, juvenile indices, or stock-recruitment relationship used by xteNCAM. |",
      "| N | Population abundance in xteNCAM follows cohort survival within a broader state-space model. | Cohort abundance after recruitment follows deterministic survival, with first-year abundance initialized using exponential survivorship. | No additional N-process deviations are estimated in this simplified translation. The exponential initialization is closer to the accepted model's initial-age structure than freely estimating all initial abundances. |",
      "| F | Fishing mortality varies across ages and years with strongly correlated process variation; the accepted assessment estimates very high temporal correlation and substantial age correlation. | Fit independent age-specific temporal random walks in F and summarize mean F over ages 5–14. | The random walk approximates the strong temporal persistence in xteNCAM, but does not reproduce correlation among ages. A closer two-dimensional AR1 approximation did not converge reliably. |",
      "| M | Natural mortality varies through time around a baseline level and includes correlated age-year process variation plus a Capelin-to-cod biomass effect. | Use the median reported accepted-assessment Mbar as a fixed baseline and estimate an AR1 M process from 1984 onward, with neighbouring ages coupled into blocks; summarize mean M over ages 5–14. | Estimating both the overall M level and its process was unstable. The accepted Capelin effect and full age-year M structure cannot be reproduced, so the reported Mbar level is used as a replication-oriented plug-in calibration. |",
      "| Catch | xteNCAM uses commercial catch-age composition together with reported landings treated as bounded information on total removals. | Fit the reported catch-at-age numbers directly using tinyAM's lognormal observation model and omit the bounded landings component. | The catch likelihood and treatment of total removals differ substantially between models. Zero catch cells are treated as missing by tinyAM's log-scale observation model. |",
      "| Index | The fall RV survey covers ages 2–14, with separate catchability for ages 2–5 and common catchability for ages 6–14. xteNCAM also allows offshore RV availability to vary during the Smith Sound aggregation period. | Use the fall RV survey only, with separate q for ages 2–5 and shared q for ages 6–14. Add age-block-specific q adjustments during 1995–2007, when more than 10 kt of cod were estimated in Smith Sound. Survey timing is approximated as 0.75. | The Smith Sound term is a simplified proxy for xteNCAM's age- and year-varying offshore-availability process. Sentinel, Smith Sound, juvenile, and other survey information used by the accepted assessment is otherwise omitted. |",
      "| Weights | xteNCAM uses annual age-specific biological weights, including beginning-of-year stock weights and separate catch weights. | Use the reported annual stock weight-at-age values for population biomass calculations. | tinyAM does not reproduce all source-specific uses of separate weight series within the accepted likelihood. |",
      "| Maturity | Female maturity-at-age varies by calendar year; the underlying maturity model includes a cohort effect. | Use the reported calendar-year-by-age female maturity values directly. | The reported table already contains the resulting annual maturity-at-age values, so no cohort-year shift is applied. tinyAM has no explicit sex structure. |",
      "| SSB | The accepted assessment calculates spawning biomass from the age-0+ population using its full mortality, maturity, and population structure. | Calculate SSB from modeled ages 2–14 using the reported weights and maturity-at-age. | Ages 0–1 and several source-model processes are absent, so SSB is an approximate rather than exact reproduction of the accepted quantity. |",
      "",
      "The source report shows age-specific population and mortality outputs graphically but does not provide the complete numerical surfaces needed for age-by-age comparison. The tinyAM model is therefore an illustrative abstraction of the accepted assessment rather than an exact reproduction. Comparisons are limited to numerical aggregate outputs available in the assessment database and should be interpreted in light of the structural differences described above."
    )
  )
}
