translate_stock <- function(source) {
  assessment_id <- source$assessment$assessment_id[[1]]
  years <- 1982:2016
  ages <- 0:7
  survey_names <- c("NEFSC Bigelow spring trawl", "NEFSC Bigelow fall trawl")

  inputs <- source$inputs
  inputs <- inputs[!(inputs$type == "catch" &
                       inputs$measure == "numbers_at_age" &
                       inputs$fleet != "Published total"), , drop = FALSE]
  obs <- database_to_tam_obs(
    assessment_id, inputs, years = years, ages = ages,
    weight_survey = "", surveys = survey_names,
    assumptions = source$assumptions
  )

  settings <- list(
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~ factor(age), mean_ages = 4),
    M_settings = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ 1,
                          fill_missing = FALSE)
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 7,
    obs = obs,
    settings = settings,
    comparison_scales = c(
      N = 1e-3, recruitment = 1e-3, ssb = 1e-3, F = 1, M = 1,
      F_bar = 1
    ),
    background = c(
      "### Summer flounder: accepted 2018 ASAP assessment",
      "",
      print_sources(source$assessment),
      "",

      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | F2018_BASE_V2 covers 1982-2017. | Fit 1982-2016. | The reported maturity-at-age series ends in 2016, so no terminal-year value is invented. |",
      "| Ages | Combined sexes, true ages 0-7+, with recruitment at age 0. | Ages 0-7, age 7 as the plus group. | The reported source age grouping is preserved. |",
      "| N | ASAP estimates abundance at age and recruitment, with its own process and initial-state treatment. | Exponential initial abundance with an IID abundance process. | This is a simpler process structure than ASAP. |",
      "| F | Four catch fleets have separate selectivity and time-varying F; fully recruited F is at age 4. | One aggregate fishery with an AR1 F process around estimated age-specific means; report Fbar at age 4. | tinyAM does not reproduce four fleet-specific selectivity or the ASAP time blocks. |",
      "| M | Fixed M-at-age for ages 0-7+; 0.26, 0.26, 0.26, 0.25, 0.25, 0.25, 0.25, 0.24 per year. | Use the published fixed age schedule. | The numerical age-specific input is preserved. |",
      "| Catch | The accepted run uses four fleet catches and age compositions. The published total catch-at-age is also reported. | Use the published total numbers-at-age once, with one observation SD. | The fleet and aggregate age tables do not match in every cell; using the published aggregate avoids adding both representations. |",
      "| Index | Many fishery-independent series are included; Bigelow spring and fall have age-specific indices. | Use only the available Bigelow spring/fall age indices with separate survey-age q values. | Other accepted indices are not currently in the database; timing uses seasonal midpoints 0.25 and 0.75. |",
      "| Weights and maturity | Table A90 gives 2013-2017 mean November SSB weights; Table A86 gives a moving-window maturity ogive through 2016. | Apply the published SSB weights as a static schedule and use the reported maturity years. | Static weights approximate an unavailable annual SSB-weight surface; fitting ends where maturity observations end. |",
      "| SSB | The accepted run reports annual SSB from its ASAP population and biological inputs. | Calculate SSB with tinyAM's annual population recursion and the translated maturity and weight schedules. | Differences in weights, timing, catch/index coverage, and population processes remain. |",
      "",
      "The 2025 management-track assessment is newer and supports current advice, but its public materials did not provide the age-specific surfaces required for this comparison. The 2018 record is the latest detailed accepted assessment represented in the database."
    )
  )
}
