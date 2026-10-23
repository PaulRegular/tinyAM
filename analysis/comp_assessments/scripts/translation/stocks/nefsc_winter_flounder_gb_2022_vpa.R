## Observations ----

assessment_id <- source$assessment$assessment_id[[1L]]
years <- 1982:2021
ages <- 1:7
surveys <- c("NEFSC spring", "DFO spring", "NEFSC fall")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  surveys = surveys,
  assumptions = source$assumptions
)

zero_index_cells <- sum(obs$index$obs == 0, na.rm = TRUE)
obs$index <- obs$index[is.finite(obs$index$obs) & obs$index$obs > 0, , drop = FALSE]
obs$index$q_key <- interaction(obs$index$survey, obs$index$age, drop = TRUE)

## Comparisons and background ----

age_plus_group <- 7
comparison_aggregates <- "ssb"
comparison_age_groups <- list(F = list(`4-6` = 4:6))
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "Annual SSB calculated from tinyAM abundance, the source maturity and stock-weight surfaces, and tinyAM spawning timing.",
    reason = "The VPA file contains MFSPAWN = 0.25 0.25, but its two values are not explained. The source and tinyAM spawning calculations cannot be confirmed as equivalent."
  ),
  recruitment = list(
    status = "matched",
    definition = "Age-1 abundance in the same calendar year.",
    reason = "Both reports identify recruitment at age 1."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Fishing mortality over ages 4-6.",
    reason = "The official report gives fully selected F for ages 4-6 but does not provide numerical F-at-age; a common population-weighted definition cannot be reconstructed."
  )
)

background <- c(
  "### Georges Bank winter flounder: 2022 VPA/ADAPT assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main limitation |",
  "|---|---|---|---|",
  "| Years and ages | VPA/ADAPT 3.1 RunSpring22 covers 1982-2021, ages 1-6 and 7+. | Fit 1982-2021 using ages 1-6 and age 7 as the plus group. | The full native age-year grid is preserved. |",
  "| Catch | Native VPA catch numbers-at-age are available for all years and ages. | Use the combined numbers-at-age input once. | The report documents missing Canadian discard data and expert-filled source gaps; the native model input does not recover those source decisions. |",
  "| N and F | VPA estimates annual population states; numerical N-at-age and F-at-age tables were not recovered. | Use an IID abundance process and an AR1 F process. | These process choices are tinyAM approximations; no age-specific official comparison is available. |",
  "| M | Fixed M = 0.3 per year for ages 1-7. | Use the same fixed age schedule. | This is an input consistency check, not an estimated-M comparison. |",
  paste0("| Surveys | The native file has NEFSC spring, DFO spring, and NEFSC fall index streams, labelled NUMBER and 1-JAN. The multiplier is not stated. | Use the three source series, with q by survey and age; retain only positive index cells for the lognormal observation model. | The source file does not distinguish true zero observations from missing cells. The excluded ", zero_index_cells, " zero cells remain in the canonical database. |"),
  "| Biology and SSB | Annual stock weights and maturity are available. The native BIOMASS and SSB weight surfaces are identical. | Use the annual stock-weight and maturity surfaces. | The meaning of the native MFSPAWN pair is unresolved, so SSB is approximate. |",
  "| Reported outputs | Table 1 reports SSB, age-1 recruitment, and FFull for 2012-2021; projection tables report 2022-2025. | Compare those available aggregate series only. | The report does not give numerical N-at-age, F-at-age, or historical uncertainty. |",
  "| Newer assessment | The 2025 WHAM assessment was accepted as the best scientific information available, but only summary tables were recovered. | Keep its inputs and outputs in a separate 2025 assessment record. | Do not combine 2025 outputs with the detailed 2022 VPA inputs. |",
  "",
  "The 2022 VPA source file labels every index stream at 1-JAN. tinyAM uses that native timing rather than inferring seasonal midpoints. Zero index cells are excluded only from the tinyAM lognormal fit because the source encoding does not identify whether they are missing or true zero observations. Catch observations and canonical source rows are not altered."
)

## Model ----

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~ factor(age), mean_ages = 4:6),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(q_form = ~ 0 + q_key, sd_form = ~ survey),
    silent = silent
  )
}
