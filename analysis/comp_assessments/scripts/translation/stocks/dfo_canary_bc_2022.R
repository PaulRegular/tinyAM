# BC Canary Rockfish (Run24 / 24.01): exploratory tinyAM translation.
# Relies on the existing, source-validated helper. No source model outputs
# are used to construct catch or survey observations.
source("analysis/comp_assessments/R/bc_canary_obs.R", local = TRUE)

## Observations ----
years <- 1977:2022
source_ages <- 1:60                 # Preserve original sampled ages and all biology.
ages <- 1:40                        # Trial a 40+ model, with hidden biological ages 40:60+.

catch <- .canary_female_numbers(source$inputs, "catch")
index <- .canary_female_numbers(source$inputs, "index")

# Select official sex-specific stock weight and effective reproductive weight.
weight <- source$inputs[
  source$inputs$type == "weight" &
    source$inputs$sex == "female" &
    !is.na(source$inputs$season) &
    source$inputs$season == "begin_year", , drop = FALSE
]
maturity <- source$inputs[
  source$inputs$measure == "maturity_fecundity_at_age" &
    source$inputs$sex == "female", , drop = FALSE
]
maturity$value <- maturity$value /
  weight$value[match(maturity$age, weight$age)]
maturity$type <- "maturity"
maturity$measure <- "maturity_at_age"
maturity$basis <- maturity$unit <- "proportion"
maturity$notes <- paste(
  "Effective female maturity from source reproductive weights divided by",
  "January population weight; retains the executable age-10 spawning cutoff.",
  "Female abundance is modeled directly; no additional 0.5 factor."
)

# Preserve the accepted model's supplied, year-specific index log SDs.
# The existing helper reconstructs female age-specific index observations but
# does not bring these aggregate-index uncertainty records along with them.
# The translator joins them by survey and year as obs$index$relative_sd.
precision <- source$inputs[
  source$inputs$type == "index" &
    source$inputs$measure == "log_index_sd" &
    source$inputs$survey %in% unique(index$survey) &
    source$inputs$year %in% years, , drop = FALSE
]
if (!nrow(precision)) {
  stop("Missing native index log SDs for the retained Canary surveys.", call. = FALSE)
}

translation_inputs <- rbind(catch, index, weight, maturity, precision)
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  translation_inputs,
  years = years,
  ages = source_ages,              # Retain full biology for shortened plus group.
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0
obs$pop_covariates$M_prior_group <- factor("female_stock")

obs$index$age_group <- cut_ages(obs$index$age, c(seq(min(obs$index$age), 30, by = 5), max(obs$index$age)))

# Inspect the index SD join before fitting: the source index SD describes the
# aggregate survey observation, NOT independent age-composition uncertainty.
if (!"relative_sd" %in% names(obs$index) ||
    any(!is.finite(obs$index$relative_sd) | obs$index$relative_sd <= 0)) {
  stop("Source index log SDs were not matched to all reconstructed age rows.",
       call. = FALSE)
}

## Background and comparisons ----
age_plus_group <- max(ages)
comparison_aggregates <- "ssb"
comparison_outputs <- source$outputs[source$outputs$sex == "female", ]
comparison_definitions <- list(
  ssb = list(
    status = "approximate",
    definition = "January 1 female mature biomass, ages 1-40+ (biology retained through 60+)",
    reason = paste(
      "Age 0 is immature; the 40+ group retains older biology internally.",
      "Official annual SSB uses the public 500-draw subset; tinyAM is likelihood-based."
    )
  ),
  M = list(
    status = "matched",
    definition = "Female constant instantaneous M, repeated across model ages/years",
    reason = paste(
      "One estimated female M is compared with one official posterior median,",
      "not independent age-year estimates. The natural-normal prior and",
      "full posterior integration are approximated."
    )
  )
)

background <- c(
  "### British Columbia Canary Rockfish: official 2022 Run24 / 24.01",
  "", print_sources(source$assessment), "",
  "| Component | Official assessment | Current tinyAM translation | Principal difference / limitation |",
  "|---|---|---|---|",
  "| Years / ages | 1935–2022; ages 0–60+, two sexes and age-0 recruitment. | Fit 1977–2022, female ages 1–40+, while retaining biological inputs through 60+. | A shorter initial history and pooled 40+ state; tinyAM retains older-age weights and maturity internally but does not reproduce SS3's separate age states. |",
  "| Population | Sex-specific deterministic cohort survival, equilibrium initialization and early recruitment deviations. | Female-only cohort survival (`N process = 'off'`), exponential initial survivorship (`init = 'exp'`). | Male population and the native early history are not modeled. Joint-sex observations are converted before retaining the female component; no additional 0.5 female factor is applied. |",
  "| Recruitment | Beverton–Holt at age 0, fixed log SD 0.9, constrained deviations and bias corrections. | Age-1 constant-mean IID log recruitment, fixed SD 0.9. | Omits juvenile survival, the stock–recruit relationship and the source's deviation and terminal-year rules. |",
  "| Fishing mortality | Two separate annual fleet intensities (trawl and other) with distinct, time-invariant, approximately asymptotic double-normal selectivities; Pope catch approximation. | One combined F-at-age surface with a rising logistic mean and age/year AR1 residuals (`process = 'ar1', mu_form = ~ logistic(age)`); Baranov catches. | Flexible F departures may approximate changing aggregate age patterns arising from shifts in fleet contributions, but this mechanism has not been verified. Unlike SS3, tinyAM does not fit fleet-specific F or condition annual F on all native catches; the independent age–year deviations are replaced by an AR1 covariance. |",
  "| Natural mortality | Constant female M estimated with a natural-scale Normal prior (mean 0.06, SD 0.018; bounds 0.02–0.2), independently of male M. | A single estimated female log-M effect centred at log(0.06) with fixed SD 0.3. | The lognormal prior approximates the source CV; native bounds and posterior integration are not reproduced. The published female M is a reference output, not an observation or fixed input. |",
  "| Catch | Two biomass catch fleets; joint-sex trawl age compositions in 36 years (1977–2017). | Joint-sex source proportions and sex-specific weights reconstruct female catch-at-age where age samples exist, pooling total trawl and other catch biomass. | The unsampled fleet is assigned the trawl composition in sampled years; no age composition is invented in other years. In particular, catches in 2018–2022 do not directly constrain F, unlike in SS3. Native composition likelihood and ageing error are omitted. |",
  "| Survey / q | Six biomass surveys plus commercial relative CPUE; sampled QCS, WCVI and NMFS have profiled scalar q and static double-normal age selectivity. | Female index-at-age for those three surveys; survey-specific rising logistic effective q on the log link (`~ 0 + survey + logistic(age, by = survey)`). | Aggregate-only indices are omitted. The fitted effective q is not factored into a separately normalized scalar q and selectivity; no time-varying q or `q < 1` restriction is imposed. |",
  "| Index uncertainty | Final accepted run supplies year-specific, reweighted aggregate-index log SDs and separate age-composition uncertainty. | Native log SDs are joined and checked in `obs$index$relative_sd`, but the current fit uses `sd_form = ~ 0 + survey` WITHOUT `sd_supplied`; survey-specific lognormal SDs are estimated. | The supplied SDs are retained as data but do not enter this candidate's likelihood. Reconstructed age rows share each survey year's biomass total, whereas tinyAM treats them as conditionally independent. |",
  "| Biology / SSB | Fixed female growth/weight and effective reproductive weights with a first-mature-age cutoff at 10; January female spawning biomass. | Original female beginning-year weight and effective maturity through age 60; 40+ internal biology; no 0.5 multiplier; one-year projection. | Shortened age states and single-sex model remain approximations, although the input biological functions and spawning-time convention are preserved. |",
  "", "Interpretation of the latest fitted dashboard: both trajectories show a historical decline, but the accepted assessment turns upward in the 2000s whereas tinyAM continues declining into the mid-2010s and recovers only modestly toward 2022. The latest fit therefore captures much of the broad magnitude/pattern but underestimates recent SSB. The user also reports improved agreement for female M. These statements are qualitative, based on the supplied dashboard and report, not recomputed summary statistics or independent fit diagnostics.",
  "", "The flexible age/year AR1 F process may partially accommodate temporal changes in effective selectivity caused by changing contributions of fleets with different fixed selectivity curves. This is a modeling hypothesis rather than a demonstrated cause of improvement; the source other-fishery component is comparatively small over this window. An explicit two-fleet fit would be needed to attribute differences to fleet mixing.",
  "", "Historical diagnostics reported for a superseded 60+ formulation (objective 1726.629, terminal SSB difference -59.09%, female M difference -36.51%) MUST NOT be attributed to the current 40+ AR1 candidate. Updated numerical metrics and optimizer diagnostics have not been provided with this translation. The comparison uses official public 500-draw SSB subset summaries; its credible limits are not tinyAM confidence intervals.",
  "", "The main unresolved structural issue is that catch totals without corresponding trawl age samples (including 2018–2022), aggregate-only indices/CPUE, joint-sex composition covariance, and fleet-specific selectivity/catch equations are not reproduced. These differences could contribute to terminal SSB divergence; no single cause has been established."
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
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "iid", mu_form = ~ 1, sd = 0.9),
    F_settings = list(process = "ar1", mu_form = ~logistic(age)),
    M_settings = list(
      process = "off",
      mu_form = ~ 0 + iid(M_prior_group, sd = 0.3),
      mu_supplied = ~ I(0.06)
    ),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(
      q_form = ~ 0 + survey + logistic(age, by = survey),
      q_link = "log",
      sd_form = ~ 0 + survey
    ),
    proj_settings = list(n_proj = 1),
    silent = silent
  )
}
