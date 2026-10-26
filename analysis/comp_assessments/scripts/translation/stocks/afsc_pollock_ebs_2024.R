## Observations ----
years <- 1964:2024
ages <- 1:15
surveys <- c(
  "NMFS bottom-trawl VAST",
  "NMFS acoustic-trawl",

  # The official model treats ATS age 1 as a separate recruitment index.
  # Include that source stream here so database_to_tam_obs() imports it;
  # it is combined back into the main acoustic-trawl survey below for the
  # simpler tinyAM age-specific index representation.
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

## Treat very small reconstructed catches (< 0.02) as zeros to align with
## rounded zero cells in the official report. Zero cells do not enter the
## positive-observation lognormal catch likelihood.
obs$catch$obs[obs$catch$obs < 0.02] <- 0


## Apply spawning time used in the assessment.
obs$maturity$spawn_time <- 0.25

## The official assessment removes ATS age 1 from the ordinary age-composition
## likelihood and fits it separately as a recruitment index. For tinyAM, retain
## that information but treat it as age 1 of the same acoustic-trawl survey.
## The source-specific 2024 age-1 observation remains absent because it failed
## the official model's uncertainty/exclusion criterion during database import.
obs$index$survey[
  obs$index$survey == "NMFS acoustic-trawl age-1 index"
] <- "NMFS acoustic-trawl"

## Survey catchability is modelled below using separate survey baselines
## and a static logistic age-selectivity curve for each survey. The earlier
## survey-by-age q blocks are not used in this final specification.


## Background and comparisons ----

age_plus_group <- 15


comparison_age_groups <- list(N = list(`10+` = 10:15), biomass_at_age = list(`3+` = 3:15))

comparison_aggregates <- "ssb"

comparison_definitions <- list(ssb = list(
  status = "matched",
  definition = "Female SSB at year fraction 0.25, ages 1-15, using stock weights and 0.5-scaled maturity",
  reason = ""
))

background <- c(
  "### Eastern Bering Sea pollock: official 2024 Model 23.0", "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
  "| Years | Population model covers 1964–2024. | Fit 1964–2024. | The modeled period is retained. |",
  "| Ages | Ages 1–15, with age 15 as the model plus group; published N tables report ages 10–15 as 10+. | Ages 1–15, with age 15 as the plus group. | Reported 10+ abundance is compared with the sum of tinyAM ages 10–15. |",
  "| Recruitment and N | Recruitment varies annually; older cohorts are propagated through F and M without annual survival deviations. Initial ages 2–15 have a shared log mean and penalized age deviations. | Recruitment retains tinyAM's default log-recruitment random walk; older ages follow deterministic cohort survival (`N process = off`) with exponential initial abundance (`init = exp`). | Removing IID N deviations substantially improved scale agreement. The initial-age parameterization and recruitment regularization are not identical to the official model. |",
  "| F | Annual F level multiplied by time-varying, age-specific fishery selectivity, with penalties on selectivity changes. | Age-specific log F follows a correlated random walk (`cor_rw`) across years, with AR1 correlation among age-specific annual innovations; no fixed logistic F curve is imposed. | The correlation allows coherent F changes across ages but does not reproduce the official separable fishing-intensity and selectivity parameterization or its smoothing penalties. |",
  "| M | Fixed age-specific M: 0.9 at age 1, 0.45 at age 2, 0.3 at ages 3–15. | Use the same supplied M vector with `M process = off`. | Fixed M matches by construction. |",
  "| Catch | Separate total catch biomass and annually weighted fishery age-composition likelihoods; no 2024 fishery age composition. | Reconstruct catch-at-age in numbers using source catch weights and fit positive cells with independent lognormal errors and one common catch SD (`sd_form = ~ 1`). Reconstructed observations below 0.02 are set to zero. | A single age-specific catch likelihood replaces the official biomass-plus-composition components; zero cells do not contribute to the lognormal likelihood. |",
  "| Bottom-trawl survey | VAST biomass and age composition are fitted separately, including a full across-year covariance for the biomass index. Selectivity uses a time-varying logistic form and an age-1 adjustment. | Reconstruct a mid-year age-specific index using source survey weights; use a survey-specific q baseline and static logistic age curve. | Independent lognormal index observations replace the biomass covariance and age-composition likelihood; time-varying selectivity is omitted. |",
  "| Acoustic-trawl survey | ATS biomass and ages 2–15 composition are fitted separately with penalized, time-varying age selectivity. Age 1 is treated as a separate recruitment index; its 2024 observation is excluded. | Combine available ATS age-1 observations with other ATS age-specific observations, using a separate acoustic q baseline and static logistic age curve. Retain the source's 2024 age-1 exclusion. | The unified index and static logistic curve simplify the source likelihood and time-varying, nonparametric selectivity. |",
  "| Survey q | Separate BTS and ATS catchability parameters combined with distinct, partly time-varying selectivity structures. | `q_form = ~ 0 + survey + logistic(age, by = survey)` provides independent q baselines and static logistic age curves for each survey. | The previous discrete survey-by-age q blocks are no longer used; age and temporal selectivity details remain simplified. |",
  "| Index error | Biomass indices and age compositions use distinct likelihood components, variance structures, and sampling weights. | Independent lognormal age-specific index errors with one estimated log-SD per survey (`sd_form = ~ 0 + survey`). | Neither the VAST covariance nor native age-composition effective weights are reproduced. |",
  "| Weights, maturity and timing | Annual stock, catch, BTS and ATS weights-at-age are used; maturity is multiplied by 0.5 for female biomass, and SSB is evaluated at year fraction 0.25. | Retain the translated stock and source-specific observation weights, scale maturity by 0.5, and set `spawn_time = 0.25`. | The female fraction and spawning-time survival adjustment are now aligned. |",
  "| SSB | Female SSB at year fraction 0.25 over ages 1–15. | Female SSB computed at year fraction 0.25 over ages 1–15 using supplied stock weights, 0.5-scaled maturity and modeled F and M. | The comparison now uses a matching biomass definition; any estimate differences reflect the fitted states and remaining assessment-model differences, not a beginning-year SSB convention. |",
  "", "The final simplified tinyAM specification uses deterministic older-age cohort survival, default random-walk recruitment, correlated age-specific F random walks, fixed source M, one constant catch-observation SD, and static survey-specific logistic q curves. Earlier versions used IID N deviations; removing them was the main change behind improved scale agreement with the official assessment.",
  "", "The official acoustic-trawl assessment treats age 1 as a separate recruitment index. tinyAM includes available age-1 observations in the acoustic age-specific index but retains the 2024 age-1 exclusion. The official BTS and ATS biomass-plus-composition likelihoods and their time-varying selectivity structures are deliberately not reproduced.",
  "", "Historical fishery CPUE and acoustic-vessel-of-opportunity indices remain omitted because they lack the age-composition information used by this age-specific translation.",
  "", "Common output comparisons use official N ages 1–9 directly and reported 10+ against the sum of tinyAM ages 10–15; recruitment at age 1; total abundance over ages 1–15; and age-3+ biomass against the same age range in tinyAM. Female SSB now uses aligned spawning timing and maturity conventions. Official F-at-age was reconstructed from pinned model parameters and source code, but its uncertainty is unavailable. Fixed M matches by construction; the results are not a replication of the full operational likelihood."
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
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(
      q_form = ~ 0 + survey + logistic(age, by = survey),
      sd_form = ~ 0 + survey
    ),
    silent = silent
  )
}
