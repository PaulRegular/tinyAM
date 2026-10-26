## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 1970:2024
ages <- 1:10
surveys <- c(
  "Shelikof winter acoustic",
  "NMFS bottom trawl",
  "ADF&G crab/groundfish trawl",
  "Summer acoustic"
)

catch_numbers <- source$inputs[
  source$inputs$type == "catch" & source$inputs$measure == "numbers_at_age", ,
  drop = FALSE
]
if (!nrow(catch_numbers) || !setequal(catch_numbers$age, 1:15)) {
  cli::cli_abort("GOA pollock needs the published catch numbers for ages 1-15.")
}
translation_inputs <- source$inputs[
  !(source$inputs$type == "catch" & source$inputs$measure == "proportion_at_age"), ,
  drop = FALSE
]

obs <- database_to_tam_obs(
  assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  index_weight_source = "source",
  maturity_multiplier = 0.5,
  surveys = surveys,
  assumptions = source$assumptions
)

shelikof <- obs$index$survey == "Shelikof winter acoustic"
obs$index <- obs$index[!(shelikof & obs$index$age == 3), , drop = FALSE]

environmental <- source$inputs[
  source$inputs$type == "covariate" &
    source$inputs$measure == "environmental_covariate", ,
  drop = FALSE
]
environmental_years <- as.integer(environmental$year)
environmental_values <- as.numeric(environmental$value)
environmental_effect <- environmental_values[
  match(obs$index$year, environmental_years)
]
is_shelikof <- obs$index$survey == "Shelikof winter acoustic"
if (anyNA(environmental_effect[is_shelikof])) {
  cli::cli_abort("Shelikof age-composition years are missing the source environmental covariate.")
}
obs$index$environmental_effect <- ifelse(is_shelikof, environmental_effect, 0)

# Only ADF&G catchability receives an annual random-walk effect.
adfg <- "ADF&G crab/groundfish trawl"
obs$index$adfg_indicator <- as.numeric(obs$index$survey == adfg)

# Reverse the age coordinate for the two acoustic surveys, producing descending
# logistic age responses; the trawl surveys retain rising logistic responses.
# Both acoustic surveys share one curve, both trawl surveys share another;
# ~0 + survey below estimates separate survey-specific q baselines.
descending_surveys <- c("Shelikof winter acoustic", "Summer acoustic")
obs$index$sel_age <- ifelse(
  obs$index$survey %in% descending_surveys,
  min(ages) + max(ages) - obs$index$age,
  obs$index$age
)
obs$index$survey_type <- ifelse(
  obs$index$survey %in% descending_surveys,
  "acoustic", "trawl"
)

catch_weights <- source$inputs[source$inputs$type == "catch_weight" &
  source$inputs$measure == "weight_at_age", ]
catch_totals <- source$inputs[source$inputs$type == "catch" &
  source$inputs$measure == "total_biomass", ]
catch_reporting <- list(
  weights = data.frame(
    year = catch_weights$year, age = catch_weights$age,
    weight = mapply(
      .translation_weight_to_kg,
      catch_weights$value, catch_weights$unit
    )
  ),
  totals = data.frame(
    year = catch_totals$year,
    yield = mapply(
      .translation_biomass_to_kg,
      catch_totals$value, catch_totals$unit
    )
  )
)

background <- c(
  "### Gulf of Alaska pollock: official 2024 Model 23d",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM candidate | Remaining difference |",
  "|---|---|---|---|",
  "| Years and ages | Western/Central/West Yakutat, 1970-2024; ages 1-10+, recruitment at age 1. | Same years, ages, and age-10 plus group. | Model period and age range align. |",
  "| Recruitment and N | Recruitment deviations around a common mean with fixed log-SD 1.3; deterministic older-age survival; initial ages linked to recruitment and M. | IID recruitment about a constant log mean with SD fixed at 1.3; N process off and exponential initialization. | Recruitment variance and temporal structure are closer; first-year recruitment and exact initial-age construction still differ. |",
  "| F | Annual fishing intensity times double-logistic fishery selectivity, with penalized changes to ascending-limb parameters. | Rising logistic mean log-F-at-age with age/year AR1 deviations; reported mean F ages 3-10. | Static single-logistic mean plus AR1 residuals differs from time-varying double-logistic selectivity. |",
  "| M | Fixed externally supplied age-specific M, scalar fixed at 1. | Same supplied M vector; M process off. | Mortality values retained. |",
  "| Catch | Total catch biomass plus age compositions using Dirichlet-multinomial likelihood; young and terminal ages accumulated. | Published catch numbers at ages 1-15 summed to modeled 10+, fitted with a common lognormal catch-at-age SD. | Observation decomposition, aggregation and likelihood differ; published zeros are missing in tinyAM. |",
  "| Survey selectivity | Shelikof descending logistic with ages 1-2 excluded; NMFS bottom trawl and summer acoustic double logistic; ADF&G rising logistic. | Two shared logistic age responses on logit-q: descending for Shelikof/summer through reversed age, ascending for the two trawl surveys; a distinct q baseline for each survey. | Curves are shared by survey type, and single logistics approximate the source double-logistic responses. |",
  "| Survey q over time | Shelikof q linked to a latent environmental AR1 process; ADF&G annual log-q RW; bottom trawl and summer acoustic q constant, with a prior on bottom-trawl q. | Observed Shelikof environmental covariate enters logit-q; a single ADF&G-specific rw(year) enters logit-q; other survey baselines are constant. | No latent environmental measurement model, source year-specific RW SDs, bottom-trawl prior or native log-q link. |",
  "| Index observations | Separate survey biomass indices and age compositions with survey timing, age-reading error and source-specific likelihoods and weights. | Age-specific indices reconstructed from biomass/composition and survey weights, retaining source timing; supplied log-SDs plus one estimated scale per survey. | Independent lognormal age-specific observations replace native biomass/composition likelihoods; aggregate-only and pooled Shelikof age-1-3 rows omitted. |",
  "| Weights, maturity and SSB | Annual spawning weights, survival to year fraction 0.21, and female fraction 0.5. | Stock weights, survey-specific weights for index reconstruction, maturity multiplied by 0.5; beginning-year tinyAM SSB. | Source spawning weights and survival-to-spawning adjustment are not reproduced in this candidate. |",
  "",
  "The translation retains the 1970-2024 model years and supplied age-specific natural mortality. Its candidate specification uses deterministic older-age survival, IID recruitment (fixed log-SD 1.3), an AR1 fishing-mortality residual around a logistic age pattern, and bounded survey q with type-shared logistic age patterns. Only ADF&G q receives a random walk; the Shelikof observed environmental covariate enters directly. This is a parsimonious approximation, not a reproduction of the official likelihood.",
  "",
  "Because q_link = 'logit', the age curves and temporal effects enter the logit-q predictor, unlike the official log-q parameterization. The older paired-age q blocks and separate ADF&G annual factor effects are no longer used in the fit. Weak q identification and sensitivity to the link should be checked after refitting.",
  "",
  "The translation uses thousand fish, average weights in kg, and biomass in tonnes. Original catch totals are used for the yield panel; predicted catch numbers are multiplied by catch weights. The published 1975-2023 catch table is rounded to 0.01 million fish and printed zeros are missing under tinyAM's lognormal likelihood. Predictions are conditional medians without lognormal mean correction. Direct official F-at-age estimates are unavailable; the database includes a clearly labeled, partial reconstruction from rounded N-at-age and fixed M for ages 1-8 in 1970-2023, without uncertainty.",
  "",
  "SSB remains definition-sensitive: the official calculation uses spawning-time biomass at fraction 0.21 with spawning weights, while the shared-biology comparison uses source N with translated stock weights and maturity. The current package-default fit converged (optimizer code 0; objective 1356.997; maximum absolute gradient 0.000219; successful uncertainty estimation; positive-definite Hessian), although residual checks flag location, spread and outliers. For 1970-2024, mean absolute percentage differences are 27.7% for total abundance, 46.5% for N-at-age, 44.2% for recruitment, 37.1% for shared-biology SSB and 30.2% for biomass. Their terminal differences are +32.4%, +77.0%, -2.0%, +8.5% and +41.1%; trend correlations are 0.949, 0.943, 0.957, 0.881 and 0.905. The partial F-at-age comparison is approximate because it derives F from rounded N-at-age and fixed M; mean absolute percentage difference is 82.9%, terminal difference -21.6%, and trend correlation 0.625. Fixed M matches by construction, and official Fbar and age-specific biomass are unavailable. No source-start retry or assumption changes were needed."
)


## Background and comparisons ----

age_plus_group <- 10
comparison_outputs <- source$outputs

comparison_definitions <- list(
  F = list(
    status = "approximate",
    definition = "Approximate total F-at-age derived for ages 1-8 from rounded official N-at-age and fixed M",
    reason = "F is reconstructed from the official cohort equation using N-at-age rounded to the nearest million fish. Only rows remaining positive across the rounding bounds are retained; age 9, the age-10+ plus group, and near-zero cells are unavailable. No uncertainty is available."
  )
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
    rec_settings = list(process = "iid", mu_form = ~1, sd = 1.3),
    F_settings = list(process = "ar1", mu_form = ~logistic(age), mean_ages = 3:10),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption, mean_ages = 3:10),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~ 0 + survey + logistic(sel_age, by = survey_type) +
        environmental_effect +
        rw(year, by = adfg_indicator),
      q_link = "logit",
      sd_form = ~ 0 + survey,
      sd_supplied = ~ relative_sd
    ),
    silent = silent
  )
}
