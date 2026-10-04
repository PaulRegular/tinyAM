translate_stock <- function(source) {
  assessment_id <- source$assessment$assessment_id[[1]]
  years <- 1970:2024
  ages <- 1:10
  surveys <- c(
    "Shelikof winter acoustic",
    "NMFS bottom trawl",
    "ADF&G crab/groundfish trawl",
    "Summer acoustic"
  )

  obs <- database_to_tam_obs(
    assessment_id,
    source$inputs,
    years = years,
    ages = ages,
    weight_survey = "",
    index_weight_source = "source",
    maturity_multiplier = 0.5,
    surveys = surveys,
    assumptions = source$assumptions
  )

  obs$catch$obs[obs$catch$age == 2] <- NA_real_
  obs$catch$age_factor <- factor(obs$catch$age, levels = ages)
  shelikof <- obs$index$survey == "Shelikof winter acoustic"
  obs$index <- obs$index[!(shelikof & obs$index$age == 3), , drop = FALSE]

  environmental <- source$inputs[
    source$inputs$type == "covariate" &
      source$inputs$measure == "environmental_covariate",
    , drop = FALSE
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

  adfg <- "ADF&G crab/groundfish trawl"
  adfg_years <- sort(unique(obs$index$year[obs$index$survey == adfg]))
  annual_q_terms <- character()
  if (length(adfg_years) > 1L) {
    for (year in adfg_years[-1L]) {
      term <- paste0("adfg_q_", year)
      obs$index[[term]] <- as.numeric(obs$index$survey == adfg &
                                       obs$index$year == year)
      annual_q_terms <- c(annual_q_terms, term)
    }
  }
  obs$index$q_key <- interaction(obs$index$survey, obs$index$age,
                                 drop = TRUE, lex.order = TRUE)

  q_terms <- c("q_key", "environmental_effect", annual_q_terms)
  q_form <- stats::reformulate(q_terms, intercept = FALSE)

  background <- c(
    "## Gulf of Alaska pollock: accepted 2024 Model 23d",
    "",
    "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
    "|---|---|---|---|",
    "| Years | 1970-2024; the accepted age-structured model is the Western/Central/West Yakutat stock. | Fit 1970-2024. | The accepted model period and stock area are retained. |",
    "| Ages | Ages 1-10+, with recruitment at age 1. | Ages 1-10, with age 10 as the plus group. | The terminal group is retained. |",
    "| N | Recruitment varies with fixed SD 1.3; older fish survive deterministically, and initial ages are tied to first-year recruitment and M. | Exponential initial abundance, no older-age process, and tinyAM's estimated random-walk recruitment process. | tinyAM cannot fix recruitment SD or reproduce the source's exact initial-state construction. |",
    "| F | One fishery with double-logistic selectivity; ascending selectivity parameters change annually with penalties. | Age-specific mean log F with an AR1 process over ages and years. | This is a smooth age-time approximation, not the source's selectivity parameterization. |",
    "| M | Fixed external age-specific M; the accepted model fixes its scalar at 1. | The same supplied M vector with the M process off. | The fixed mortality input is retained. |",
    "| Catch | Total catch biomass plus number compositions; the first composition bin combines ages 1-2 and age 10 is 10+. | Reconstructed catch numbers-at-age, fitted with independent lognormal errors; the grouped ages 1-2 bin is not fitted as age 2. | tinyAM has no grouped catch prediction or separate total-catch and Dirichlet-multinomial likelihood. |",
    "| Index | Four active biomass indices with periodic age compositions; Shelikof age-1/2 indices are disabled and its first composition bin pools ages 1-3. Shelikof q includes an environmental effect and ADF&G q varies annually. | Age-specific numbers reconstructed with matching survey weights; the pooled Shelikof ages 1-3 bin is omitted rather than treated as age 3, and remaining reported ages are retained. Age-specific q includes the observed Shelikof covariate and ADF&G annual effects on composition years. | Aggregate-only years, latent environmental dynamics, q penalties, supplied aggregate-index SDs, grouped-age predictions, and the source composition likelihood are not represented. |",
    "| Weights and maturity | Annual stock, catch, and survey weights; a constant maturity vector and female fraction 0.5. | Annual stock weights, source survey/catch weights for conversions, and source maturity multiplied by 0.5. | Separate weight purposes are retained during conversion. |",
    "",
    "Catch and survey compositions are converted to numbers-at-age with their corresponding annual total and weights. The accepted source reports population numbers, recruitment, and SSB in million fish and thousand tonnes; comparison scales convert tinyAM outputs to those units.",
    "",
    "The transcribed assessment outputs do not include an age-specific F surface, so the accepted F values are unavailable in the comparison dashboard."
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 10,
    obs = obs,
    comparison_scales = c(
      N = 1e-6, recruitment = 1e-6, abundance = 1e-6,
      ssb = 1e-6, biomass = 1e-6, biomass_at_age = 1e-6
    ),
    settings = list(
      N_settings = list(process = "off", init = "exp"),
      F_settings = list(process = "ar1", mu_form = ~ 0 + age_factor,
                        mean_ages = 3:10),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption, mean_ages = 3:10),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = q_form,
                            sd_form = ~ 0 + survey,
                            fill_missing = TRUE)
    ),
    background = background
  )
}
