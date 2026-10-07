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

  catch_numbers <- source$inputs[
    source$inputs$type == "catch" & source$inputs$measure == "numbers_at_age",
    , drop = FALSE
  ]
  if (!nrow(catch_numbers) || !setequal(catch_numbers$age, 1:15)) {
    cli::cli_abort("GOA pollock needs the published catch numbers for ages 1-15.")
  }
  translation_inputs <- source$inputs[
    !(source$inputs$type == "catch" & source$inputs$measure == "proportion_at_age"),
    , drop = FALSE
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
  obs$index$adfg_year <- factor(
    ifelse(obs$index$survey == adfg, obs$index$year, adfg_years[1]),
    levels = adfg_years
  )
  trawl <- obs$index$survey %in% c("NMFS bottom trawl", adfg)
  obs$index$q_order <- ifelse(trawl, obs$index$age, 11 - obs$index$age)
  q_form <- ~ survey + mono(q_order, by = survey) + environmental_effect + adfg_year

  catch_weights <- source$inputs[source$inputs$type == "catch_weight" &
                                  source$inputs$measure == "weight_at_age", ]
  catch_totals <- source$inputs[source$inputs$type == "catch" &
                                 source$inputs$measure == "total_biomass", ]
  catch_reporting <- list(
    weights = data.frame(year = catch_weights$year, age = catch_weights$age,
                        weight = mapply(.translation_weight_to_kg,
                                        catch_weights$value, catch_weights$unit)),
    totals = data.frame(year = catch_totals$year,
                       yield = mapply(.translation_biomass_to_kg,
                                      catch_totals$value, catch_totals$unit))
  )

  background <- c(
    "### Gulf of Alaska pollock: accepted 2024 Model 23d",
    "",
    "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
    "|---|------|------|------|",
    "| Years | 1970-2024; the accepted age-structured model is the Western/Central/West Yakutat stock. | Fit 1970-2024. | The accepted model period and stock area are retained. |",
    "| Ages | Ages 1-10+, with recruitment at age 1. | Ages 1-10, with age 10 as the plus group. | The terminal group is retained. |",
    "| N | Recruitment varies with fixed SD 1.3; older fish survive deterministically, and initial ages are tied to first-year recruitment and M. | Exponential initial abundance, no older-age process, and tinyAM's estimated random-walk recruitment process. | tinyAM cannot fix recruitment SD or reproduce the source's exact initial-state construction. |",
    "| F | One fishery with double-logistic selectivity; ascending selectivity parameters change annually with penalties. | Age-specific mean log F with an AR1 process over ages and years. | This is a smooth age-time approximation, not the source's selectivity parameterization. |",
    "| M | Fixed external age-specific M; the accepted model fixes its scalar at 1. | The same supplied M vector with the M process off. | The fixed mortality input is retained. |",
    "| Catch | Total catch biomass plus number compositions; the accepted likelihood pools ages 1-2 and ages 10+, while the detailed report table gives catch numbers separately for ages 1-15. | Published catch numbers at ages 1-15 with a common log-SD; ages 10-15 are summed into tinyAM age 10+. The yield panel compares original total catch biomass with all-age predictions using catch weights. | The report table is rounded to 0.01 million fish; tinyAM does not reproduce the source's pooled composition, age-reading error, total-catch or Dirichlet-multinomial likelihood. |",
    "| Index | Four active biomass indices with periodic age compositions; Shelikof age-1/2 indices are disabled and its first composition bin pools ages 1-3. Shelikof q includes an environmental effect and ADF&G q varies annually. | Age-specific numbers reconstructed with matching survey weights; the pooled Shelikof ages 1-3 bin is omitted rather than treated as age 3, and remaining reported ages are retained. Catchability rises with age for the two trawls and declines with age for both acoustic surveys, using independent monotone steps that can form plateaus. It includes the observed Shelikof covariate and ADF&G annual effects on composition years. | Monotone steps approximate the source logistic limbs; there is no source rule that explicitly pools older ages. Aggregate-only years, latent environmental dynamics, q penalties and the bottom-trawl q prior, supplied aggregate-index SDs, grouped-age predictions, age-reading error, and the source composition likelihood are not represented. |",
    "| Weights and maturity | Annual stock, catch, and survey weights; a constant maturity vector and female fraction 0.5. | Annual stock weights, source survey/catch weights for conversions, and source maturity multiplied by 0.5. | Separate weight purposes are retained during conversion; tinyAM SSB uses stock weights at the start of the year, while accepted SSB uses spawning weights after survival to year fraction 0.21. |",
    "",
    "Catch and survey compositions are converted to numbers-at-age with their corresponding annual total and weights. The accepted source reports population numbers, recruitment, and SSB in million fish and thousand tonnes; comparison scales convert tinyAM outputs to those units.",
    "",
    "The detailed report table gives catch numbers by age 1-15 for 1975-2023. Values are rounded to 0.01 million fish, so printed zeros are below the table's reporting precision; tinyAM's lognormal catch likelihood treats zeros as missing. The accepted model's native composition still pools ages 1-2 and 10+, and its age-reading-error matrix is not represented in tinyAM. Yield predictions sum conditional median catches; they do not include a lognormal mean correction. The accepted F surface is unavailable."
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 10,
    obs = obs,
    catch_reporting = catch_reporting,
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
