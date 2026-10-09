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
  obs$index$adfg_year <- droplevels(cut_years(
    ifelse(obs$index$survey == adfg, obs$index$year, adfg_years[1]),
    seq(min(adfg_years), max(adfg_years))
  ))
  obs$index$q_age_block <- cut_ages(obs$index$age, c(1, 3, 5, 7, 9, 10))
  levels(obs$index$q_age_block) <- c("1-2", "3-4", "5-6", "7-8", "9-10", "9-10")

  obs$index$q_key <- interaction(
    obs$index$survey,
    obs$index$q_age_block,
    drop = TRUE
  )

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
    print_sources(source$assessment),
    "",

    "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
    "|---|------|------|------|",
    "| Years | Western/Central/West Yakutat stock, 1970-2024. | Same stock and period. | Retained. |",
    "| Ages | Ages 1-10+, with recruitment at age 1. | Ages 1-10, with age 10 as the plus group. | The terminal group is retained. |",
    "| N | Variable recruitment with fixed SD 1.3, deterministic older-age survival, and initial abundance based on recruitment and M. | Random-walk recruitment with estimated SD, deterministic older-age survival, and exponential initial abundance. | Recruitment and initialization differ. |",
    "| F | One fishery with double-logistic selectivity and penalized annual changes. | Independent temporal random walks in log F by age; summaries use ages 3-10. | Replaces the source selectivity model. |",
    "| M | Fixed external age-specific M; the accepted model fixes its scalar at 1. | The same supplied M vector with the M process off. | The fixed mortality input is retained. |",
    "| Catch | Total biomass and age compositions, with ages 1-2 and 10+ pooled in the likelihood. | Detailed published numbers at ages 1-15, summed to 10+, with a common log-SD. | Uses the finer reporting table and a lognormal likelihood rather than the source composition likelihood. |",
    "| Index | Shelikof winter acoustic, NMFS bottom trawl, ADF&G trawl, and summer acoustic biomass indices with age compositions and survey selectivity curves. | Numbers-at-age reconstructed with matching survey weights and source timing. Each survey has paired-age q blocks (1-2, 3-4, 5-6, 7-8, 9-10) and a logit link. Shelikof environmental and ADF&G annual effects act on logit-q; SD combines supplied log-SDs with one estimated level per survey. | Simplifies selectivity and observation errors. Aggregate-only years and the pooled Shelikof ages 1-3 bin are omitted. Source age-reading error, latent environmental dynamics, q penalties and priors are not represented. |",
    "| Weights and maturity | Annual stock, catch, survey and spawning weights; constant maturity and female fraction 0.5. | Stock weights, matching catch/survey weights, and maturity multiplied by 0.5. | tinyAM SSB uses start-year stock weights; accepted SSB uses spawning weights and survival to year fraction 0.21. |",
    "",
    "tinyAM fits individual fish with weights in kg. Accepted N and recruitment remain in million fish, and accepted SSB and biomass in thousand tonnes; comparisons convert tinyAM outputs to these units.",
    "",
    "The logit link restricts q below one; it does not reproduce the source bottom-trawl q prior. Bottom-trawl q for ages 5-8 approaches this boundary and is weakly estimated.",
    "",
    "Catch numbers for 1975-2023 are rounded to 0.01 million fish; printed zeros are treated as missing by the lognormal likelihood. The yield panel uses original catch biomass and predicted catches multiplied by catch weights. Predictions are conditional medians, without a lognormal mean correction. Accepted F-at-age is unavailable."
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
      F_settings = list(process = "rw", mu_form = NULL,
                        mean_ages = 3:10),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption, mean_ages = 3:10),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 0 + q_key + environmental_effect + adfg_year,
                            q_link = "logit",
                            sd_form = ~ 0 + survey,
                            sd_supplied = ~relative_sd,
                            fill_missing = FALSE)
    ),
    background = background
  )
}
