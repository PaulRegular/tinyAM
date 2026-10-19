## Observations ----
years <- 1978:2023
ages <- 2:11
inputs <- source$inputs

catch_input <- inputs[
  inputs$type == "catch" & inputs$measure == "numbers_at_age" &
    inputs$season == "fall" & inputs$fleet == "Fixed gear", ,
  drop = FALSE
]
catch_values <- aggregate(value ~ year + age, catch_input, sum)
catch_input <- inputs[rep(NA_integer_, nrow(catch_values)), , drop = FALSE]
catch_input$assessment_id <- source$assessment$assessment_id
catch_input$type <- "catch"
catch_input$measure <- "numbers_at_age"
catch_input$basis <- "numbers"
catch_input$fleet <- "Fixed gear, regions summed"
catch_input$survey <- catch_input$sex <- catch_input$region <- ""
catch_input$season <- "fall"
catch_input$year <- catch_values$year
catch_input$year_basis <- "calendar_year"
catch_input$age <- catch_values$age
catch_input$value <- catch_values$value
catch_input$unit <- "thousand fish"
catch_input$source_type <- "translation_assumption"
catch_input$source_reference <- "DFO Research Document 2024/058, Table 6"
catch_input$transformation <- "Summed fixed-gear catch-at-age across North, Middle, and South regions."
catch_input$notes <- paste(
  "The source has three regional populations. The 2024 mobile catch table",
  "is excluded because it repeats fixed-gear values."
)

index_input <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age", ,
  drop = FALSE
]
cpue <- index_input$survey == "Fall fixed-gear CPUE"
index_input$survey[cpue] <- paste(
  index_input$survey[cpue], index_input$region[cpue], sep = " - "
)

weight_input <- inputs[
  inputs$type == "catch_weight" & inputs$measure == "weight_at_age" &
    inputs$season == "fall" & is.finite(inputs$value) & inputs$value > 0, ,
  drop = FALSE
]
pooled_weight <- aggregate(value ~ year + age, weight_input, mean)
weight_matrix <- matrix(
  NA_real_, length(years), length(ages),
  dimnames = list(year = years, age = ages)
)
weight_matrix[cbind(
  match(pooled_weight$year, years), match(pooled_weight$age, ages)
)] <- pooled_weight$value

for (j in seq_along(ages)) {
  observed <- is.finite(weight_matrix[, j])
  weight_matrix[, j] <- stats::approx(
    years[observed], weight_matrix[observed, j], xout = years, rule = 2
  )$y
}

stock_weight <- weight_matrix
stock_weight[-1, -1] <- sqrt(
  weight_matrix[-length(years), -length(ages)] * weight_matrix[-1, -1]
)

weight_grid <- expand.grid(year = years, age = ages)
weight_rows <- inputs[rep(NA_integer_, nrow(weight_grid)), , drop = FALSE]
weight_rows$assessment_id <- source$assessment$assessment_id
weight_rows$type <- "weight"
weight_rows$measure <- "weight_at_age"
weight_rows$basis <- "kg_per_fish"
weight_rows$fleet <- weight_rows$survey <- weight_rows$sex <- weight_rows$region <- ""
weight_rows$season <- "fall"
weight_rows$year <- weight_grid$year
weight_rows$year_basis <- "calendar_year"
weight_rows$age <- weight_grid$age
weight_rows$value <- as.vector(stock_weight)
weight_rows$unit <- "kg per fish"
weight_rows$source_type <- "translation_assumption"
weight_rows$source_reference <- paste(
  "DFO Research Document 2024/058, Tables 7 and 11;",
  "DFO Research Document 2022/068, p. 6"
)
weight_rows$transformation <- paste(
  "Arithmetic mean of available fixed/mobile and regional source weights",
  "for each year-age; gaps interpolated by age. For ages above 2 and years",
  "after 1978, beginning-year weight is the geometric mean of age a-1 in",
  "year t-1 and age a in year t. Edge cells retain same-year source weight."
)
weight_rows$notes <- paste(
  "Fit-only pooled stock weights; the source model fits three regions and",
  "the aggregation of regional/gear weights is not fully specified."
)

maturity_grid <- expand.grid(year = years, age = ages)
maturity_rows <- inputs[rep(NA_integer_, nrow(maturity_grid)), , drop = FALSE]
maturity_rows$assessment_id <- source$assessment$assessment_id
maturity_rows$type <- "maturity"
maturity_rows$measure <- "maturity_at_age"
maturity_rows$basis <- "proportion"
maturity_rows$fleet <- maturity_rows$survey <- maturity_rows$sex <- maturity_rows$region <- ""
maturity_rows$season <- "fall"
maturity_rows$year <- maturity_grid$year
maturity_rows$year_basis <- "calendar_year"
maturity_rows$age <- maturity_grid$age
maturity_rows$value <- as.numeric(maturity_grid$age >= 4)
maturity_rows$unit <- "proportion"
maturity_rows$source_type <- "translation_assumption"
maturity_rows$source_reference <- "DFO Research Document 2022/068, p. 10"
maturity_rows$transformation <- "Knife-edge maturity: ages 2-3 immature; ages 4-11+ mature."
maturity_rows$notes <- ""

fit_inputs <- rbind(catch_input, index_input, weight_rows, maturity_rows)
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  fit_inputs,
  years = years,
  ages = ages,
  sampling_times = c(
    "Fall fixed-gear CPUE - North" = 0.8,
    "Fall fixed-gear CPUE - Middle" = 0.8,
    "Fall fixed-gear CPUE - South" = 0.8,
    "4Tmno acoustic survey" = 0.8,
    "DFO multispecies RV survey" = 0.75
  ),
  surveys = c(
    "Fall fixed-gear CPUE - North",
    "Fall fixed-gear CPUE - Middle",
    "Fall fixed-gear CPUE - South",
    "4Tmno acoustic survey",
    "DFO multispecies RV survey"
  ),
  assumptions = source$assumptions
)

obs$maturity$spawn_time <- 212 / 365
obs$pop_covariates$M_process_center <- 0.2
obs$pop_covariates$age_blocks <- cut_ages(
  obs$pop_covariates$age, c(2, 7, 11)
)

is_cpue <- grepl("Fall fixed-gear CPUE", obs$index$survey, fixed = TRUE)
first_cpue_year <- min(obs$index$year[is_cpue])
obs$index$cpue_year <- ifelse(is_cpue, obs$index$year, first_cpue_year)
obs$index$cpue_active <- as.numeric(is_cpue)


## Background and comparisons ----
age_plus_group <- 11
comparison_outputs <- source$outputs[
  source$outputs$measure %in% c(
    "numbers_at_age", "fishing_mortality_at_age", "biomass_at_age",
    "mature_biomass_at_age", "recruitment", "Fbar", "total_numbers",
    "total_biomass", "SSB"
  ) &
    (is.na(source$outputs$region) | !nzchar(source$outputs$region)) &
    !(source$outputs$measure %in% c("total_numbers", "total_biomass") &
        !is.na(source$outputs$age_group)), ,
  drop = FALSE
]
comparison_outputs$age[comparison_outputs$measure == "SSB"] <- NA_integer_
comparison_outputs$age_group[comparison_outputs$measure == "SSB"] <- NA_character_
comparison_aggregates <- "ssb"

comparison_definitions <- list(
  biomass_at_age = list(
    status = "approximate",
    definition = "Official August 1 biomass at age compared with tinyAM annual biomass at age",
    reason = "The official output is August 1 biomass; tinyAM biomass-at-age is a start-of-year quantity."
  ),
  biomass = list(
    status = "approximate",
    definition = "Official August 1 total biomass compared with tinyAM annual biomass",
    reason = "The official output is August 1 biomass; tinyAM annual biomass uses its model reporting convention."
  ),
  ssb = list(
    status = "matched",
    definition = "August 1 mature biomass at ages 4-11+",
    reason = "Both use the ages 4+ maturity convention and August 1 timing; tinyAM pools three regional populations and uses approximate pooled weights."
  ),
  F_bar = list(
    status = "matched",
    definition = "January 1 abundance-weighted F at ages 5-10",
    reason = "The source and tinyAM calculations use ages 5-10; tinyAM represents the three source regions as one population."
  )
)

background <- c(
  "### Southern Gulf fall-spawning herring: detailed assessment to 2023",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years | Population estimates cover 1978-2023; CPUE, acoustic and RV indices have shorter series. | Fit 1978-2023 with all available indices. | The population period follows the official assessment; indices enter only in observed years. |",
  "| Ages | Ages 2-11+, with age 11 as the plus group. | Retain ages 2-11+. | The model age range and plus group are preserved. |",
  "| N | Recruitment enters at age 2; older fish survive with cohort structure. | Random-walk recruitment, IID N deviations and exponential initial abundance. | tinyAM does not reproduce the source's initial-state integration and likelihood. |",
  "| F | Three regional populations with fishery selectivity and time-varying F. | One population with a correlated random walk in age-specific F. | Regional populations and source selectivity are pooled or simplified. |",
  "| M | Estimated log-M random walks for ages 2-6 and 7-11+, with a 0.2 initial mean and 0.075 increment SD in the documented methods. | Two age-blocked AR1 processes centered on 0.2 with estimated process variation. | Stationary AR1 dynamics do not reproduce the source random walks or fixed increment SD. |",
  "| Catch | Fixed- and mobile-gear catch-at-age by region. The 2024 mobile table repeats fixed-gear values. | Sum fixed-gear catch-at-age across regions; the suspect mobile table is excluded. | This is a documented subset of removals, not a reconstructed total catch. |",
  "| Indices | Regional fixed-gear CPUE, fall acoustic age index, and September multispecies RV age index. | Preserve three regional CPUE series plus acoustic and RV series with separate q and SD levels. | The three regional CPUE series index one pooled tinyAM population. |",
  "| Timing | CPUE and acoustic observations are collected in fall; RV sampling is in September. | Use approximate fractions 0.8 for CPUE/acoustic and 0.75 for RV; SSB is timed to August 1 (212/365). | Exact annual survey dates are not supplied in the source tables. |",
  "| Weights and maturity | Beginning-year weights derive from fixed/mobile commercial weights; ages 4+ are mature. | Average available gear-region weights, interpolate gaps, apply the adjacent age-year geometric mean, and use knife-edge maturity at age 4. | Regional/gear pooling and edge-year weights are approximate. |",
  "| Comparisons | January 1 N and F, August 1 biomass and SSB, recruitment at age 2. | Compare N, F, recruitment, August 1 SSB, and Fbar; biomass-at-age and total biomass are marked approximate because of timing. | No numerical age-specific M output or source uncertainty is available. |",
  "",
  "The official assessment fits North, Middle, and South populations separately. tinyAM uses one pooled population. Catch includes only fixed-gear observations because the 2024 mobile-gear table duplicates fixed-gear catch values. The regional CPUE series remain separate and receive separate catchability levels, but they are all linked to the pooled population.",
  "",
  "The acoustic age-index values are in thousand fish, as clarified by DFO Research Document 2016/060, Table 16. CPUE and RV observations retain their number-per-net-haul and mean-number-per-tow scales. The observation likelihoods, source selectivity, and M process are simplified.",
  "",
  "Weights are pooled from available regional and gear-specific input rows because the source model uses regional weight surfaces and tinyAM has one population. The official fitted August 1 biomass-at-age is not directly comparable with tinyAM's annual biomass-at-age; August 1 SSB is compared using the source's ages 4+ maturity convention.",
  "",
  "The 2026 assessment is recorded separately as summary-only. Its revised acoustic-biomass series is not mixed into this 2024 translation."
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
    rec_settings = list(process = "rw"),
    F_settings = list(process = "cor_rw", mean_ages = 5:10),
    M_settings = list(
      process = "ar1", mu_form = NULL,
      mu_supplied = ~M_process_center,
      age_breaks = c(2, 7, 11)
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~0 + survey + mono(age, by = survey) +
        rw(cpue_year, by = cpue_active),
      q_link = "log",
      sd_form = ~0 + survey
    ),
    silent = silent
  )
}
