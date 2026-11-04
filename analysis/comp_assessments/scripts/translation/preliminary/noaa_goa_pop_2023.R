## Observations ----
years <- 1990:2023
ages <- 2:25
survey <- "GOA bottom trawl survey"
inputs <- source$inputs
inputs <- inputs[inputs$measure != "proportion_at_length", ]
# Use the model population weights; preserve supporting survey/fishery weights separately.
population_weight <- inputs$type == "weight" & (is.na(inputs$survey) | inputs$survey == "")
population_weights <- inputs[population_weight, ]
maturity <- source$outputs[source$outputs$type == "biology" &
  source$outputs$measure == "maturity_at_age" & source$outputs$year == 1961, ]
template <- inputs[which(population_weight)[1], ]
maturity_rows <- template[rep(1, nrow(maturity)), ]
maturity_rows$type <- "maturity"
maturity_rows$measure <- "maturity_at_age"
maturity_rows$basis <- "proportion"
maturity_rows$age <- maturity$age
maturity_rows$value <- maturity$value
maturity_rows$unit <- "proportion"
maturity_rows$source_reference <- maturity$source_reference
maturity_rows$transformation <- "Official fitted female maturity ogive used as a fixed biological plug-in; its uncertainty is not propagated."
inputs <- rbind(inputs, maturity_rows)
catch_weight <- inputs$type == "catch_weight"
inputs$season[catch_weight] <- "annual"
inputs <- rbind(inputs[!catch_weight, ], .translation_repeat_biology(
  inputs[catch_weight, ], years, NULL, "Catch weights"))
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, weight_survey = "", maturity_multiplier = 0.5,
  sampling_times = stats::setNames(0, survey), assumptions = source$assumptions)
# Preserve ages 26–29+ for tinyAM's internal plus-group biology.
older_weight <- .translation_repeat_biology(population_weights[population_weights$age > max(ages), ],
  years, NULL, "Older population weights")
if (nrow(older_weight)) {
  obs$weight <- rbind(obs$weight, data.frame(year = older_weight$year,
    age = older_weight$age, obs = older_weight$value / 1000, unit = "kg"))
  older_maturity <- .translation_repeat_biology(maturity_rows[maturity_rows$age > max(ages), ],
    years, NULL, "Older maturity")
  obs$maturity <- rbind(obs$maturity, data.frame(year = older_maturity$year,
    age = older_maturity$age, obs = older_maturity$value * 0.5))
  older <- older_weight[c("year", "age")]
  obs$catch <- rbind(obs$catch, transform(older, obs = NA_real_, unit = "thousand fish"))
  obs$pop_covariates <- rbind(obs$pop_covariates, older)
}
N_settings <- list(process = "off", init = "exp")
rec_settings <- list(process = "iid", mu_form = ~1)
F_settings <- list(process = "off", mu_form = ~logistic(age) + factor(year), mean_ages = 5:15)
M_settings <- list(process = "off", mu_form = ~1, mu_supplied = NULL)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~logistic(age), sd_form = ~1)
age_plus_group <- 25L
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", reason = "Observed age compositions include age error; tinyAM treats their bins as population ages."),
  ssb = list(status = "approximate", reason = "Fitted maturity is a plug-in and plus-group biology is approximated; native SSB changes timing in 2023."))

## Background ----
background <- c(
  "### Gulf of Alaska Pacific ocean perch: preliminary comparison with Official Model 20.1", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Coverage | 1961–2023, population ages 2–29+; observed age bins end at 25+. | 1990–2023, model ages 2–25+ with older source biology retained. | Only years with matched age compositions supply reconstructed observations. The first catch-composition year is 1990. |",
  "| N and recruitment | Deterministic survival with penalized recruitment deviations and equilibrium initial ages. | Deterministic survival, IID recruitment about a mean and exponential initial abundance. | Initial mortality and recruitment penalty details differ. |",
  "| F and q | Parametric selectivity, with changing fishery blocks and survey selectivity. | Static logistic F and q age curves, with annual F effects. | Time-varying selectivity and source shape restrictions are simplified. |",
  "| M | One estimated constant with a lognormal prior. | Estimate one constant without that prior. | No reported M estimate is silently fixed as an input. |",
  "| Catch and index | Biomass totals, age compositions with age-error matrices and sample weights, and fishery length data. | Reconstruct number-at-age only where totals, number compositions and matching weights overlap. Static fishery weights apply to the annual catch stream. | Independent lognormal age errors omit age error, shared-total dependence, source sample weights and length likelihoods. Using the age-25 weight for the observed 25+ bin is approximate. |",
  "| Survey timing | Native prediction is q times beginning-year N, selectivity and weights, without within-year survival. | Use effective fraction 0. | This follows the source equation, not the calendar date of the trawl survey. |",
  "| Biology | Model weights in g/fish; estimated female maturity with a 0.5 female fraction. | Convert weights to kg and use the reported maturity ogive as a fixed plug-in, multiplied by 0.5. | Maturity uncertainty is not propagated. Supporting survey and fishery weights are used only for the corresponding conversions. |",
  "| SSB | Beginning-year SSB through 2022, survival to one-third year in 2023. | Common-definition beginning-year mature biomass. | Native SSB is retained separately; equivalence is not claimed. |",
  "", "Default initialization is attempted first. Official N may initialize one numerical retry after explicit summation of ages 25–29+. F is formula-based and retains package-default starts. No observations are derived from fitted N or F."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages, N_settings = N_settings,
    rec_settings = rec_settings, F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent)
  fit <- fit_default_first(fit_model, retry_start = function() {
    dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages,
      N_settings = N_settings, rec_settings = rec_settings, F_settings = F_settings,
      M_settings = M_settings, catch_settings = catch_settings, index_settings = index_settings)
    states <- source
    n <- states$outputs[states$outputs$measure == "numbers_at_age", ]
    n$age <- pmin(n$age, max(ages))
    pooled <- aggregate(value ~ year + age, n, sum)
    n <- n[!duplicated(paste(n$year, n$age)), ]
    n$value <- pooled$value[match(paste(n$year, n$age), paste(pooled$year, pooled$age))]
    # F is formula-based here; only N initializes, and all coefficients retain defaults.
    states$outputs <- n
    source_start_parameters(states, dat)
  }, attempt_file = file.path(.assessment_cache(source$assessment$assessment_id,
    translation_tier), "initialization_attempts.rds"))
}
