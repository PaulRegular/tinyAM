## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 1983:2025
ages <- 1:4
inputs <- source$inputs

catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age" &
  inputs$age %in% ages, ]
annual_catch <- stats::aggregate(value ~ year + age, catch, sum)
annual_rows <- function(rows, values) {
  out <- rows[match(paste(values$year, values$age), paste(rows$year, rows$age)), ]
  out$value <- values$value
  out$season <- "annual"
  out
}
catch_rows <- annual_rows(catch, annual_catch)
catch_rows$transformation <- "Analysis only: sum first- and second-half catch numbers at the same year and age."

m <- inputs[inputs$type == "M" & inputs$age %in% ages, ]
m_rows <- annual_rows(m, stats::aggregate(value ~ year + age, m, sum))
m_rows$basis <- "per_year"
m_rows$unit <- "per year"
m_rows$transformation <- "Analysis only: add the two integrated half-year M hazards at the same year and age."

biology <- inputs[inputs$type %in% c("weight", "maturity") &
  inputs$season == "1" & inputs$age %in% ages, ]
biology$season <- "annual"
index <- inputs[inputs$type == "index" & inputs$age %in% ages, ]
# These are scaled indices, not absolute fish counts to normalize to thousands.
index$unit <- ifelse(index$survey == "Dredge", "native survey index (unit unresolved)",
  "source index scale per effort unit (RTM; absolute multiplier not stated)")
translated_inputs <- rbind(catch_rows, m_rows, biology, index)
obs <- database_to_tam_obs(assessment_id, translated_inputs,
  years = years, ages = ages, sampling_times = c(Dredge = .9375, RTM = 0),
  surveys = c("Dredge", "RTM"), assumptions = source$assumptions)
obs$maturity$spawn_time <- 0

effort <- inputs[inputs$type == "covariate" & inputs$measure == "effort" &
  inputs$season == "annual", ]
stopifnot(!anyDuplicated(effort$year), setequal(effort$year, years), all(effort$value > 0))
pop <- obs$pop_covariates
pop$log_effort <- log(effort$value[match(pop$year, effort$year)]) - mean(log(effort$value))
period <- cut(pop$year, c(1982, 1998, 2021, 2025), labels = c("early", "middle", "late"))
pop$F_mean_group <- interaction(period, pmin(pop$age, 3), drop = TRUE)
obs$pop_covariates <- pop
obs$catch$catch_sd_group <- factor(ifelse(obs$catch$age <= 2, "young", "old"))
obs$index$index_sd_group <- factor(ifelse(obs$index$survey == "Dredge", "Dredge_1",
  ifelse(obs$index$age == 1, "RTM_1", "RTM_2_3")))

cw <- inputs[inputs$type == "catch_weight" & inputs$age %in% ages, ]
key <- function(x) paste(x$year, x$age, x$season)
catch$weight <- cw$value[match(key(catch), key(cw))] / 1000
catch$biomass_kg <- catch$value * 1e6 * catch$weight
annual_mass <- stats::aggregate(biomass_kg ~ year + age, catch, sum)
annual_mass$count <- annual_catch$value[match(paste(annual_mass$year, annual_mass$age),
  paste(annual_catch$year, annual_catch$age))] * 1e6
stopifnot(all(is.finite(annual_mass$biomass_kg)), all(annual_mass$count > 0))
catch_reporting <- list(
  weights = data.frame(year = annual_mass$year, age = annual_mass$age,
    weight = annual_mass$biomass_kg / annual_mass$count),
  totals = stats::setNames(stats::aggregate(biomass_kg ~ year, annual_mass, sum), c("year", "yield")))

## Background and comparisons ----
age_plus_group <- 4
comparison_outputs <- source$outputs[!(source$outputs$measure == "fishing_mortality_at_age" &
  source$outputs$season != "annual"), ]
comparison_aggregates <- c("abundance", "biomass", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", definition = "January abundance at ages 1-4+.",
    reason = "An annual model approximates the two-season assessment; July age-0 states are omitted."),
  F = list(status = "approximate", definition = "Annual integrated F at ages 1-4+.",
    reason = "The annual F mean uses effort and source age/period groups, with independent residual states."),
  F_bar = list(status = "approximate", definition = "January-N-weighted annual F over ages 1-2.",
    reason = "Recalculated from official annual F and January N, rather than comparing native arithmetic Fbar."),
  M = list(status = "matched", definition = "Sum of the supplied integrated half-year M hazards at each age.",
    reason = "Fixed-M input consistency only; M is not estimated in either fit."),
  M_bar = list(status = "approximate", definition = "January-N-weighted fixed annual M over ages 1-4+.",
    reason = "Any difference comes from N weights, not M estimation. This is not evidence of estimated-M agreement."),
  recruitment = list(status = "approximate", definition = "January age-1 abundance in the same calendar year.",
    reason = "Official January N at age 1 replaces native July age-0 recruitment for this comparison; no cohort-year shift."),
  ssb = list(status = "approximate", definition = "January SSB over ages 1-4+ under shared first-half biology.",
    reason = "Same timing and biology, but annual rather than seasonal population dynamics and rounded source inputs."),
  abundance = list(status = "approximate", definition = "January total abundance over ages 1-4+.",
    reason = "The later July age-0 recruits are not added to January abundance."),
  biomass = list(status = "approximate", definition = "January biomass over ages 1-4+ under shared first-half weights.",
    reason = "Reported weights are rounded; source native aggregates are retained separately."))
background <- c(
  "### Sandeel Area 1r: official 2026 SMS-effort assessment", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | 1983-2025, two half-years, age 0 enters in July, January ages 1-4+. | Calendar-year model, ages 1-4+. | tinyAM has annual dynamics. July age-0 catch is excluded from the source likelihood; Dredge age 0 is omitted here. |",
  "| N and recruitment | Free initial January abundance, deterministic older-cohort survival; weak hockey-stick penalty on July recruitment, fixed breakpoint 105,809 t and weight 0.05. | Free initial abundance, no older-age N process, IID age-1 entry about an estimated mean. | January entrants are a different recruitment state. The weighted hockey-stick model is not reproduced. |",
  "| F | Proportional to seasonal standardized effort; age/season coefficients change in 1999 and 2022; ages 3-4+ share selectivity. | IID F about age-group-by-period means and a fitted log annual-effort effect. | The complete published annual effort is a proxy; the seasonal proportional relationship and exact latent-state sharing are not reproduced. |",
  "| M | Fixed supplied half-year integrated hazards. | Sum the two seasonal hazards at each year and age, with no M process. | Annual cumulative mortality is preserved; within-year exposure is simplified. Fixed M is an input check. |",
  "| Catch | Positive half-year numbers at ages 1-4+, with age/season SD groups. | Sum half-year counts, convert million to thousand fish, estimate SDs for ages 1-2 and 3-4+. | An annual lognormal observation is not the product of two seasonal likelihoods. Reported eligible catch cells number 320 versus 318 in the model diagnostic; exact selection flags are unavailable. |",
  "| Surveys | Dredge ages 0-1 averages abundance over annual times 0.875-1; RTM ages 1-3 uses configured time 0. | Dredge age 1 at midpoint 0.9375; RTM ages 1-3 at time 0. Separate survey-age q; RTM ages 2-3 share SD. | Midpoint timing approximates the Dredge integral. Original scaled index values are preserved; q is not restricted to one. |",
  "| Biology and reporting | Seasonal weights/maturity; January SSB. | First-half biology, spawning time zero; catch totals use observed-catch-weighted seasonal catch weights. | Catch reporting is ages 1-4+ only, rather than SAG all-age landings. No source estimates alter observations. |",
  "| Comparisons | Native N0 is July recruitment, older N is January; annual F sums both seasons. | Compare January ages 1-4+, same-year age-1 recruitment and common-definition aggregates. | Native July recruitment and native aggregates remain source context. Derived aggregates have no invented uncertainty. |",
  "", "This is a limited annual approximation, not a reproduction of the seasonal model. Two seasonal effort rows, exact catch flags and unrounded native arrays are unavailable. The precise production software version is also unverified. The fit starts with package defaults; one aligned source-start retry is permitted only after numerical failure."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "iid", mu_form = ~0 + F_mean_group + log_effort, mean_ages = 1:2),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + catch_sd_group),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + index_sd_group),
  silent = silent, ...)
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages)
    f <- .translation_start_surface(comparison_outputs, "mortality", "fishing_mortality_at_age", years, ages)
    list(log_r0 = log(n[1L, 1L]), log_r = log(n[-1L, 1L]),
      log_n0 = log(n[1L, -1L]), log_f = log(f))
  }, attempt_file = file.path("analysis/comp_assessments/results/cache", assessment_id,
    "initialization_attempts.rds"))
}
