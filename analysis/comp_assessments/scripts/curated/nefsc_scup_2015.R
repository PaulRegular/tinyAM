## Observations ----
years <- 1984:2014
ages <- 0:7
surveys <- c("NEFSC fall trawl survey", "CTDEEP fall trawl survey")

source_inputs <- source$inputs
ctdeep <- source_inputs$type %in% "index" & source_inputs$survey %in% surveys[2L]
# The native 10+ bin contributes to the modeled 7+ index. Use the age series only.
ctdeep_plus <- ctdeep & source_inputs$measure %in% "index_by_age_group" &
  source_inputs$age %in% 10L
source_inputs$measure[ctdeep_plus] <- "numbers_at_age"
source_inputs <- source_inputs[!(ctdeep & source_inputs$measure %in% "total_numbers"), ]

catch_age <- source_inputs[
  source_inputs$type == "catch" & source_inputs$measure == "numbers_at_age",
  , drop = FALSE
]
if (!all(catch_age$unit %in% "thousand fish")) {
  cli::cli_abort("Scup source catch-at-age must be in thousand fish.")
}
catch_age$value <- as.numeric(catch_age$value)
catch_age$value[catch_age$fleet == "Recreational discards"] <-
  catch_age$value[catch_age$fleet == "Recreational discards"] * 0.15
# The shared converter prepares the biological inputs and indices with one catch series.
# Replace its pooled catch table below with the four fitted fleet series.
catch_total <- aggregate(value ~ year + age, catch_age, sum)
catch_rows <- catch_age[rep(1L, nrow(catch_total)), , drop = FALSE]
catch_rows$year <- catch_total$year
catch_rows$age <- catch_total$age
catch_rows$value <- catch_total$value
catch_rows$fleet <- "Total dead removals"
catch_rows$unit <- "thousand fish"
catch_rows$source_type <- "reconstructed_source_input"
catch_rows$source_reference <- paste(
  "NEFSC Reference Document 15-08, Tables A13, A16, A18 and A20, pp. 100, 103, 105, 107"
)
catch_rows$transformation <- paste(
  "Sum the four source catch-at-age series; multiply recreational discards by their reported 15% mortality rate."
)
catch_rows$notes <- "Analysis-only total removals series; commercial discards are fully mortal and recreational discards have 15% mortality."

source_weights <- source_inputs[
  source_inputs$type == "weight" & source_inputs$measure == "weight_at_age",
  , drop = FALSE
]
catch_for_weight <- source_inputs[
  source_inputs$type == "catch" & source_inputs$measure == "numbers_at_age" &
    as.numeric(source_inputs$age) >= 7,
  , drop = FALSE
]
catch_for_weight$value <- as.numeric(catch_for_weight$value)
weight_key <- paste(source_weights$year, source_weights$age)
catch_for_weight$weight <- as.numeric(source_weights$value[
  match(paste(catch_for_weight$year, catch_for_weight$age), weight_key)
])
if (anyNA(catch_for_weight$weight)) stop("Missing source weight for plus-group reconstruction.")
catch_for_weight$weighted <- catch_for_weight$value * catch_for_weight$weight
plus_weight <- aggregate(cbind(value, weighted) ~ year, catch_for_weight, sum)
plus_weight$weight <- plus_weight$weighted / plus_weight$value
if (any(!is.finite(plus_weight$weight))) {
  available <- is.finite(plus_weight$weight)
  plus_weight$weight[!available] <- stats::approx(
    plus_weight$year[available], plus_weight$weight[available],
    xout = plus_weight$year[!available], rule = 1
  )$y
}
if (any(!is.finite(plus_weight$weight))) stop("Could not reconstruct a complete 7+ weight series.")

fit_inputs <- source_inputs[
  !(source_inputs$type == "catch" & source_inputs$measure == "numbers_at_age") &
    !(source_inputs$type == "index" & source_inputs$measure == "numbers_at_age" &
      as.numeric(source_inputs$value) <= 0),
  , drop = FALSE
]
fit_inputs <- rbind(fit_inputs, catch_rows)

obs <- database_to_tam_obs(
  source$assessment$assessment_id[[1L]], fit_inputs,
  years = years, ages = ages, surveys = surveys,
  sampling_times = setNames(c(0.8, 0.75), surveys),
  assumptions = source$assumptions
)

obs$weight$obs[obs$weight$age == 7] <- plus_weight$weight[
  match(obs$weight$year[obs$weight$age == 7], plus_weight$year)
]
obs$maturity$spawn_time <- 0.5
catch_age$age <- pmin(as.integer(catch_age$age), max(ages))
fleet_catch <- aggregate(value ~ year + age + fleet, catch_age, sum)
if (nrow(fleet_catch) != length(years) * length(ages) * 4L ||
    anyDuplicated(fleet_catch[c("year", "age", "fleet")])) {
  stop("The four source catch fleets require complete, unique year-age grids.")
}
obs$catch <- data.frame(
  year = as.integer(fleet_catch$year), age = as.integer(fleet_catch$age),
  fleet = fleet_catch$fleet, obs = fleet_catch$value, unit = "thousand fish"
)
translation <- attr(obs, "translation")
translation$catch_method <- paste(
  "Keep numbers-at-age separate by source fleet, sum ages 7–10+ into age 7+,",
  "and apply 15% mortality to recreational discards."
)
translation$source_provenance <- rbind(
  .translation_source_provenance(catch_age, "catch", translation$catch_method),
  translation$source_provenance[translation$source_provenance$component != "catch", ]
)
attr(obs, "translation") <- translation
obs$catch$age_factor <- factor(obs$catch$age, levels = ages)
obs$index$q_key <- interaction(obs$index$survey, factor(obs$index$age, levels = ages),
  drop = TRUE, lex.order = TRUE
)

source_age_series <- function(measure, name) {
  x <- source$outputs[
    !is.na(source$outputs$measure) & source$outputs$measure == measure,
    c("year", "age", "age_group", "value"), drop = FALSE
  ]
  x$year <- as.integer(x$year)
  x$age <- suppressWarnings(as.integer(x$age))
  plus <- !is.na(x$age_group) & grepl("^[0-9]+\\+$", x$age_group)
  x$age[plus] <- as.integer(sub("\\+$", "", x$age_group[plus]))
  x$value <- as.numeric(x$value)
  names(x)[names(x) == "value"] <- name
  x[c("year", "age", name)]
}

source_n <- source_age_series("numbers_at_age", "N")
source_f <- source_age_series("fishing_mortality_at_age", "F")
source_m_values <- source_inputs[
  source_inputs$type == "M" & source_inputs$measure == "natural_mortality_at_age",
  , drop = FALSE
]
source_m <- data.frame(
  year = rep(years, each = length(ages)),
  age = rep(ages, times = length(years)),
  M = rep(as.numeric(source_m_values$value[match(ages, as.integer(source_m_values$age))]),
    times = length(years))
)
source_w <- obs$weight[c("year", "age", "obs")]
names(source_w)[3L] <- "weight"
source_p <- obs$maturity[c("year", "age", "obs")]
names(source_p)[3L] <- "maturity"
source_p$spawn_time <- obs$maturity$spawn_time
source_biology <- Reduce(function(x, y) merge(x, y, by = c("year", "age")),
  list(source_n, source_f, source_m, source_w, source_p)
)
source_biology$biomass <- source_biology$N * source_biology$weight
source_biology$ssb <- source_biology$biomass * source_biology$maturity *
  exp(-(source_biology$F + source_biology$M) * source_biology$spawn_time)

derived_comparison_rows <- function(measure, value, notes) {
  rows <- source$outputs[rep(NA_integer_, nrow(source_biology)), , drop = FALSE]
  rows$assessment_id <- source$assessment$assessment_id[[1L]]
  rows$type <- "biomass"
  rows$measure <- measure
  rows$year <- source_biology$year
  rows$age <- as.character(source_biology$age)
  rows$age_group <- ifelse(source_biology$age == max(ages), paste0(max(ages), "+"), "")
  rows$value <- value
  rows$se <- rows$lwr <- rows$upr <- NA_real_
  rows$unit <- "metric tons"
  rows$source_type <- "derived_source_output"
  rows$source_reference <- paste(
    "Derived from NEFSC Reference Document 15-08, Tables A24, A55-A57,",
    "MAT3, and fixed M assumption"
  )
  rows$notes <- notes
  rows
}

comparison_outputs <- rbind(
  source$outputs,
  derived_comparison_rows(
    "biomass_at_age", source_biology$biomass,
    "Common-definition total biomass by age: January 1 source N times translated annual weight; age 7 uses the translated 7+ weight."
  ),
  derived_comparison_rows(
    "mature_biomass_at_age", source_biology$ssb,
    "Common-definition SSB by age: source N, translated weight and maturity, with source F and M survival to the June spawning midpoint."
  )
)

## Background and comparisons ----
age_plus_group <- 7
comparison_definitions <- list(
  M = list(
    status = "matched",
    definition = "Fixed M = 0.20 per year, checked against the supplied assessment input",
    reason = "This is an input consistency check, not an estimated model comparison."
  ),
  ssb = list(
    status = "approximate",
    definition = "Common-definition SSB by age using source N and F/M with translated weight and maturity inputs",
    reason = "The translated 7+ weight and June spawning midpoint are approximations; the reported aggregate official SSB remains available as a native output."
  ),
  biomass = list(
    status = "approximate",
    definition = "Total biomass calculated from source N-at-age and translated weights",
    reason = "Age 7+ weight is reconstructed from the source ages 7-10 and fishery catch numbers."
  ),
  biomass_at_age = list(
    status = "approximate",
    definition = "Age-specific biomass calculated from source N-at-age and translated weights",
    reason = "Age 7+ weight is reconstructed from the source ages 7-10 and fishery catch numbers."
  )
)
background <- c(
  "### Scup: official 2015 SAW 60 ASAP assessment", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|---|---|---|",
  "| Years and ages | S60_BASE_18 modeled 1963-2014, true ages 0-7+. Detailed catch-at-age starts in 1984. | Fit 1984-2014, ages 0-7 with age 7 as the plus group. | Restricts the comparison to the years with reported catch-at-age, annual weights, and age-specific population outputs. |",
  "| Catch | Four source fleets; commercial discards are fully mortal and recreational discards have 15% mortality. | Retain the four source series after applying the reported recreational-discard mortality. | tinyAM fits lognormal numbers at age with a common SD. The source uses fleet-specific observation and composition likelihoods. Zero catch cells are retained as source values. |",
  "| NEFSC fall index | Stratified mean numbers per tow at ages 0–8. HBB values from 2009 onward are calibrated to ALB equivalents. The official model fits total indices and age compositions. | Fit positive age-specific numbers per tow at 0.8 of the year, with ages 7 and 8 combined into age 7+. | The age-specific lognormal likelihood approximates the official model's separate total-index and composition likelihoods. |",
  "| CTDEEP fall index | Table A43 reports Long Island Sound stratified mean numbers per tow, ages 0–10+, 1984–2014, with no 2010 survey. | Fit the age series with ages 7–10+ pooled to 7+, a separate q at each age, and one survey log-SD. | The 0.75 sampling fraction approximates the September–October midpoint. The published all-age total is used to check transcription. |",
  "| N, recruitment, and F | ASAP estimates annual abundance at age, age-0 recruitment, and fishing mortality with cohort survival between years. | Deterministic N survival, IID recruitment around a constant mean, and an AR1 F process with fleet-specific age means. | F deviations are independent across fleets with shared process parameters. The source uses time-blocked selectivity. |",
  "| M | Fixed M = 0.20 per year for all ages and years. | Use the same fixed M surface. | Fixed M is checked for input consistency, not treated as an estimated model match. |",
  "| Weights and maturity | Annual catch-number-weighted weights at ages 0-10; annual sex-combined MAT3 maturity. Spawning occurs May-August and peaks in June. | Retain annual weights and maturity; calculate age-7+ weight by weighting ages 7-10 with total source catch numbers; use a June midpoint at 0.5 of the year. | The same total-fishery weight table is used for all fleets. The plus-group weight is an explicit translation calculation. In 1999 no catch occurred at ages 7-10, so the missing plus-group weight is linearly interpolated between adjacent years for this fit only. |",
  "| Biomass, SSB, and recruitment | Annual aggregate SSB and age-0 recruitment are reported; 2014 SSB includes a 90% MCMC interval. | Derive source biomass-at-age from N and translated weights, and source SSB-at-age from N, translated weights and maturity with official F/M survival to the June spawning midpoint. | Biomass and SSB comparisons are approximate because the 7+ weight and spawning calculation are reconstructed; the directly reported aggregate SSB and its available uncertainty remain in the source record. |",
  "",
  "The fit approximates the 2015 S60_BASE_18 assessment with two of its survey series and age-specific lognormal observation likelihoods. The historical comparison uses the 2015 source record."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "iid", mu_form = ~1),
    F_settings = list(process = "ar1", mu_form = ~0 + fleet:age_factor, mean_ages = 3:7),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey),
    silent = silent
  )
}
