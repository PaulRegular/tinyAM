## Observations ----
years <- 1971:2018
ages <- 3:12 # Accepted model starts at age 2, but lack data on catch numbers for age 2; starting at age 3
rv <- "DFO September RV survey"
mobile <- "Mobile Sentinel August survey"

rv_weights <- source$inputs[
  source$inputs$type == "weight" &
    source$inputs$measure == "weight_at_age" &
    source$inputs$survey == rv &
    source$inputs$year %in% years &
    source$inputs$age %in% 2:11,
  , drop = FALSE
]
missing_years <- setdiff(years, unique(rv_weights$year))
if (!all(missing_years %in% c(1980, 1985))) {
  cli::cli_abort("Only the documented 1980 and 1985 RV weight gaps can be interpolated.")
}
interpolated_weights <- do.call(rbind, lapply(missing_years, function(year) {
  rows <- rv_weights[rv_weights$year == year - 1, , drop = FALSE]
  next_year <- rv_weights[rv_weights$year == year + 1, , drop = FALSE]
  if (!setequal(rows$age, 2:11) || !setequal(next_year$age, 2:11)) {
    cli::cli_abort("Adjacent RV weight years must contain ages 2-11 before interpolation.")
  }
  rows <- rows[match(2:11, rows$age), , drop = FALSE]
  next_year <- next_year[match(2:11, next_year$age), , drop = FALSE]
  rows$year <- year
  rows$value <- (rows$value + next_year$value) / 2
  rows$source_type <- "translation_assumption"
  rows$transformation <- "Linear interpolation of RV weights for this fit only."
  rows$notes <- "No official RV weight-at-age was tabulated for this year."
  rows
}))
rv_weights <- rbind(rv_weights, interpolated_weights)
plus_weight <- rv_weights[rv_weights$age == 11, , drop = FALSE]
plus_weight$age <- 12
plus_weight$source_type <- "translation_assumption"
plus_weight$transformation <- "Age-11 RV weight carried forward to the tinyAM 12+ group."
plus_weight$notes <- "Fit-only proxy; RV weights are not reported for ages 12+."

landings_at_age <- !is.na(source$inputs$type) &
  source$inputs$type == "catch" &
  !is.na(source$inputs$measure) &
  source$inputs$measure == "landings_numbers_at_age"
if (!any(landings_at_age)) {
  cli::cli_abort("The source-reported landings-at-age series is missing.")
}
catch_proxy <- source$inputs[landings_at_age, , drop = FALSE]
catch_proxy$measure <- "numbers_at_age"
catch_proxy$source_type <- "translation_assumption"
catch_proxy$transformation <- paste(
  "Used source landings-at-age as a fit-only proxy for catch-at-age;",
  "no scaling to total catch was possible because the fitted age",
  "composition was not tabulated."
)
catch_proxy$notes <- paste(
  "This is not the accepted model's catch-at-age input. The source",
  "model fits total catch biomass and age proportions separately."
)
inputs <- rbind(
  source$inputs[!landings_at_age, , drop = FALSE], catch_proxy,
  interpolated_weights, plus_weight
)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  inputs,
  years = years,
  ages = ages,
  weight_survey = rv,
  sampling_times = c(setNames(0.75, rv), setNames(0.625, mobile)),
  surveys = c(rv, mobile),
  exclude_index_years = setNames(list(c(1980, 1985, 2003)), rv),
  assumptions = source$assumptions
)
obs$weight$M_prior_mean <- ifelse(obs$weight$age <= 4, 0.65, 0.15)
comparison_outputs <- source$outputs[
  source$outputs$measure %in% c(
    "numbers_at_age", "fishing_mortality_at_age", "SSB", "recruitment"
  ), , drop = FALSE
]
terminal_m <- source$outputs[
  source$outputs$measure == "natural_mortality_at_age", , drop = FALSE
]
if (nrow(terminal_m)) {
  terminal_m <- do.call(rbind, lapply(seq_len(nrow(terminal_m)), function(i) {
    row <- terminal_m[i, , drop = FALSE]
    group_ages <- switch(as.character(row$age_group),
                         "5-8" = 5:8, "9+" = 9:12, integer())
    if (!length(group_ages)) return(NULL)
    row <- row[rep(1L, length(group_ages)), , drop = FALSE]
    row$age <- group_ages
    row$age_group <- NA_character_
    row$notes <- paste(row$notes,
                       "The reported group value is expanded across its ages for comparison.")
    row
  }))
  comparison_outputs <- rbind(comparison_outputs, terminal_m)
}
obs$catch$age_blocks <- cut_ages(
  obs$catch$age,
  c(3, 5, 7, 9, 12)
)


## Background and comparisons ----

age_plus_group <- 12


comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001)

background <- c("### Southern Gulf cod: accepted detailed assessment to 2018", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Ages | The accepted model uses ages 2-12+, with recruitment at age 2. | Fit ages 3-12+, with recruitment at age 3. | Age-2 catch numbers are unavailable in the main published landings series, so age 2 is excluded rather than fitted without catch observations. |",
    "| N | Recruitment enters at age 2 and depends on spawning biomass two years earlier, with autocorrelated variation. | Recruitment enters at age 3, with deterministic cohort survival and a recruitment random walk. | The recruitment age and underlying recruitment dynamics differ from the accepted model. |",
    "| F | The source estimates fully recruited F and period-specific logistic selectivity. | Use a year- and age-correlated AR1 F process with separate mean log F for age blocks 3-4, 5-6, 7-8, and 9-12+. | The blockwise mean and correlated deviations approximate changing fishing mortality without reproducing the source selectivity model. |",
    "| M | The source estimates log-M random walks for ages 2-4, 5-8, and 9+, with priors on initial M and fixed increment SD 0.075. | Use random walks for ages 3-4, 5-8, and 9-12+, with initial states informed by supplied M values and a warm-start fit. | Age 2 is excluded; tinyAM does not reproduce the source initial-M priors or fixed increment SD. |",
    "| Catch | The source fits annual catch biomass and proportions-at-age for ages 2-12+. Published landings numbers-at-age cover ages 3-12+. | Fit published landings numbers-at-age for ages 3-12+ as a proxy for catch-at-age. | The fitted source catch compositions are unavailable, and the tinyAM likelihood differs from the source biomass-plus-composition likelihood. |",
    "| Index | The source uses RV, mobile sentinel, and longline indices, including observations at age 2. | Use RV and mobile sentinel indices, restricted to modeled ages 3-12+, with approximate sampling times 0.75 and 0.625. Exclude RV 1980, 1985, and 2003. | Age 2 and the longline index are omitted; survey timing and index reconstruction are approximations. |",
    "| Comparison | The report provides SSB, numbers-at-age, fishing mortality, and age-2 recruitment estimates. | Compare SSB, numbers-at-age for ages 3+, F, and reported terminal M groups. | Direct recruitment comparisons are inappropriate because tinyAM recruitment is defined at age 3 rather than age 2. |")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_stage <- "warm_start"
  warm_fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~0 + age_blocks),
    M_settings = list(process = "iid", mu_form = NULL, mu_supplied = ~M_prior_mean, age_breaks = c(3,
        5, 9, 12), first_dev_year = 1971L),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey, fill_missing = FALSE),
    silent = silent
  )
  if (!isTRUE(warm_fit$is_converged)) cli::cli_abort("The preliminary fit did not converge.")
  start_par <- as.list(warm_fit$sdrep, "Estimate")
  fit_stage <- "fit"
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "ar1", mu_form = ~0 + age_blocks),
    M_settings = list(process = "rw", mu_form = NULL, mu_supplied = ~M_prior_mean, age_breaks = c(3,
        5, 9, 12), first_dev_year = 1971L),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
