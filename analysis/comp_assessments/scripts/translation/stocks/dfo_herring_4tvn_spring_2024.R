translate_stock <- function(source) {
  years <- 1978:2023
  ages <- 2:11
  inputs <- source$inputs

  catch <- inputs[
    inputs$type == "catch" & inputs$measure == "numbers_at_age" &
      inputs$season == "spring",
    , drop = FALSE
  ]
  catch <- aggregate(value ~ year + age, catch, sum)
  catch_rows <- data.frame(
    assessment_id = source$assessment$assessment_id,
    type = "catch", measure = "numbers_at_age", basis = "numbers",
    fleet = "", survey = "", sex = "", region = "", season = "spring",
    year = catch$year, year_basis = "calendar_year", age = catch$age,
    value = catch$value, unit = "thousand fish", sampling_time = NA_real_,
    source_type = "translation_assumption",
    source_reference = paste("DFO Research Document 2024/058, Tables 4 and 8"),
    transformation = "Summed fixed- and mobile-gear catch-at-age for each year and age.",
    notes = "Combined stream for tinyAM; the source assessment models its catch observations differently.",
    stringsAsFactors = FALSE
  )

  gear_weights <- inputs[
    inputs$type == "catch_weight" & inputs$measure == "weight_at_age" &
      inputs$season == "spring" & inputs$value > 0,
    , drop = FALSE
  ]
  combined_weights <- aggregate(value ~ year + age, gear_weights, function(x) {
    exp(mean(log(x)))
  })
  weight_grid <- expand.grid(year = years, age = ages)
  weight_grid <- merge(weight_grid, combined_weights,
                       by = c("year", "age"), all.x = TRUE, sort = FALSE)
  weight_grid <- weight_grid[order(weight_grid$year, weight_grid$age), ]
  for (age in ages) {
    rows <- weight_grid$age == age
    observed <- rows & is.finite(weight_grid$value)
    weight_grid$value[rows] <- stats::approx(
      weight_grid$year[observed], weight_grid$value[observed],
      xout = weight_grid$year[rows], rule = 2
    )$y
  }
  source_weight <- matrix(weight_grid$value, nrow = length(years), byrow = TRUE,
                          dimnames = list(year = years, age = ages))
  stock_weight <- source_weight
  for (i in 2:length(years)) {
    for (j in 2:length(ages)) {
      stock_weight[i, j] <- sqrt(source_weight[i - 1, j - 1] * source_weight[i, j])
    }
  }
  weight_rows <- expand.grid(year = years, age = ages)
  weight_rows <- data.frame(
    assessment_id = source$assessment$assessment_id,
    type = "weight", measure = "weight_at_age", basis = "kg_per_fish",
    fleet = "", survey = "", sex = "", region = "", season = "spring",
    year = weight_rows$year, year_basis = "calendar_year",
    age = weight_rows$age,
    value = as.vector(stock_weight), unit = "kg per fish",
    sampling_time = NA_real_, source_type = "translation_assumption",
    source_reference = "DFO Research Document 2024/058, Tables 5 and 9; DFO Research Document 2022/068, p. 6",
    transformation = paste(
      "Geometric mean of available fixed/mobile weights by year-age; missing",
      "values interpolated within age. Beginning-year age a uses the geometric",
      "mean of age a-1 in year t-1 and age a in year t; age 2 and 1978 use",
      "same-year weights because earlier ages/years are unavailable."
    ),
    notes = "Fit-only stock-weight approximation; the accepted fitted weight surface and gear-combination rule are not published.",
    stringsAsFactors = FALSE
  )

  add_composition_fields <- function(x) {
    x$observation_id <- ""
    x$length_bin <- NA_real_
    x$length_bin_lower <- NA_real_
    x$length_bin_upper <- NA_real_
    x$sample_size <- NA_real_
    x$age_error <- NA_real_
    x$partition <- NA_real_
    x[names(inputs)]
  }
  catch_rows <- add_composition_fields(catch_rows)
  weight_rows <- add_composition_fields(weight_rows)

  fit_inputs <- inputs[
    !(inputs$type == "catch" & inputs$measure == "numbers_at_age") &
      !(inputs$type == "weight" & inputs$measure == "weight_at_age"),
    , drop = FALSE
  ]
  fit_inputs <- rbind(fit_inputs, catch_rows, weight_rows)
  obs <- database_to_tam_obs(
    source$assessment$assessment_id,
    fit_inputs,
    years = years,
    ages = ages,
    sampling_times = c("Spring fixed-gear CPUE" = 0.25),
    surveys = "Spring fixed-gear CPUE",
    assumptions = source$assumptions
  )
  obs$weight$M_process_center <- 0.2

  settings <- list(
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "ar1", mu_form = NULL),
    M_settings = list(
      process = "ar1", mu_form = NULL, mu_supplied = ~ M_process_center,
      age_breaks = c(2, 7, 11), first_dev_year = 1978L
    ),
    catch_settings = list(sd_form = ~ 1, fill_missing = TRUE),
    index_settings = list(q_form = ~ 1, sd_form = ~ 1,
                          fill_missing = TRUE)
  )
  dat <- do.call(tinyAM::make_dat, c(
    list(obs = obs, years = years, ages = ages), settings
  ))
  start_par <- tinyAM::make_par(dat)
  source_surface <- function(measure, type, scale = 1) {
    rows <- source$outputs[
      source$outputs$measure == measure & source$outputs$type == type,
      , drop = FALSE
    ]
    surface <- matrix(NA_real_, length(years), length(ages),
                      dimnames = list(year = years, age = ages))
    surface[cbind(match(rows$year, years), match(rows$age, ages))] <-
      as.numeric(rows$value) * scale
    if (any(!is.finite(surface))) {
      cli::cli_abort("Source {measure} output does not cover the fitted year-age grid.")
    }
    surface
  }
  source_N <- source_surface("numbers_at_age", "population", 1000)
  source_F <- source_surface("fishing_mortality_at_age", "mortality")
  source_recruitment <- source$outputs[
    source$outputs$measure == "recruitment" &
      source$outputs$type == "recruitment", , drop = FALSE
  ]
  source_recruitment <- source_recruitment[
    match(years, source_recruitment$year), , drop = FALSE
  ]
  start_par$log_r0 <- log(source_recruitment$value[[1]] * 1000)
  start_par$log_r <- setNames(
    log(source_recruitment$value[-1] * 1000), as.character(years[-1])
  )
  start_par$log_f[] <- log(pmax(source_F, 1e-6))
  index_rows <- obs$index
  index_N <- source_N[cbind(match(index_rows$year, years), match(index_rows$age, ages))]
  index_F <- source_F[cbind(match(index_rows$year, years), match(index_rows$age, ages))]
  q_start <- stats::median(
    index_rows$obs / (index_N * exp(-(index_F + 0.2) * index_rows$samp_time)),
    na.rm = TRUE
  )
  start_par$log_q[] <- log(q_start)
  start_par$log_sd_r <- log(0.5)
  start_par$log_sd_f <- log(0.2)
  start_par$log_sd_m <- log(0.075)
  start_par$log_sd_catch[] <- log(0.2)
  start_par$log_sd_index[] <- log(0.2)

  list(
    years = years,
    ages = ages,
    age_plus_group = 11,
    obs = obs,
    comparison_outputs = source$outputs[
      source$outputs$measure %in% c(
        "numbers_at_age", "fishing_mortality_at_age", "biomass_at_age",
        "recruitment", "Fbar"
      ), , drop = FALSE
    ],
    comparison_scales = c(N = 1e-3, F = 1, recruitment = 1e-3,
                          biomass_at_age = 1e-3, F_bar = 1),
    start_par = start_par,
    settings = settings,
    background = c(
      "### Southern Gulf spring-spawning herring: detailed assessment to 2023",
      "",
      print_sources(source$assessment),
      "",

      "The current database record is DFO's accepted 2024 assessment. The 2026 assessment is recorded separately as summary-only because its detailed methods and numerical tables are still in preparation. The 2022 report supplies model-method context where the 2024 support document is silent.",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | The SCA reports 1978-2023 population estimates; CPUE is available for 1990-2021. | Fit 1978-2023 and retain the spring fixed-gear CPUE series at 1990-2021. | Catch and population outputs extend beyond the CPUE series. |",
      "| Ages | Ages 2-11+, with age 11 as the plus group. | Ages 2-11, with age 11 as the plus group. | The reported age range is retained. |",
      "| N | Age-2 recruitment varies annually; older first-year cohorts are reconstructed from recruitment and survival. | Exponential initial abundance with tinyAM's estimated annual recruitment variation. | tinyAM does not reproduce the source's initial-cohort estimation or fixed recruitment SD. |",
      "| F | Logistic fishery selectivity changes across three time blocks; the source estimates initial fishing mortality and observation/process terms. | Age- and year-correlated AR1 F states. | This is a compact approximation to changing selectivity; it is not the source's period-specific logistic model. |",
      "| M | Log-M follows random walks for ages 2-6 and 7-11+, with 0.2 initial-M prior means and increment SD fixed at 0.075. | Two M age blocks follow an AR1 process centered on 0.2. | tinyAM's random walk leaves its first M state unpenalized; AR1 supplies a mean-reverting initial-state distribution, but changes the source process and estimates its variation. |",
      "| Catch | Fixed- and mobile-gear catch-at-age are reported separately in thousand fish. The accepted SCA uses total catch and age-composition likelihoods. | Sum gear catches into one age-specific stream; zero ages are filled as latent observations for the lognormal likelihood. | tinyAM cannot reproduce the source total-plus-composition likelihood. |",
      "| Index | Spring fixed-gear CPUE and a fishery-independent acoustic survey are used as aggregate indices with age compositions; CPUE q follows a random walk. | Use published age-specific CPUE values at a spring timing of 0.25, with a constant q; omit acoustic sample counts. | The report does not provide compatible aggregate index values for tinyAM; acoustic sample counts are composition data, not abundance indices. |",
      "| Weights and maturity | The source uses a knife-edge maturity schedule at ages 3-4 and beginning-year weights derived from gear-specific catch weights. | Use the same maturity schedule and a fit-only weight surface derived from the published gear weights. | Gear combination for the source weight surface is not fully specified; see the explicit transformation in this recipe. |",
      "| Comparison | Numerical January 1 N, biomass-at-age and F-at-age are published, with age-2 recruitment and Fbar. | Compare these common quantities; the source's tabulated biomass is not SSB. | SSB and age-specific M are not available as numerical tables for this assessment. |",
      "",
      "The comparison is a structural approximation, not a reproduction of the accepted likelihood. Survey timing is represented by a seasonal midpoint, and the source's estimated M prior is not available as a tinyAM penalty."
    )
  )
}
