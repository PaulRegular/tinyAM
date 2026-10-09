translate_stock <- function(source) {
  years <- 1978:2023
  ages <- 2:11
  inputs <- source$inputs

  catch <- inputs[
    inputs$type == "catch" & inputs$measure == "numbers_at_age" &
      inputs$season == "spring",
    , drop = FALSE
  ]
  gear_catch <- catch[c("year", "age", "fleet", "value")]
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
  gear_weights <- merge(
    gear_weights[c("year", "age", "fleet", "value")], gear_catch,
    by = c("year", "age", "fleet"), suffixes = c("_weight", "_catch")
  )
  combined_weights <- do.call(rbind, lapply(
    split(gear_weights, interaction(gear_weights$year, gear_weights$age)),
    function(x) {
      if (!nrow(x)) return(NULL)
      value <- if (sum(x$value_catch) > 0) {
        stats::weighted.mean(x$value_weight, x$value_catch)
      } else {
        mean(x$value_weight)
      }
      data.frame(year = x$year[1], age = x$age[1], value = value)
    }
  ))
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
      "Catch-number-weighted arithmetic mean of available fixed/mobile weights",
      "by year-age (unweighted mean if both catches are zero); missing",
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
    sampling_times = c(
      "Spring fixed-gear CPUE" = 0.25,
      "4Tmno acoustic survey" = 0.75
    ),
    surveys = c("Spring fixed-gear CPUE", "4Tmno acoustic survey"),
    assumptions = source$assumptions
  )
  obs$weight$M_process_center <- 0.2
  obs$catch$age_blocks <- cut_ages(
    obs$catch$age,
    c(2, 4, 6, 8, 11)
  )

  is_cpue <- obs$index$survey == "Spring fixed-gear CPUE"

  cpue_years <- sort(unique(obs$index$year[is_cpue]))
  cpue_breaks <- floor(seq(min(cpue_years), max(cpue_years), length.out = 4))

  obs$index$cpue_period <- "baseline"
  obs$index$cpue_period[is_cpue] <- as.character(cut_years(obs$index$year[is_cpue], cpue_breaks))
  obs$index$cpue_period <- factor(obs$index$cpue_period) |> relevel(ref = "baseline")

  settings <- list(
    N_settings = list(
      process = "iid",
      init = "exp"
    ),
    F_settings = list(
      process = "ar1",
      mu_form = ~ 0 + age_blocks,
      mean_ages = 6:8
    ),
    M_settings = list(
      process = "ar1",
      mu_form = NULL,
      mu_supplied = ~ M_process_center,
      age_breaks = c(2, 7, 11),
      first_dev_year = 1978L
    ),
    catch_settings = list(
      sd_form = ~ 1,
      fill_missing = FALSE
    ),
    index_settings = list(
      q_form = ~ 0 + mono(age, by = survey) + cpue_period,
      q_link = "log",
      sd_form = ~ 0 + survey,
      fill_missing = FALSE
    )
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
  start_par$log_n[] <- log(source_N[-1, -1])
  index_rows <- dat$obs$index
  index_N <- source_N[cbind(match(index_rows$year, years), match(index_rows$age, ages))]
  index_F <- source_F[cbind(match(index_rows$year, years), match(index_rows$age, ages))]
  q_start <- index_rows$obs /
    (index_N * exp(-(index_F + 0.2) * index_rows$samp_time))
  positive <- is.finite(q_start) & q_start > 0
  start_par$log_q[] <- stats::lm.fit(
    dat$q_modmat[positive, , drop = FALSE],
    log(q_start[positive]) - drop(
      dat$q_mono_modmat[positive, , drop = FALSE] %*% start_par$dq
    )
  )$coefficients
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
        "recruitment", "Fbar", "total_numbers", "total_biomass",
        "SSB", "mature_biomass_at_age"
      ) & !(source$outputs$measure %in% c("total_numbers", "total_biomass") &
              !is.na(source$outputs$age_group)), , drop = FALSE
    ],
    comparison_scales = c(N = 1e-3, F = 1, recruitment = 1e-3,
                          abundance = 1e-3, biomass = 1e-3,
                          biomass_at_age = 1e-3, ssb = 1e-3,
                          ssb_mat = 1e-3, F_bar = 1),
    comparison_definitions = list(ssb = list(
      status = "matched",
      definition = "January 1 mature biomass at ages 4-11+, derived from reported biomass-at-age and maturity",
      reason = "This common-definition comparison is not the accepted assessment's April 1 SSB. Source uncertainty is unavailable."
    )),
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
      "| N | Age-2 recruitment varies annually, with initial-cohort deviations; older cohorts subsequently follow survival. | Exponential initial abundance, annual recruitment variation and IID process error for older ages. | The older-age N process is retained for a stable fit; it adds flexibility absent from the source survival model. Recruitment variation is estimated rather than fixed. |",
      "| F | Logistic fishery selectivity changes across three time blocks; the source estimates initial fishing mortality and observation/process terms. | Age- and year-correlated AR1 F states. | This is a compact approximation to changing selectivity; it is not the source's period-specific logistic model. |",
      "| M | Log-M follows random walks for ages 2-6 and 7-11+, with 0.2 initial-M prior means and increment SD fixed at 0.075. | Two M age blocks follow an AR1 process centered on 0.2. | tinyAM's random walk leaves its first M state unpenalized; AR1 supplies a mean-reverting initial-state distribution, but changes the source process and estimates its variation. |",
      "| Catch | Fixed- and mobile-gear catch-at-age are reported separately in thousand fish. The accepted SCA uses total catch and age-composition likelihoods. | Sum gear catches into one age-specific stream and convert to fish. Zeros and missing observations are excluded rather than filled. | tinyAM cannot reproduce the source total-plus-composition likelihood. |",
      "| Index | Spring fixed-gear CPUE and the acoustic survey use age compositions and age-aggregated biomass indices; the acoustic biomass likelihood was weighted by 3 in the 2022 method report. | Use CPUE ages 4-11 and acoustic ages 2-10 at approximate timings of 0.25 and 0.75, with separate non-decreasing age-q curves and separate observation SDs. CPUE q changes in three time blocks. | tinyAM fits age-specific lognormal indices rather than composition plus aggregate biomass. The 2022 biomass likelihoods use CPUE ages 4-10 and acoustic ages 4-8. CPUE time blocks approximate the source q random walk. |",
      "| Index scale | Table 15 omits the unit; the same historical acoustic series is explicitly labeled thousands of fish in DFO 2016/060, Table 16. | Record the table values as thousand fish in the database; the shared converter multiplies by 1,000 to obtain fish. Keep the log-q link. | Unit clarification corrects the acoustic q scale without introducing a manual stock-script multiplier or changing process settings. |",
      "| Weights and maturity | Maturity is knife-edge between ages 3 and 4. Beginning-year weights use combined fishery weights and a geometric mean across adjacent age-year cells. | Combine gear weights with a catch-number-weighted mean, then apply the adjacent-cell geometric mean and the same maturity schedule. | Gear weighting is inferred: its median difference from weights implied by published biomass/N is 0.073%. Missing cells, age 2 and 1978 still require the documented approximation. |",
      "| Comparison | January 1 N, biomass-at-age and F-at-age, age-2 recruitment and Fbar are tabulated. | Also show total abundance, total biomass and January 1 mature biomass derived from the source age tables. | Derived source totals have no uncertainty. The SSB panel compares January 1 mature biomass, not the assessment's April 1 SSB. No numerical M surface is available. |",
      "",
      "The comparison is a structural approximation, not a reproduction of the accepted likelihood. Acoustic values are converted from thousand fish to fish using the database units; the 2025 biomass download belongs to the separate 2026 record. Acoustic timing is approximated by 0.75 from the late September to early October survey window. The source's initial-M prior, fixed process SDs and acoustic likelihood weight are not reproduced.",
      "",
      "The current fit retains IID N, exponential initial abundance and AR1 M. Trials with source-like deterministic older-age survival or free/random initial abundance did not converge from the tested starts. An M random walk and smoother CPUE q converged but worsened agreement with the source outputs. The revised gear-weight calculation improves the biological representation; it does not remove the remaining differences in abundance, recruitment or mortality."
    )
  )
}
