## Observations ----
years <- 1980:2025
ages <- 1:7
surveys <- c("UK BTS", "FR GFS")

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  assumptions = source$assumptions
)
obs$index$sd_block <- factor(obs$index$survey)
obs$index$q_key <- interaction(
  obs$index$survey,
  pmin(obs$index$age, 5L),
  drop = TRUE,
  sep = "."
)

N_settings <- list(process = "rw", init = "exp")
rec_settings <- list(process = "rw", mu_form = NULL)
F_settings <- list(process = "cor_rw", mu_form = NULL, mean_ages = 3:6)
M_settings <- list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)

dat <- tinyAM::prepare_tam(
  data = obs,
  years = years,
  ages = ages,
  N_settings = N_settings,
  rec_settings = rec_settings,
  F_settings = F_settings,
  M_settings = M_settings,
  catch_settings = catch_settings,
  index_settings = index_settings
)

comparison_aggregates <- "ssb"
comparison_definitions <- list(
  F_bar = list(
    status = "non_equivalent",
    definition = "ICES SAG unweighted Fbar for ages 3-6; tinyAM reports abundance-weighted F over ages 3-6.",
    reason = "Official N-at-age and F-at-age are unavailable, so an ICES-weighted mean cannot be reconstructed."
  ),
  ssb = list(
    status = "matched",
    definition = "Reported aggregate SSB at start of year for ages 1-7; both use the same stock weights and maturity with no pre-spawning mortality.",
    reason = ""
  )
)

background <- c(
  "### Eastern English Channel plaice: Official 2026 ICES assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM translation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | Fit through 2025; ages 1-7 with age 7 as the plus group; 2026 recruitment and SSB are forecasts. | Fit 1980-2025, ages 1-7. | The forecast year is shown as official context, not fitted as a historical year. |",
  "| Recruitment and N | Age-1 recruitment and cohort abundance follow the SAM process structure. | Use package random walks for recruitment and cohort survival with exponential N0. | tinyAM does not reproduce SAM's exact age-specific variance sharing or initialization. |",
  "| F and M | SAM has an age-correlated F random walk, shares the age-6/7 F state, and uses fixed age-specific M. | Use tinyAM correlated-RW F and the same fixed M input. | tinyAM estimates separate age states and only approximates the age-correlation structure. |",
  "| Catch and indices | One aggregate catch-at-age stream; age-specific UK BTS Q3 and FR GFS Q4 indices. | Keep native age-specific observations, use separate survey q patterns and SDs. | Index scaling is not stated; SAM uses survey-specific q and observation variance sharing. |",
  "| Biology and SSB | Annual stock/catch/landings/discard weights, fixed maturity, and start-of-year SSB timing. | Use the archived biology surfaces and zero pre-spawning mortality. | No biological surfaces are substituted. |",
  "| Comparison outputs | SAG reports recruitment, SSB and Fbar, with 95% intervals. | Compare age-1 recruitment and start-of-year SSB; show Fbar as non-equivalent. | Official N-at-age, F-at-age and total abundance are unavailable. The local pinned SAM refit did not reproduce the published SAG estimates, so its state estimates are not reported as official. |",
  "",
  "The fit first uses package-default starting values. If numerical checks fail, one retry uses N, F, and survey-q values from the cached local refit of the pinned SAM model as starting values only; tinyAM estimates all parameters freely. That local refit does not reproduce SAG's published aggregate estimates and is not treated as an official output source. Catch data and model inputs come from the public SAM archive; annual catch-at-age totals agree with SAG catches to rounding. Fbar is retained as an FYI comparison only: ICES reports an unweighted ages 3-6 mean, while tinyAM reports abundance-weighted F. The archive configuration also says fbarRange=3 7, so the published reporting definition is kept and the discrepancy is noted in the source review."
)

fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) {
    tinyAM::fit_tam(
      data = obs,
      years = years,
      ages = ages,
      N_settings = N_settings,
      rec_settings = rec_settings,
      F_settings = F_settings,
      M_settings = M_settings,
      catch_settings = catch_settings,
      index_settings = index_settings,
      silent = silent,
      start_par = start_par
    )
  }
  retry_start <- function() {
    source_fit_file <- file.path(
      "analysis/comp_assessments/source_cache/ices_plaice_7d_2026",
      "fit_local_locked.RData"
    )
    if (!file.exists(source_fit_file)) {
      cli::cli_abort("The cached, pinned SAM refit needed for the one retry is unavailable.")
    }
    source_fit_env <- new.env(parent = emptyenv())
    suppressWarnings(load(source_fit_file, envir = source_fit_env))
    source_fit <- source_fit_env[["fit"]]
    if (is.null(source_fit) || !is.list(source_fit[["pl"]])) {
      cli::cli_abort("The cached SAM refit does not contain fitted process states.")
    }

    log_N <- t(source_fit[["pl"]][["logN"]])
    log_F <- t(source_fit[["pl"]][["logF"]])
    if (!identical(dim(log_N), c(length(years), length(ages))) ||
      !identical(dim(log_F), c(length(years), length(ages) - 1L)) ||
      !identical(as.integer(source_fit[["data"]][["years"]]), years)) {
      cli::cli_abort("The cached SAM N/F state surfaces do not align with this translation.")
    }
    log_F <- cbind(log_F, log_F[, ncol(log_F)])

    starts <- tinyAM::make_par(dat)
    starts[["log_r0"]] <- log_N[1L, 1L]
    starts[["log_r"]] <- log_N[-1L, 1L]
    starts[["log_n"]] <- log_N[-1L, -1L, drop = FALSE]
    starts[["log_f"]] <- log_F

    index <- obs$index
    n_at_index <- exp(log_N[cbind(
      match(index$year, years),
      match(index$age, ages)
    )])
    q_start <- tapply(
      index$obs / n_at_index,
      as.character(index$q_key),
      stats::median,
      na.rm = TRUE
    )
    names(q_start) <- paste0("q_key", names(q_start))
    if (!setequal(names(q_start), names(starts[["log_q"]])) ||
      any(!is.finite(q_start)) || any(q_start <= 0)) {
      cli::cli_abort("Could not align positive survey-q starting values to the fitted groups.")
    }
    starts[["log_q"]][names(q_start)] <- log(q_start)
    starts
  }
  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      "analysis/comp_assessments/results/cache/ices_plaice_english_channel_7d_2026",
      "initialization_attempts.rds"
    )
  )
}
