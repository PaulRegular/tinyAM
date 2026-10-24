pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

database <- read_database()
source_before <- read_assessment("nefsc_white_hake_2025", database)
result <- run_assessment("nefsc_white_hake_2025", database, fit = FALSE, silent = TRUE)
obs <- result$obs
prepared <- .assessment_stock(source_before, do_fit = FALSE, silent = TRUE)
stopifnot(is.null(prepared$error))
settings <- prepared$env$model_settings
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = 1963:2024, ages = 1:9), settings
))
starts <- prepared$env$source_starts()
reported_r <- subset(source_before$outputs, measure == "recruitment")
reported_r <- as.numeric(reported_r$value[match(1963:2024, reported_r$year)])
stopifnot(
  identical(all.vars(settings$rec_settings$mu_form), character()),
  identical(settings$rec_settings$process, "iid"),
  !any(grepl("sr_alpha|sr_beta", names(tinyAM::make_par(dat)))),
  isTRUE(all.equal(unname(exp(starts$log_r0)), reported_r[1L])),
  isTRUE(all.equal(unname(exp(starts$log_r)), reported_r[-1L])),
  !"sd_supplied" %in% names(settings$index_settings),
  is.null(prepared$env$fit)
)
prepared$env$source$outputs$unit[
  prepared$env$source$outputs$measure == "recruitment"
] <- "fish"
testthat::expect_error(prepared$env$source_starts(), "thousand fish")

stopifnot(
  result$diagnostics$status == "not_fitted",
  inherits(obs, "list"),
  tinyAM::check_obs(obs),
  isTRUE(all.equal(as.numeric(sort(unique(obs$weight$year))), 1963:2024)),
  isTRUE(all.equal(as.numeric(sort(unique(obs$weight$age))), 1:9)),
  nrow(obs$weight) == 62L * 9L,
  nrow(obs$maturity) == 62L * 9L,
  nrow(obs$catch) == 62L * 9L,
  setequal(unique(obs$catch$year[is.finite(obs$catch$obs)]), 1989:2024),
  nrow(obs$index) == 170L * 9L - 10L,
  all(tapply(obs$index$obs, obs$index$q_key, function(x) any(is.finite(x) & x > 0))),
  setequal(unique(obs$index$survey), c(
    "NEFSC spring bottom trawl survey",
    "NEFSC fall bottom trawl survey",
    "ASMFC shrimp trawl survey",
    "NEFSC spring bottom longline survey",
    "NEFSC fall bottom longline survey"
  )),
  setequal(unique(obs$index$samp_time), c(0.25, 0.5, 0.75)),
  all(obs$pop_covariates$M_assumption == 0.2),
  identical(source_before$inputs, result$source$inputs),
  is.null(result$fit),
  is.null(result$dashboard_file),
  grepl("aggregate-index and composition likelihoods", paste(result$background, collapse = " "))
)

cat("White Hake tinyAM observations, source preservation, and fit-disabled preparation passed.\n")
