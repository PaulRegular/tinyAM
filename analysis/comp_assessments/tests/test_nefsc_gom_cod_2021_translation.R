pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

id <- "nefsc_gom_cod_2021_m02"
result <- run_assessment(id, fit = FALSE, silent = TRUE)
source_before <- result$source
stock <- .assessment_stock(source_before, do_fit = FALSE, silent = TRUE)
stopifnot(
  is.null(stock$error),
  is.null(stock$env$fit),
  identical(as.integer(stock$env$years), 1982:2019),
  identical(as.integer(stock$env$ages), 1:9),
  tinyAM::check_obs(result$obs),
  nrow(result$obs$catch) == 38L * 9L,
  nrow(result$obs$weight) == 38L * 9L,
  nrow(result$obs$maturity) == 38L * 9L,
  setequal(unique(result$obs$index$survey), c("NEFSC spring", "NEFSC fall")),
  setequal(unique(result$obs$index$samp_time), c(0.25, 0.75)),
  all(abs(result$obs$catch$relative_sd - sqrt(log1p(0.05^2))) < 1e-12),
  all(result$obs$pop_covariates$M_assumption == 0.2),
  identical(stock$env$comparison_outputs$measure, c(
    rep("SSB", 8L), rep("recruitment", 8L)
  )),
  identical(source_before$outputs, result$source$outputs),
  is.null(result$fit),
  result$diagnostics$status == "not_fitted"
)
settings <- stock$env$model_settings
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = result$obs, years = 1982:2019, ages = 1:9), settings
))
stopifnot(
  nrow(dat$obs$catch) == nrow(result$obs$catch),
  identical(settings$M_settings$process, "off"),
  identical(settings$N_settings$init, "exp"),
  identical(settings$F_settings$process, "ar1")
)
cat("NEFSC 2021 Gulf of Maine cod tinyAM observations, model preparation, and fit-disabled guard passed.\n")




