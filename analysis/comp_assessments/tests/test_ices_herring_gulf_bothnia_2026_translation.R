if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run the translation check.", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

result <- run_assessment(
  "ices_herring_gulf_bothnia_2026", fit = FALSE, silent = TRUE
)
obs <- result$obs

stopifnot(
  result$diagnostics$status == "not_fitted",
  !is.null(obs),
  nrow(obs$catch) == 46L * 15L,
  setequal(unique(obs$catch$year), 1980:2025),
  setequal(unique(obs$catch$age), 1:15),
  all(is.finite(obs$catch$obs)),
  nrow(obs$index) == 19L * 15L,
  setequal(unique(obs$index$year), 2007:2025),
  setequal(unique(obs$index$age), 1:15),
  all(obs$index$survey == "BIAS acoustic survey"),
  all(obs$index$samp_time == 0.75),
  all(is.finite(obs$index$obs) & obs$index$obs > 0),
  nrow(obs$weight) == 46L * 15L,
  all(obs$weight$unit == "kg"),
  nrow(obs$maturity) == 46L * 15L,
  all(obs$maturity$spawn_time == 0),
  all(obs$pop_covariates$M_assumption > 0),
  setequal(unique(obs$pop_covariates$age), 1:15)
)

cat("Gulf of Bothnia herring observation conversion and preparation passed.\n")
