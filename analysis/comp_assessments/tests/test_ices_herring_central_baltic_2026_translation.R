pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

result <- run_assessment("ices_herring_central_baltic_2026", fit = FALSE)
obs <- result$obs
m <- result$source$inputs[
  result$source$inputs$type == "M" &
    result$source$inputs$measure == "fixed_mortality_medium",
]
key <- function(x) paste(x$year, x$age, sep = "\034")
m_index <- match(key(obs$pop_covariates), key(m))

stopifnot(
  result$diagnostics$status == "not_fitted",
  nrow(obs$catch) == 52L * 9L,
  all(is.na(obs$catch$obs[obs$catch$age == 0])),
  nrow(obs$index) == 20L * 8L,
  identical(as.integer(sort(unique(obs$index$age))), 1:8),
  all(obs$index$survey == "BIAS"),
  all(obs$index$samp_time == 0.8),
  all(as.character(obs$index$q_key) == as.character(obs$index$age)),
  !anyNA(m_index),
  isTRUE(all.equal(obs$pop_covariates$M_assumption, m$value[m_index])),
  nrow(obs$weight) == 52L * 9L,
  nrow(obs$maturity) == 52L * 9L,
  grepl("medium-M scenario", paste(result$background, collapse = " "))
)

cat("Central Baltic herring translation preparation and source M mapping passed.\n")