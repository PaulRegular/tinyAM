pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")

x <- run_assessment("dfo_american_plaice_4t_2016", fit = FALSE, silent = TRUE)
model <- .test_stock(x$source)
index <- x$obs$index
rv <- index[index$survey == "DFO Research Vessel Trawl Survey", , drop = FALSE]
sentinel <- index[index$survey == "Mobile Gear Sentinel Survey", , drop = FALSE]

stopifnot(
  identical(x$diagnostics$status, "not_fitted"),
  identical(x$translation_tier, "curated"),
  is.null(x$fit),
  identical(model$dat$F_settings$process, "cor_rw"),
  identical(model$dat$M_settings$process, "rw"),
  identical(model$dat$rec$type, "rw"),
  identical(as.integer(sort(unique(x$obs$catch$year))), 1976:2015),
  identical(as.integer(sort(unique(x$obs$weight$age))), 4:16),
  identical(as.integer(sort(unique(rv$year))), setdiff(1976:2015, 2003)),
  identical(as.integer(sort(unique(sentinel$year))), 2003:2015),
  nrow(rv) == 39L * 13L,
  nrow(sentinel) == 13L * 13L,
  all(rv$samp_time == 0.75),
  all(sentinel$samp_time == 0.625),
  length(unique(index$q_key)) == 26L,
  length(unique(x$obs$pop_covariates$M_assumption)) == 2L,
  x$obs$pop_covariates$M_assumption[
    x$obs$pop_covariates$age == 4
  ][1] == 0.6,
  x$obs$pop_covariates$M_assumption[
    x$obs$pop_covariates$age == 10
  ][1] == 0.3,
  x$obs$weight$obs[x$obs$weight$year == 2003 & x$obs$weight$age == 16] == 0,
  x$comparison$comparison_definitions$ssb$status == "matched",
  identical(x$comparison$comparison_age_groups$F[["4-9"]], 4:9),
  identical(x$comparison$comparison_age_groups$F[["10+"]], 10:16),
  identical(x$comparison$comparison_age_groups$M[["4-9"]], 4:9),
  identical(x$comparison$comparison_age_groups$M[["10+"]], 10:16)
)

cat("Curated 4T American plaice translation prepares both survey series and the mortality processes without fitting.\n")
