pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

x <- run_assessment("dfo_haddock_4x5y_2016", fit = FALSE, silent = TRUE)
catch <- x$obs$catch
index <- x$obs$index
ages_11plus <- x$source$inputs[
  x$source$inputs$type == "catch" & x$source$inputs$year == 2015 &
    x$source$inputs$age >= 11, "value"
]
catch_11plus <- catch$obs[catch$year == 2015 & catch$age == 11]

stopifnot(
  identical(x$diagnostics$status, "not_fitted"),
  is.null(x$fit),
  all(sort(unique(x$obs$weight$age)) == 1:11),
  length(unique(x$obs$weight$year)) == 31L,
  all(sort(unique(index$age)) == 1:10),
  all(sort(unique(index$year)) == 1985:2015),
  nrow(index) == 31L * 10L - 1L,
  all(index$samp_time == 0.5),
  length(unique(x$obs$pop_covariates$M_assumption)) == 4L,
  length(catch_11plus) == 1L,
  catch_11plus == sum(ages_11plus),
  identical(x$comparison$comparison_definitions$ssb$status, "approximate")
)

cat("4X5Y haddock translation prepares the intended model data without fitting.\n")
