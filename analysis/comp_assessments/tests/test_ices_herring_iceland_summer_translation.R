root <- file.path("analysis", "comp_assessments")
.libPaths(c(file.path(root, "source_cache", "R-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "ices_herring_iceland_summer_2026"
source_data <- read_assessment(assessment_id, read_database())
stock_env <- new.env(parent = environment(.assessment_stock))
stock_env$source <- source_data
stock_env$do_fit <- FALSE
stock_env$silent <- TRUE
sys.source(file.path(
  root, "scripts", "translation", "stocks", paste0(assessment_id, ".R")
), envir = stock_env)

derived <- stock_env$comparison_outputs[
  stock_env$comparison_outputs$source_type == "derived_official_output", , drop = FALSE
]
stopifnot(
  identical(stock_env$years, 1980:2023), identical(stock_env$ages, 2:15),
  nrow(stock_env$obs$catch) == 616L,
  nrow(stock_env$obs$index) == 534L,
  sum(is.na(stock_env$obs$catch$obs)) == 4L,
  all(stock_env$obs$catch$obs[!is.na(stock_env$obs$catch$obs)] > 0),
  all(stock_env$obs$index$obs > 0),
  identical(levels(stock_env$obs$pop_covariates$F_age_block),
    c(paste0("age", 2:6), "age7plus")),
  identical(levels(stock_env$obs$index$q_age),
    c("age2to3", "age4", "age5", "age6", "age7plus")),
  identical(sort(unique(derived$measure)),
    c("biomass_at_age", "mature_biomass_at_age")),
  sum(derived$measure == "biomass_at_age") == 616L,
  sum(derived$measure == "mature_biomass_at_age") == 616L,
  all(is.na(derived$se)),
  !exists("start_par", envir = stock_env, inherits = FALSE),
  !is.null(stock_env$prepared)
)

cat("Icelandic summer herring translation preparation passed.\n")
