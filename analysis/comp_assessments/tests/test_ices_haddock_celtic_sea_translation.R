source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")

assessment_id <- "ices_haddock_celtic_sea_7b_k_2026"
source_record <- read_assessment(assessment_id)
translated <- .test_stock(source_record)
obs <- translated$obs
dat <- translated$dat

catch_weights <- source_record$inputs[
  source_record$inputs$type == "catch" &
    source_record$inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
weight_key <- function(x) paste(x$year, x$age, sep = "\r")
weight_index <- match(weight_key(obs$catch), weight_key(catch_weights))

tinyAM::check_obs(obs)
stopifnot(
  identical(translated$years, 1993:2025),
  identical(translated$ages, 0:8),
  as.integer(translated$age_plus_group) == 8L,
  nrow(obs$catch) == 33L * 9L,
  sum(obs$catch$obs == 0, na.rm = TRUE) == 35L,
  nrow(obs$index) == 23L * 8L - 1L,
  all(obs$index$samp_time == 0.875),
  !any(obs$index$year == 2005L & obs$index$age == 7L),
  all(is.finite(obs$index$obs) & obs$index$obs > 0),
  all(obs$catch$relative_sd == 1 / sqrt(catch_weights$value[weight_index])),
  setequal(unique(obs$index$sd_age), c("ages0_5", "ages6_7")),
  dat$N_settings$process == "iid",
  dat$N_settings$init == "exp",
  dat$rec_settings$process == "rw",
  dat$F_settings$process == "cor_rw",
  identical(dat$F_settings$mean_ages, 3:5),
  dat$M_settings$process == "off",
  length(unique(obs$index$q_key)) == 8L,
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 0] == 1.09),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$age == 8] == 0.36)
)

cat("ICES Celtic Sea haddock translation preparation passed.\n")
