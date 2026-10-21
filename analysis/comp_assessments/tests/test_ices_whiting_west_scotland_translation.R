source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")
source("analysis/comp_assessments/R/read_database.R")

database <- read_database()
assessment_id <- "ices_whiting_west_scotland_2025"
source_data <- read_assessment(assessment_id, database)
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat
index <- obs$index
weights <- source_data$inputs[
  source_data$inputs$type == "index" &
    source_data$inputs$measure == "relative_precision_weight",
  , drop = FALSE
]
key <- function(x) paste(x$survey, x$year, x$age, sep = "\034")
weight_key <- key(weights)
matched <- match(key(index), weight_key)
has_weight <- !is.na(matched)

stopifnot(
  tinyAM::check_obs(obs),
  identical(translated$years, 1981:2025),
  identical(translated$ages, 0:7),
  identical(translated$age_plus_group, 7L),
  nrow(obs$catch) == 360L,
  sum(is.finite(obs$catch$obs)) == 327L,
  identical(levels(obs$catch$sd_block), c("age0", "ages1_6", "age7_plus")),
  nrow(index) == 472L,
  identical(as.integer(table(factor(
    index$survey, levels = c("WCIBTS.Q1", "SCO.Q1", "SWC.Q4")
  ))), c(156L, 84L, 232L)),
  all(index$samp_time[index$survey == "WCIBTS.Q1"] == 0.125),
  all(index$samp_time[index$survey == "SCO.Q1"] == 0.125),
  isTRUE(all.equal(
    unique(index$samp_time[index$survey == "SWC.Q4"]),
    0.8530996, tolerance = 1e-8
  )),
  length(unique(index$q_key)) == 17L,
  !anyNA(index$relative_sd),
  all(is.finite(index$relative_sd) & index$relative_sd > 0),
  isTRUE(all.equal(
    index$relative_precision_weight[has_weight],
    as.numeric(weights$value[matched[has_weight]])
  )),
  isTRUE(all.equal(
    index$relative_sd[has_weight],
    1 / sqrt(as.numeric(weights$value[matched[has_weight]]))
  )),
  all(index$relative_sd[!has_weight] == 1),
  all(obs$pop_covariates$M_assumption > 0),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$rec_settings$process, "rw"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$F_settings$mean_ages, 1:3),
  identical(dat$M_settings$process, "off"),
  is.null(translated$start_par),
  !"start_par" %in% names(translated$fit_call),
  is.null(translated$fit)
)

cat("West of Scotland whiting translation structure passed validation.\n")
