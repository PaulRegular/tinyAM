source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))

assessment_id <- "ices_norway_pout_north_sea_2026_benchmark"
record <- read_assessment(assessment_id, read_database())
fit_env <- new.env(parent = baseenv())
load(file.path(
  root, "source_cache", assessment_id, "native", "sum.RData"
), fit_env)
fit <- fit_env$sesamsum
dat <- fit$data
outputs <- record$outputs
predictions <- outputs[
  outputs$measure %in% c("predicted_catch", "predicted_index"), , drop = FALSE
]
predictions$value <- as.numeric(predictions$value)
survey_names <- c(
  "IBTS Q1 (north of 57N)",
  "EGFS Q3",
  "SGFS Q3 (age-0 through 2012; ages 1-3+ thereafter)",
  "IBTS other countries Q3",
  "SGFS Q3 age 0 (2013 onward)"
)
expected <- data.frame(
  type = ifelse(dat$fleet == 1L, "catch", "index"),
  measure = ifelse(dat$fleet == 1L, "predicted_catch", "predicted_index"),
  fleet = ifelse(dat$fleet == 1L, "Total catch", ""),
  survey = ifelse(dat$fleet == 1L, "", survey_names[pmax(dat$fleet - 1L, 1L)]),
  season = paste0(
    "Q", as.integer(round((dat$t1 - floor(dat$t1 + 1e-8)) * 4) + 1)
  ),
  year = as.character(floor(dat$t1 + 1e-8)),
  age = as.character(dat$ageFrom),
  age_group = ifelse(
    dat$ageFrom == dat$ageTo,
    ifelse(dat$ageFrom == 3L, "3+", ""),
    paste0(dat$ageFrom, "-", dat$ageTo, "+")
  ),
  value = exp(as.numeric(fit$rep$predObs)),
  stringsAsFactors = FALSE
)
key <- function(x) {
  fields <- c("type", "measure", "fleet", "survey", "season", "year", "age", "age_group")
  values <- lapply(x[fields], function(value) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    value
  })
  do.call(paste, c(values, sep = "\r"))
}
actual <- predictions[c(
  "type", "measure", "fleet", "survey", "season", "year", "age", "age_group", "value"
)]
actual$year <- as.character(actual$year)
actual$age <- as.character(actual$age)
matched <- match(key(expected), key(actual))

stopifnot(
  fit$opt$convergence == 0L,
  length(fit$rep$predObs) == 1092L,
  nrow(predictions) == 1092L,
  sum(predictions$measure == "predicted_catch") == 668L,
  sum(predictions$measure == "predicted_index") == 424L,
  !anyDuplicated(key(actual)),
  !anyNA(matched),
  isTRUE(all.equal(actual$value[matched], expected$value, tolerance = 1e-12)),
  all(predictions$source_type == "native_model"),
  all(is.na(predictions$se)),
  all(is.na(predictions$lwr)),
  all(is.na(predictions$upr)),
  all(predictions$unit[predictions$measure == "predicted_catch"] == "million fish"),
  all(predictions$unit[predictions$measure == "predicted_index"] == "native survey-index units"),
  sum(predictions$age_group == "2-3+", na.rm = TRUE) == 26L,
  all(grepl("fit\\$rep\\$predObs", predictions$source_reference)),
  all(grepl("no bias correction", predictions$notes, fixed = TRUE))
)

cat("ICES Norway pout native prediction records passed.\n")
