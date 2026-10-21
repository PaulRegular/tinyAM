source("analysis/comp_assessments/tests/helper_stock.R")
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/read_database.R")

database <- read_database()
assessment_id <- "ices_whiting_north_sea_2026"
assessment <- database$assessments[
  database$assessments$assessment_id == assessment_id, , drop = FALSE
]
outputs <- database$outputs[
  database$outputs$assessment_id == assessment_id, , drop = FALSE
]

state_measures <- c(
  "numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age"
)
for (measure in state_measures) {
  x <- outputs[outputs$measure == measure, , drop = FALSE]
  stopifnot(
    nrow(x) == 49L * 9L,
    all(is.finite(x$value)), all(is.finite(x$se)),
    all(is.finite(x$lwr)), all(is.finite(x$upr)),
    all(x$lwr < x$value & x$value < x$upr),
    max(abs(log(x$lwr / x$value) + stats::qnorm(.975) * x$se / x$value)) < 5e-7,
    max(abs(log(x$upr / x$value) - stats::qnorm(.975) * x$se / x$value)) < 5e-7
  )
}

q <- outputs[outputs$measure == "q", , drop = FALSE]
catch_predictions <- outputs[outputs$measure == "predicted_catch", , drop = FALSE]
index_predictions <- outputs[outputs$measure == "predicted_index", , drop = FALSE]
output_key <- do.call(paste, c(lapply(outputs[c(
  "type", "measure", "fleet", "survey", "sex", "region", "season",
  "year", "age", "age_group"
)], function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value
}), sep = "\034"))
stopifnot(
  identical(assessment$outputs_status, "complete"),
  !anyDuplicated(output_key),
  nrow(q) == 13L,
  all(is.finite(q$value)), all(is.finite(q$se)),
  all(q$lwr < q$value & q$value < q$upr),
  nrow(catch_predictions) == 432L,
  nrow(index_predictions) == 509L,
  sum(index_predictions$survey == "IBTS-Q1") == 264L,
  sum(index_predictions$survey == "IBTS-Q3") == 245L,
  all(catch_predictions$unit == "thousand fish"),
  all(index_predictions$unit == "native survey-index units"),
  all(vapply(c("recruitment", "Fbar", "total_biomass"), function(measure) {
    any(outputs$measure == measure & outputs$year == 2026L &
      is.finite(outputs$lwr) & is.finite(outputs$upr))
  }, logical(1)))
)

cat("North Sea whiting native output coverage and uncertainty scales passed.\n")
