source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/read_database.R")

database <- read_database()
assessment_id <- "ices_whiting_west_scotland_2025"
assessment <- database$assessments[
  database$assessments$assessment_id == assessment_id, , drop = FALSE
]
inputs <- database$inputs[
  database$inputs$assessment_id == assessment_id, , drop = FALSE
]
outputs <- database$outputs[
  database$outputs$assessment_id == assessment_id, , drop = FALSE
]
series <- function(x, type, measure) {
  x[x$type == type & x$measure == measure, , drop = FALSE
  ]
}
key <- function(x) do.call(paste, c(lapply(x[c(
  "type", "measure", "fleet", "survey", "sex", "region", "season",
  "year", "age", "age_group"
)], function(value) {
  value <- as.character(value)
  value[is.na(value)] <- ""
  value
}), sep = "\034"))

catch <- series(inputs, "catch", "numbers_at_age")
index <- series(inputs, "index", "numbers_at_age")
weights <- series(inputs, "index", "relative_precision_weight")
n_age <- series(outputs, "population", "numbers_at_age")
f_age <- series(outputs, "mortality", "fishing_mortality_at_age")
m_age <- series(outputs, "mortality", "natural_mortality_at_age")
pred_catch <- series(outputs, "catch", "predicted_catch")
pred_index <- series(outputs, "index", "predicted_index")
q <- series(outputs, "catchability", "q")

stopifnot(
  nrow(assessment) == 1L,
  isTRUE(as.logical(assessment$is_current)),
  isTRUE(as.logical(assessment$is_applied)),
  assessment$inputs_status == "complete",
  assessment$outputs_status == "complete",
  assessment$terminal_year == 2025L,
  assessment$estimate_terminal_year == 2025L,
  nrow(catch) == 327L,
  identical(as.integer(range(as.integer(catch$year))), c(1981L, 2024L)),
  setequal(unique(as.integer(catch$age)), 0:7),
  all(catch$unit == "thousand fish"),
  nrow(index) == 472L,
  sum(index$survey == "WCIBTS.Q1") == 156L,
  sum(index$survey == "SCO.Q1") == 84L,
  sum(index$survey == "SWC.Q4") == 232L,
  identical(sort(unique(index$survey)), c("SCO.Q1", "SWC.Q4", "WCIBTS.Q1")),
  all(index$unit == "native survey index units"),
  all(is.finite(index$sampling_time)),
  all(weights$value > 0),
  nrow(weights) == 316L,
  nrow(n_age) == 360L,
  nrow(f_age) == 360L,
  nrow(m_age) == 360L,
  all(is.finite(n_age$value)), all(is.finite(n_age$se)),
  all(is.finite(f_age$value)), all(is.finite(f_age$se)),
  all(is.na(m_age$se)),
  nrow(pred_catch) == 327L,
  nrow(pred_index) == 472L,
  nrow(q) == 20L,
  !anyDuplicated(key(outputs)),
  isTRUE(all.equal(
    n_age$value[n_age$year == 2025L & n_age$age == 0L],
    385196.6, tolerance = 1e-5
  )),
  isTRUE(all.equal(
    outputs$value[outputs$measure == "SSB" & outputs$year == 2025L],
    19371.91, tolerance = 1e-5
  )),
  isTRUE(all.equal(
    outputs$value[outputs$measure == "Fbar" & outputs$year == 2025L],
    0.03022427, tolerance = 1e-6
  ))
)

cat("West of Scotland whiting source inputs and outputs passed validation.\n")
