root <- file.path("analysis", "comp_assessments")
inputs <- read.csv(file.path(root, "database", "inputs.csv"),
                   stringsAsFactors = FALSE)
inputs <- inputs[inputs$assessment_id == "dfo_cod_3pn4rs_2025", ]
index <- inputs[inputs$type == "index", ]
weights <- inputs[inputs$type == "catch_weight", ]

stopifnot(
  nrow(index) == 30L * 11L,
  setequal(index$year, 1995:2024),
  setequal(index$age, 1:11),
  all(table(index$year, index$age) == 1L),
  all(index$survey == "Sentinel mobile"),
  all(index$unit == "mean numbers per tow"),
  all(index$source_reference == "DFO Research Document 2026/010 Table 28"),
  all(index$sampling_time == 0.54),
  all(grepl("approximation", index$notes)),
  index$value[index$year == 1997 & index$age == 1] == 0,
  index$value[index$year == 2024 & index$age == 11] == 0.087,
  all(grepl("11\\+", index$notes[index$age == 11])),
  nrow(weights) == 503L,
  all(weights$unit == "kg"),
  all(weights$basis == "kg_per_fish"),
  all(weights$source_reference == "DFO Research Document 2026/010 Table 24"),
  weights$value[weights$year == 1974 & weights$age == 2] == 0,
  weights$value[weights$year == 2024 & weights$age == 11] == 3.41,
  !any(weights$year == 2006 & weights$age %in% 2:3),
  all(is.finite(inputs$value)),
  !any(inputs$type %in% c("weight", "maturity"))
)

# Catch weights and raw maturity samples must not make this stock fit-ready.
source(file.path(root, "R", "run_assessment.R"))
source_data <- read_assessment("dfo_cod_3pn4rs_2025", read_database())
before <- source_data
error <- tryCatch(database_to_tam_obs(
  "dfo_cod_3pn4rs_2025", source_data$inputs, years = 1995:2024, ages = 2:11
), error = identity)
stopifnot(inherits(error, "error"),
          grepl("weight|maturity", conditionMessage(error)),
          identical(source_data, before),
          source_data$assessment$inputs_status == "partial")
# Isolate maturity readiness with test-only complete weights, never source data.
weight_fixture <- source_data$inputs[rep(1L, 30L * 10L), ]
weight_fixture$type <- "weight"
weight_fixture$measure <- "weight_at_age"
weight_fixture$basis <- "kg_per_fish"
weight_fixture$unit <- "kg"
weight_fixture$year <- rep(1995:2024, each = 10L)
weight_fixture$age <- rep(2:11, times = 30L)
weight_fixture$value <- 1
weight_fixture$source_reference <- "Test fixture, not assessment biology"
error <- tryCatch(database_to_tam_obs(
  "dfo_cod_3pn4rs_2025", rbind(source_data$inputs, weight_fixture),
  years = 1995:2024, ages = 2:11
), error = identity)
stopifnot(inherits(error, "error"),
          grepl("maturity", conditionMessage(error)),
          identical(source_data, before))
cat("3Pn4RS source-input tests passed.\n")
