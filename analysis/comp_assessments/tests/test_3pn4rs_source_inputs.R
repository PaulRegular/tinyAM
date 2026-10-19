root <- file.path("analysis", "comp_assessments")
inputs <- read.csv(file.path(root, "database", "inputs.csv"),
  stringsAsFactors = FALSE
)
inputs <- inputs[inputs$assessment_id == "dfo_cod_3pn4rs_2025", ]
index <- inputs[inputs$type == "index", ]
weights <- inputs[inputs$type == "catch_weight", ]
stock_weights <- inputs[inputs$type == "weight", ]
maturity <- inputs[inputs$type == "maturity", ]
M_inputs <- inputs[inputs$type == "M", ]

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
  nrow(stock_weights) == 40L * 10L,
  setequal(stock_weights$year, 1985:2024),
  setequal(stock_weights$age, 2:11),
  all(stock_weights$source_type == "digitized"),
  all(stock_weights$source_reference == "DFO Research Document 2026/010 Figure 90"),
  all(stock_weights$unit == "kg"),
  all(grepl("blue fitted stock-weight curve", stock_weights$notes)),
  nrow(maturity) == 42L * 16L,
  setequal(maturity$year, 1983:2024),
  setequal(maturity$age, 0:15),
  all(maturity$source_type == "digitized"),
  all(maturity$source_reference == "DFO Research Document 2026/010 Figure 70"),
  all(maturity$value >= 0 & maturity$value <= 1),
  all(vapply(split(maturity, maturity$year), function(x) {
    all(diff(x$value[order(x$age)]) >= 0)
  }, logical(1))),
  sum(M_inputs$source_type == "official_document") == 41L * 8L,
  all(M_inputs$value[M_inputs$age == 4 & M_inputs$year >= 1984] == 0.45),
  all(M_inputs$value[M_inputs$age >= 5 & M_inputs$year >= 1984] == 0.15),
  all(grepl("process means, not estimated M surfaces", M_inputs$notes[
    M_inputs$source_type == "official_document"
  ])),
  all(is.finite(inputs$value)),
  sum(inputs$type == "weight") == 400L,
  sum(inputs$type == "maturity") == 672L
)

source(file.path(root, "R", "run_assessment.R"))
source_data <- read_assessment("dfo_cod_3pn4rs_2025", read_database())
before <- source_data
obs <- database_to_tam_obs(
  "dfo_cod_3pn4rs_2025", source_data$inputs,
  years = 1995:2024, ages = 2:11,
  surveys = "Sentinel mobile",
  sampling_times = c("Sentinel mobile" = 0.54),
  assumptions = source_data$assumptions
)
stopifnot(
  nrow(obs$catch) == 30L * 10L,
  nrow(obs$index) == 30L * 10L,
  nrow(obs$weight) == 30L * 10L,
  nrow(obs$maturity) == 30L * 10L,
  all(is.finite(obs$pop_covariates$M_assumption)),
  identical(tinyAM::check_obs(obs), TRUE),
  identical(source_data, before),
  source_data$assessment$inputs_status == "partial"
)
cat("3Pn4RS source-input tests passed.\n")
