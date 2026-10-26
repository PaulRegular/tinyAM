source("analysis/comp_assessments/tests/helper_stock.R")

root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))

assessment_id <- "ices_haddock_iceland_2025"
source_data <- read_assessment(assessment_id)
outputs <- source_data$outputs
n_at_age <- outputs[
  outputs$measure == "numbers_at_age" &
    outputs$year %in% 1979:2025 &
    !is.na(outputs$age),
  , drop = FALSE
]
total_numbers <- outputs[
  outputs$measure == "total_numbers" & outputs$year %in% 1979:2025,
  , drop = FALSE
]
expected <- aggregate(value ~ year, n_at_age, sum)
observed <- total_numbers[match(expected$year, total_numbers$year), ]

stopifnot(
  nrow(n_at_age) == 47L * 12L,
  nrow(total_numbers) == 47L,
  identical(sort(unique(total_numbers$year)), 1979:2025),
  all(total_numbers$type == "population"),
  all(total_numbers$unit == "thousand fish"),
  all(total_numbers$source_type == "derived_source_output"),
  all(is.na(total_numbers$age)),
  all(is.na(total_numbers$se)),
  all(is.na(total_numbers$lwr)),
  all(is.na(total_numbers$upr)),
  isTRUE(all.equal(observed$value, expected$value, tolerance = 1e-10)),
  all(grepl("age:1", total_numbers$notes, fixed = TRUE)),
  all(grepl("group:12+", total_numbers$notes, fixed = TRUE)),
  all(grepl("No standard error or interval was combined", total_numbers$notes, fixed = TRUE))
)

cat("Icelandic haddock 2025 derived total-abundance outputs passed.\n")
