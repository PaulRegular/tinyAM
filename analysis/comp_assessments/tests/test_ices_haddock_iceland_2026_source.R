source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_haddock_iceland_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
original <- source_data
inputs <- source_data$inputs
outputs <- source_data$outputs
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", ]
catch_weight <- inputs[inputs$type == "catch_weight", ]
stock_weight <- inputs[inputs$type == "weight", ]
maturity <- inputs[inputs$type == "maturity", ]
smb <- inputs[!is.na(inputs$survey) & inputs$survey == "IS-SMB", ]
smh <- inputs[!is.na(inputs$survey) & inputs$survey == "IS-SMH", ]
n_at_age <- outputs[outputs$measure == "numbers_at_age", ]
f_at_age <- outputs[outputs$measure == "fishing_mortality_at_age", ]
summary <- outputs[outputs$measure %in% c("SSB", "reference_biomass", "recruitment"), ]

stopifnot(
  nrow(source_data$stock) == 1L,
  source_data$stock$authority_stock_id == "had.27.5a",
  source_data$stock$ocean == "Northeast Atlantic",
  nrow(source_data$assessment) == 1L,
  isTRUE(source_data$assessment$is_current),
  source_data$assessment$assessment_year == 2026L,
  source_data$assessment$terminal_year == 2026L,
  source_data$assessment$estimate_terminal_year == 2026L,
  source_data$assessment$inputs_status == "partial",
  source_data$assessment$outputs_status == "partial",
  nrow(catch) == 47L * 12L,
  setequal(unique(catch$year), 1979:2025),
  all(catch$unit == "thousand fish"),
  catch$value[catch$year == 1979 & catch$age == 1] == 785,
  catch$value[catch$year == 2025 & catch$age == 12] == 73,
  nrow(catch_weight) == 47L * 12L,
  max(catch_weight$year) == 2025L,
  !any(catch_weight$year == 2026L),
  nrow(stock_weight) == 48L * 12L,
  max(stock_weight$year) == 2026L,
  nrow(maturity) == 48L * 12L,
  nrow(smb) == 467L,
  min(smb$year) == 1985L,
  max(smb$year) == 2026L,
  nrow(smh) == 319L,
  min(smh$year) == 1995L,
  max(smh$year) == 2025L,
  all(smb$unit == "survey index (unit unresolved)"),
  all(smh$unit == "survey index (unit unresolved)"),
  all(smb$sampling_time == 0.20),
  all(smh$sampling_time == 0.80),
  all(inputs$value[inputs$measure == "natural_mortality_at_age"] == 0.2),
  nrow(n_at_age) == 48L * 12L,
  max(n_at_age$year) == 2026L,
  all(n_at_age$age_group[n_at_age$age == 12] == "12+"),
  nrow(f_at_age) == 47L * 12L,
  max(f_at_age$year) == 2025L,
  all(f_at_age$age_group[f_at_age$age == 12] == "12+"),
  nrow(summary) == 48L * 3L,
  as.numeric(summary$value[summary$measure == "SSB" & summary$year == 2026]) == 139691,
  as.numeric(summary$lwr[summary$measure == "SSB" & summary$year == 2026]) == 112156,
  as.numeric(summary$upr[summary$measure == "SSB" & summary$year == 2026]) == 173986,
  as.numeric(summary$value[summary$measure == "recruitment" & summary$year == 2026]) == 235664,
  !any(outputs$measure %in% "Fbar"),
  length(unique(source_data$assumptions$setting)) >= 10L
)

other_assessment <- database$assessments[
  database$assessments$assessment_id == "ices_haddock_iceland_2025", ]
stopifnot(nrow(other_assessment) == 1L, !isTRUE(other_assessment$is_current))
review <- paste(readLines(file.path(root, "source_reviews", "ices_haddock_iceland_2026.md"), warn = FALSE), collapse = " ")
stopifnot(
  grepl("4,000 g", review, fixed = TRUE),
  grepl("role in the assessment projections is not specified", review, fixed = TRUE),
  grepl("Fbar", review, fixed = TRUE),
  grepl("interval coverage", review, fixed = TRUE),
  grepl("not used to fill gaps", review, fixed = TRUE)
)

cat("Icelandic haddock 2026 source records passed.\n")
