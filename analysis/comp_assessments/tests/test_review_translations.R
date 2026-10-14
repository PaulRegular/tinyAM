pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

db <- read_database()
current <- db$assessments$assessment_id[db$assessments$is_current %in% TRUE]
record <- read.csv("analysis/comp_assessments/results/translation_review.csv",
  na.strings = c("", "NA")
)
stopifnot(
  setequal(record$assessment_id, current),
  !anyDuplicated(record[c("assessment_id", "candidate")]),
  !any(record$decision == "pending_review"),
  length(unique(record$baseline_revision)) == 1L,
  !anyNA(record$baseline_revision),
  sum(record$decision == "blocked") == 2L,
  all(is.na(record$recruitment_percent[record$recruitment_status %in% "non_equivalent"]))
)

fixture <- data.frame(
  assessment_id = "fixture",
  component = c("recruitment", "recruitment", "recruitment", "index", "index", "F", "M", "M"),
  setting = c(
    "stock_recruitment", "stock_recruitment", "deviation_distribution",
    "observation_correlation", "observation_correlation", "process", "process", "process_sd"
  ),
  value = c(
    "Beverton-Holt", "hockey-stick", "AR1", "Independent", "AR1 across ages",
    "Random walk with age-correlated increments", "AR1", "Fixed 0.075"
  ),
  source_reference = "offline fixture", notes = ""
)
audit <- audit_assumptions("fixture", fixture)
stopifnot(
  identical(
    audit$tinyam_support,
    c(
      "partially_supported", "unsupported", "partially_supported", "supported",
      "unsupported", "partially_supported", "partially_supported", "partially_supported"
    )
  ),
  grepl("supplied SD", audit$audit_notes[8]),
  !any(grepl("basic random walk|does not fit.*autocorrelated", audit$audit_notes))
)
cat("Historical review record and current-capability audit tests passed.\n")
