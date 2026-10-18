pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

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
cat("Current-capability audit classifications passed.\n")
