root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_committed_assessment.R"))
source(file.path(root, "R", "database_to_tiny_obs.R"))
source(file.path(root, "R", "database_to_tiny_M.R"))
source(file.path(root, "R", "audit_assumptions.R"))

assessment_id <- "afsc_pollock_ebs_2024"
source_data <- read_committed_assessment(assessment_id)
years <- 1964:2024
ages <- 1:15
surveys <- c("NMFS bottom-trawl VAST", "NMFS acoustic-trawl",
             "NMFS acoustic-trawl age-1 index")

obs <- database_to_tiny_obs(
  assessment_id, source_data$inputs, years = years, ages = ages,
  surveys = surveys, maturity_reference_year = 1964,
  maturity_multiplier = 0.5
)
M <- database_to_tiny_M(assessment_id, source_data$inputs,
                         source_data$assumptions, years = years, ages = ages)
audit <- audit_assumptions(assessment_id, source_data$assumptions)

translation_decisions <- data.frame(
  component = c("source revision", "model years", "model ages", "catch",
                "BTS", "ATS ages 2-15", "ATS age 1", "CPUE/AVO",
                "maturity", "M", "model settings"),
  decision = c(
    source_data$commit,
    "1964-2024",
    "Ages 1-15; the native model plus group is age 15",
    "Translate total biomass plus number compositions with source catch weights; 2024 composition is missing",
    "Translate biomass plus number composition using stock weights as an explicit approximation",
    "Assign full source ATS biomass to ages 2-15 using stock weights, as an approximation; age 1 is fitted separately",
    "Retain the direct age-1 index; the source excludes 2024 by its uncertainty rule",
    "Omit the total-only indices because tinyAM requires age-specific indices",
    "Repeat the documented 1964 maturity vector and multiply by 0.5 for the source female-SSB convention",
    paste("No numerical M surface is recorded in inputs.csv.", M$notes,
          "The documented baseline vector's role alongside predation mortality is unresolved."),
    "Not generated until the source M treatment and remaining process/likelihood mappings are resolved"
  ),
  stringsAsFactors = FALSE
)

out_dir <- file.path(root, "results", "processed_inputs", assessment_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (name in c("catch", "index", "weight", "maturity")) {
  write.csv(obs[[name]], file.path(out_dir, paste0(name, ".csv")),
            row.names = FALSE, na = "")
}
write.csv(audit, file.path(out_dir, "assumption_audit.csv"),
          row.names = FALSE, na = "")
write.csv(translation_decisions,
          file.path(out_dir, "translation_decisions.csv"),
          row.names = FALSE, na = "")
write.csv(data.frame(status = M$status, notes = M$notes),
          file.path(out_dir, "M_source_status.csv"), row.names = FALSE)
writeLines(source_data$commit, file.path(out_dir, "database_revision.txt"))
cat("Wrote EBS pollock translation to", out_dir, "from", source_data$commit, "\n")
