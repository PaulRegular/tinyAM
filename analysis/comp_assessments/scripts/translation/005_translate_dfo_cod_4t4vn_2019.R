root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "database_to_tiny_obs.R"))
source(file.path(root, "R", "database_to_tiny_M.R"))
source(file.path(root, "R", "audit_assumptions.R"))

read_table <- function(name) {
  read.csv(file.path(root, "database", name), stringsAsFactors = FALSE,
           na.strings = c("", "NA"), check.names = FALSE)
}

assessment_id <- "dfo_cod_4t4vn_2019"
years <- 1986:2018
ages <- 2:11
sampling_times <- c("DFO September RV survey" = 0.75,
                    "Mobile Sentinel August survey" = 0.625,
                    "Longline Sentinel survey" = 0.67)

inputs <- read_table("inputs.csv")
assumptions <- read_table("assumptions.csv")

obs <- database_to_tiny_obs(
  assessment_id, inputs, years = years, ages = ages,
  weight_survey = "DFO September RV survey",
  sampling_times = sampling_times,
  surveys = names(sampling_times),
  exclude_index_years = list("DFO September RV survey" = 2003)
)

m <- database_to_tiny_M(assessment_id, inputs, assumptions,
                        years = years, ages = ages)
m_start <- c("2" = 0.65, "3" = 0.65, "4" = 0.65,
             "5" = 0.15, "6" = 0.15, "7" = 0.15, "8" = 0.15,
             "9" = 0.15, "10" = 0.15, "11" = 0.15)
obs$weight$M_assumption <- unname(m_start[as.character(obs$weight$age)])

fit_settings <- list(
  years = years,
  ages = ages,
  N_settings = list(process = "off", init = "exp"),
  F_settings = list(process = "rw", mu_form = NULL),
  M_settings = list(process = "rw", mu_form = NULL,
                    mu_supplied = ~ M_assumption,
                    age_breaks = c(2, 5, 9, 11),
                    first_dev_year = min(years)),
  catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
  index_settings = list(q_form = ~ 0 + q_key,
                        sd_form = ~ 0 + survey,
                        fill_missing = FALSE)
)

audit <- audit_assumptions(assessment_id, assumptions)
longline_timing <- audit$component == "index" &
  audit$survey == "Longline Sentinel survey" &
  audit$setting == "sampling_time"
audit$tinyam_support[longline_timing] <- "partially_supported"
audit$audit_notes[longline_timing] <- paste(
  "The report gives a July-October sampling window. The translation uses its",
  "approximate midpoint (0.67); exact within-year timing is unavailable."
)

translation_decisions <- data.frame(
  component = c("model years", "model ages", "catch", "survey indices",
                "survey timing", "natural mortality", "recruitment", "fishing mortality"),
  choice = c(
    "1986-2018",
    "Ages 2-11, with age 11 treated as the tinyAM plus age",
    "Use available landed numbers-at-age, convert thousand fish to fish, and leave age 2 missing",
    "Reconstruct numbers-at-age from aggregate biomass and age proportions; use survey weights, with RV weights for longline",
    "Use 0.75 for September RV, 0.625 for August mobile, and 0.67 for the July-October longline midpoint",
    "Use an age-blocked RW for 2-4, 5-8, and 9-11; source prior means are starting values only",
    "Disable extra cohort residuals; tinyAM recruitment remains a random walk rather than the source AR recruitment-rate process",
    "Use an age-year RW as the closest available approximation to period-specific fleet selectivity"
  ),
  reason = c(
    "The available September RV weights have gaps in 1980 and 1985; no values are interpolated.",
    "The recorded weight-at-age surface ends at age 11, while the source assessment uses 12+.",
    "The source assessment also fits total catch and age composition separately; the available landings table is incomplete for that likelihood.",
    "The source model fits aggregate index and composition likelihoods separately; tinyAM uses age-specific lognormal observations. RV 2003 is excluded because the report says it was not used in population-model fitting.",
    "Month-level timing is only approximate; longline dates are available only as a July-October window.",
    "The source fixes RW innovation SD at 0.075 and places priors on starting M levels; tinyAM estimates the RW SD and has no matching starting-level prior.",
    "The source has autocorrelated recruitment-rate deviations; tinyAM uses a basic random walk for recruitment.",
    "The source uses logistic selectivity curves for four time periods; tinyAM does not reproduce that structure directly."
  ),
  stringsAsFactors = FALSE
)

if (requireNamespace("tinyAM", quietly = TRUE)) tinyAM::check_obs(obs)

grid_key <- function(x, columns) {
  do.call(paste, c(lapply(x[columns], as.character), sep = "\r"))
}
expected_grid <- expand.grid(year = years, age = ages)
for (name in c("catch", "weight", "maturity")) {
  if (anyDuplicated(grid_key(obs[[name]], c("year", "age"))) ||
      !setequal(grid_key(obs[[name]], c("year", "age")),
                grid_key(expected_grid, c("year", "age")))) {
    .translation_abort(paste(name, "does not cover the requested year-age grid."))
  }
}
if (anyNA(obs$weight$obs) || anyNA(obs$maturity$obs) ||
    any(!is.finite(obs$index$obs)) || anyNA(obs$index$samp_time) ||
    any(obs$index$samp_time < 0 | obs$index$samp_time > 1) ||
    anyNA(obs$index$q_key)) {
  .translation_abort("Translated observation tables contain invalid values or index fields.")
}

out_dir <- file.path(root, "results", "processed_inputs", assessment_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(obs$catch, file.path(out_dir, "catch.csv"), row.names = FALSE, na = "")
write.csv(obs$index, file.path(out_dir, "index.csv"), row.names = FALSE, na = "")
write.csv(obs$weight, file.path(out_dir, "weight.csv"), row.names = FALSE, na = "")
write.csv(obs$maturity, file.path(out_dir, "maturity.csv"), row.names = FALSE, na = "")
settings_lines <- capture.output(dput(fit_settings))
writeLines(sub("[[:space:]]+$", "", settings_lines),
           file.path(out_dir, "fit_settings.R"))
write.csv(audit, file.path(out_dir, "assumption_audit.csv"),
          row.names = FALSE, na = "")
write.csv(attr(obs, "translation")$source_provenance,
          file.path(out_dir, "source_provenance.csv"),
          row.names = FALSE, na = "")
write.csv(translation_decisions, file.path(out_dir, "translation_decisions.csv"),
          row.names = FALSE, na = "")
write.csv(data.frame(status = m$status, notes = m$notes),
          file.path(out_dir, "M_source_status.csv"), row.names = FALSE)

cat("Wrote Southern Gulf cod translation to ", out_dir, "\n", sep = "")
