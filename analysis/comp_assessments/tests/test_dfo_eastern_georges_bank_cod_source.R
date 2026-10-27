source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")
db <- read_database()
ids <- c("dfo_cod_eastern_georges_bank_2025_framework",
  "dfo_cod_eastern_georges_bank_2025_summary", "dfo_cod_eastern_georges_bank_2026_summary")
records <- lapply(ids, read_assessment, database = db)
framework <- records[[1L]]
inputs <- framework$inputs
totals <- subset(inputs, type == "catch")
m <- subset(inputs, type == "M" & measure == "natural_mortality_at_age")
weights <- subset(inputs, type == "catch_weight")
survey_indices <- subset(inputs, type == "index" & measure == "total_biomass")
q <- framework$outputs
value <- function(x, year, fleet, age = NULL) {
  rows <- x$year == year & if (nzchar(fleet)) x$fleet == fleet else is.na(x$fleet)
  if (!is.null(age)) rows <- rows & x$age == age
  as.numeric(x$value[rows])
}
index_value <- function(x, year, survey) {
  as.numeric(x$value[x$year == year & x$survey == survey])
}
stopifnot(
  framework$stocks$authority_stock_id == "5Zejm (Eastern Georges Bank)",
  framework$stocks$charbonneau_id == "NEFSC-TRAC-CERT_morhua_GB",
  framework$stocks$ocean == "Northwest Atlantic",
  framework$assessment$is_current,
  framework$assessment$terminal_year == 2023,
  all(vapply(records, function(x) x$assessment$is_applied, logical(1))),
  identical(vapply(records, function(x) x$assessment$assessment_year, integer(1)), c(2025L,2025L,2026L)),
  identical(vapply(records, function(x) x$assessment$terminal_year, integer(1)), c(2023L,2024L,2025L)),
  all(vapply(records, function(x) all(unlist(x$assessment[c("inputs_status","outputs_status","assumptions_status")]) == "partial"), logical(1))),
  !records[[2L]]$assessment$is_current, !records[[3L]]$assessment$is_current,
  nrow(inputs) == 1410L, nrow(totals) == 138L, nrow(weights) == 720L,
  nrow(m) == 414L, setequal(m$year, 1978:2023), setequal(m$age, 1:9),
  all(table(m$age, m$year) == 1L), all(m$value > 0), all(m$unit == "per_year"),
  all(m$source_type == "digitized"), all(grepl("Figure 11", m$source_reference)),
  all(grepl("No smoother or interpolation", m$transformation)),
  all(grepl("9+", m$notes[m$age == 9], fixed = TRUE)),
  value(m, 1980, "", 1) == .61, value(m, 2001, "", 1) == 1.55,
  value(m, 2019, "", 9) == .273, value(m, 2023, "", 9) == .255,
  setequal(unique(totals$year), 1978:2023),
  all(table(totals$fleet) == 46L), length(unique(totals$fleet)) == 3L,
  all(totals$unit == "metric tonnes"), all(is.na(totals$age)),
  value(totals,1978,"Canada reported") == 8875,
  value(totals,1978,"US_Reported EGB") == 5502,
  value(totals,2007,"US_Misattributed EGB") == 233,
  value(totals,2023,"US_Reported EGB") == 31,
  all(as.numeric(subset(totals, fleet == "US_Misattributed EGB" &
    (year < 2005 | year > 2016))$value) == 1),
  setequal(unique(weights$year),1978:2022), setequal(unique(weights$age),1:16),
  all(weights$unit == "kg per fish"),
  value(weights,1978,"Canada reported",16) == 19.7,
  value(weights,1984,"Canada reported",1) == 0,
  value(weights,2022,"Canada reported",9) == 12.8,
  nrow(survey_indices) == 138L,
  setequal(survey_indices$survey, c("DFO winter RV", "NMFS spring RV", "NMFS fall RV")),
  all(table(survey_indices$survey) == 46L),
  all(vapply(split(survey_indices$year, survey_indices$survey), function(x) setequal(x, 1978:2023), logical(1))),
  all(survey_indices$unit == "metric tonnes"), all(survey_indices$basis == "biomass"),
  all(survey_indices$region == "Eastern Georges Bank"),
  all(survey_indices$source_type == "official_table"), all(is.na(survey_indices$age)),
  identical(as.integer(tapply(!is.na(survey_indices$value), survey_indices$survey, sum)), c(37L,44L,44L)),
  index_value(survey_indices, 1986, "DFO winter RV") == 18633,
  index_value(survey_indices, 1978, "NMFS spring RV") == 22731,
  index_value(survey_indices, 1978, "NMFS fall RV") == 12213,
  is.na(index_value(survey_indices, 2022, "DFO winter RV")),
  is.na(index_value(survey_indices, 2023, "NMFS spring RV")),
  is.na(index_value(survey_indices, 2023, "NMFS fall RV")),
  !any(inputs$measure == "numbers_at_age"),
  all(is.na(inputs$sampling_time)),
  nrow(q) == 3L, all(q$measure == "q"), all(is.na(q$year)), all(is.na(q$age)),
  identical(as.numeric(q$value), c(.434,.223,.109)),
  setequal(q$survey, c("DFO winter RV","NMFS spring RV","NMFS fall RV")),
  grepl("NTS", subset(framework$assumptions, setting == "Canadian_age_data_selection")$value),
  grepl("not AR1 recruitment", subset(framework$assumptions, setting == "process" & component == "recruitment")$notes)
)
summary_2025 <- records[[2L]]
summary_2026 <- records[[3L]]
ssb_2025 <- subset(summary_2025$outputs, measure == "SSB")
ssb_2026 <- subset(summary_2026$outputs, measure == "SSB")
all_outputs <- do.call(rbind, lapply(records, `[[`, "outputs"))
stopifnot(
  nrow(summary_2025$inputs) == 10L, nrow(summary_2026$inputs) == 20L,
  setequal(summary_2025$inputs$year, 2020:2024),
  setequal(summary_2026$inputs$year, 2016:2025),
  value(summary_2026$inputs,2025,"Canada") == 297,
  value(summary_2026$inputs,2025,"US") == 20,
  value(summary_2026$inputs,2023,"US") == 32,
  all(as.numeric(ssb_2025$value) == c(7.96,10.9)), all(ssb_2025$unit == "kt"),
  ssb_2026$value == 8469, ssb_2026$year == 2025,
  grepl("median", ssb_2026$notes),
  !any(ssb_2026$year == 2022),
  !any(all_outputs$measure %in% c("Fbar","numbers_at_age","fishing_mortality_at_age","natural_mortality_at_age")),
  sum(all_outputs$measure == "fishing_mortality_intensity") == 2L,
  all(is.na(all_outputs$se)), all(is.na(all_outputs$lwr)), all(is.na(all_outputs$upr)),
  all(nzchar(all_outputs$source_reference)),
  !any(grepl("sensitivity|Run 5|Run 6|Run 9", q$source_reference)),
  file.exists("analysis/comp_assessments/source_reviews/ck_nefsc_trac_cert_morhua_gb.md")
)
before <- framework
error <- tryCatch(database_to_tam_obs(ids[1], framework$inputs,
  years = 1978:2023, ages = 1:9, assumptions = framework$assumptions), error = identity)
stopifnot(inherits(error,"error"), identical(framework,before))
cat("EGB cod source identity, run separation, input selection, native units, digitized fixed M and unresolved translation gaps passed.\n")
