source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")
db <- read_database()
ids <- c("dfo_cod_eastern_georges_bank_2025_framework",
  "dfo_cod_eastern_georges_bank_2025_summary", "dfo_cod_eastern_georges_bank_2026_summary")
records <- lapply(ids, read_assessment, database = db)
framework <- records[[1L]]
inputs <- framework$inputs
totals <- subset(inputs, type == "catch")
weights <- subset(inputs, type == "catch_weight")
q <- framework$outputs
value <- function(x, year, fleet, age = NULL) {
  rows <- x$year == year & x$fleet == fleet
  if (!is.null(age)) rows <- rows & x$age == age
  as.numeric(x$value[rows])
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
  nrow(inputs) == 858L, nrow(totals) == 138L, nrow(weights) == 720L,
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
  !any(inputs$type %in% c("weight","maturity","M","index")),
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
cat("EGB cod source identity, run separation, input selection, native units and unresolved translation gaps passed.\n")
