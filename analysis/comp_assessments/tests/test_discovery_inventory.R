register <- read.csv('analysis/comp_assessments/discovery/candidates.csv', stringsAsFactors = FALSE, na.strings = c('', 'NA'))
stocks <- read.csv('analysis/comp_assessments/database/stocks.csv', stringsAsFactors = FALSE)
stopifnot(nrow(register) > 0L, !anyNA(register$candidate_id), !anyDuplicated(register$candidate_id),
  all(register$organization %in% c('NOAA', 'ICES', 'DFO')),
  all(register$match_status %in% c('unreviewed', 'matched', 'new', 'ambiguous')),
  all(register$screening_status %in% c('not_started', 'in_progress', 'screened')),
  all(register$eligibility %in% c('unresolved', 'eligible', 'ineligible')),
  all(register$import_status %in% c('not_started', 'in_progress', 'imported', 'blocked')),
  all(is.na(register$matched_stock_id) | register$matched_stock_id %in% stocks$stock_id))
ices_matches <- register[register$organization == 'ICES' & register$match_status == 'matched', ]
if (nrow(ices_matches)) {
  stock_index <- match(ices_matches$matched_stock_id, stocks$stock_id)
  stopifnot(all(tolower(ices_matches$official_stock_id) == tolower(stocks$authority_stock_id[stock_index])))
}
eligible <- register$eligibility == 'eligible'
for (column in c('catch_age_evidence', 'survey_age_evidence', 'weights_evidence',
                 'maturity_evidence', 'mortality_evidence', 'timing_evidence',
                 'output_evidence', 'detailed_source_url', 'review_date')) {
  stopifnot(!any(is.na(register[[column]][eligible]) | !nzchar(register[[column]][eligible])))
}
cat('Discovery register structure and eligibility evidence checks passed.\n')