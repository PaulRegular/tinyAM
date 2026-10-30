# Summarize source-list rows without treating unscreened rows as new stocks.
root <- file.path('analysis', 'comp_assessments')
x <- read.csv(file.path(root, 'discovery', 'candidates.csv'), stringsAsFactors = FALSE, na.strings = c('', 'NA'))
e <- new.env(parent = emptyenv())
load(file.path(root, 'source_cache', 'discovery_noaa_stocksmart', 'stock_assessment_summary.rda'), envir = e)
noaa <- e$stock_assessment_summary
noaa <- noaa[order(noaa$stock_id, -as.numeric(noaa$assessment_year), noaa$assessment_id, na.last = TRUE), ]
noaa <- noaa[!duplicated(noaa$stock_id), ]
region <- rep(NA_character_, nrow(x))
is_noaa <- x$organization == 'NOAA'
region[is_noaa] <- noaa$science_center[match(x$official_stock_id[is_noaa], as.character(noaa$stock_id))]
is_dfo <- x$organization == 'DFO'
region[is_dfo] <- sub('^DFO:([^:]+):.*$', '\\1', x$candidate_id[is_dfo])
is_ices <- x$organization == 'ICES'
region[is_ices] <- x$assessment_area[is_ices]
region[!is.na(region)] <- trimws(gsub('[[:space:]]+', ' ', region[!is.na(region)]))
region[is.na(region) | !nzchar(trimws(region))] <- 'Unspecified'
key <- data.frame(organization = x$organization, source_region = region, stringsAsFactors = FALSE)
metrics <- data.frame(inventory_rows = 1L, matched = as.integer(x$match_status == 'matched'),
  new = as.integer(x$match_status == 'new'), ambiguous = as.integer(x$match_status == 'ambiguous'),
  identity_unreviewed = as.integer(x$match_status == 'unreviewed'),
  eligible = as.integer(x$eligibility == 'eligible'), excluded = as.integer(x$eligibility == 'ineligible'),
  eligible_new = as.integer(x$match_status == 'new' & x$eligibility == 'eligible'),
  eligibility_unresolved = as.integer(x$eligibility == 'unresolved'),
  screening_not_started = as.integer(x$screening_status == 'not_started'),
  screening_in_progress = as.integer(x$screening_status == 'in_progress'),
  screening_complete = as.integer(x$screening_status == 'screened'))
out <- aggregate(metrics, by = key, FUN = sum)
out <- out[order(out$organization, out$source_region), ]
out$matched_units <- 0L
out$crosswalked_candidate_units <- 0L
out$identity_unlinked_rows <- 0L
for (i in seq_len(nrow(out))) {
  rows <- x$organization == out$organization[i] & region == out$source_region[i]
  matched_ids <- unique(x$matched_stock_id[rows & !is.na(x$matched_stock_id) & nzchar(x$matched_stock_id)])
  candidate_ids <- unique(x$candidate_unit_id[rows & !is.na(x$candidate_unit_id) & nzchar(x$candidate_unit_id)])
  linked <- (!is.na(x$matched_stock_id) & nzchar(x$matched_stock_id)) |
    (!is.na(x$candidate_unit_id) & nzchar(x$candidate_unit_id))
  out$matched_units[i] <- length(matched_ids)
  out$crosswalked_candidate_units[i] <- length(candidate_ids)
  out$identity_unlinked_rows[i] <- sum(rows & !linked)
}
stopifnot(sum(out$inventory_rows) == nrow(x),
  all(out$inventory_rows == out$matched + out$new + out$ambiguous + out$identity_unreviewed),
  all(out$inventory_rows == out$eligible + out$excluded + out$eligibility_unresolved),
  all(out$inventory_rows == out$screening_not_started + out$screening_in_progress + out$screening_complete))
write.csv(out, file.path(root, 'discovery', 'inventory_progress.csv'), row.names = FALSE, na = '')
print(aggregate(out[3:ncol(out)], list(organization = out$organization), sum), row.names = FALSE)
