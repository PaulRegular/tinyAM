root <- file.path('analysis', 'comp_assessments')
read_table <- function(name) read.csv(file.path(root, name), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'), check.names = FALSE)
canonical_paths <- c(stocks = 'database/stocks.csv',
  assessments = 'database/assessments.csv', assumptions = 'database/assumptions.csv',
  inputs = 'database/inputs.csv', outputs = 'database/outputs.csv')
canonical_before <- tools::md5sum(file.path(root, canonical_paths))

source(file.path(root, 'scripts', 'database',
  '533_screen_noaa_gulf_smoothhound_sedar39_2015.R'))

register <- read_table('discovery/candidates.csv')
row <- register[register$candidate_id == 'NOAA:14230', , drop = FALSE]
cache <- file.path(root, 'source_cache', 'discovery_noaa_smoothhound_2026')
file <- 'sedar39_gulf_smoothhound_2015.pdf'
manifest <- read.csv(file.path(cache, 'sources.csv'), stringsAsFactors = FALSE)
source_row <- manifest[manifest$file == file, , drop = FALSE]
expected_size <- 13513678
expected_hash <- 'E7E24243E13BE0179BF1376B5D2D7C23A962D513B1151EC0CC3E66F7C41EA8C9'

stopifnot(nrow(row) == 1L, row$official_stock_id == '14230',
  row$scientific_name == 'Mustelus spp.', row$assessment_area == 'Gulf of Mexico',
  row$candidate_unit_id == 'noaa_gulf_smoothhound_complex',
  is.na(row$matched_stock_id), row$match_status == 'new',
  row$screening_status == 'screened', row$eligibility == 'ineligible',
  row$import_status == 'not_started', row$review_date == '2026-11-02',
  grepl('no age-composition information', row$catch_age_evidence, fixed = TRUE),
  grepl('state-space surplus-production', row$decision_reason, fixed = TRUE),
  grepl('fully accepted', row$notes, fixed = TRUE),
  identical(canonical_before, tools::md5sum(file.path(root, canonical_paths))),
  nrow(source_row) == 1L, source_row$bytes == expected_size,
  source_row$sha256 == expected_hash, source_row$retrieved_on == '2026-11-02',
  nzchar(source_row$source_url), nzchar(source_row$pinned_revision),
  file.exists(file.path(cache, file)),
  file.info(file.path(cache, file))$size == expected_size,
  file.exists(file.path(root, 'discovery', 'source_reviews',
    'noaa_gulf_smoothhound_sedar39_2015.md')))

signature <- readBin(file.path(cache, file), 'raw', n = 5L)
stopifnot(identical(rawToChar(signature), '%PDF-'))

progress <- read_table('discovery/inventory_progress.csv')
noaa <- progress[progress$organization == 'NOAA', , drop = FALSE]
se <- noaa[noaa$source_region == 'SEFSC', , drop = FALSE]
units <- read_table('discovery/inventory_units.csv')
stopifnot(sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 644L, sum(progress$identity_unreviewed) == 120L,
  sum(progress$eligible) == 165L, sum(progress$excluded) == 573L,
  sum(progress$eligible_new) == 95L,
  sum(progress$eligibility_unresolved) == 167L,
  sum(progress$screening_complete) == 885L,
  sum(progress$screening_not_started) == 20L,
  sum(noaa$new) == 369L, sum(noaa$identity_unreviewed) == 22L,
  sum(noaa$excluded) == 241L, sum(noaa$eligibility_unresolved) == 107L,
  sum(noaa$screening_complete) == 412L, sum(noaa$screening_not_started) == 20L,
  nrow(se) == 1L, se$inventory_rows == 90L, se$new == 78L,
  se$identity_unreviewed == 9L, se$excluded == 45L,
  se$eligibility_unresolved == 27L, se$screening_complete == 81L,
  se$screening_not_started == 9L, nrow(units) == 654L,
  sum(units$inventory_entries) == 693L,
  sum(!is.na(units$candidate_unit_id) & nzchar(units$candidate_unit_id)) == 554L,
  sum(progress$identity_unlinked_rows) == 212L)

cat('NOAA Gulf smoothhound source-screening checks passed.\n')
