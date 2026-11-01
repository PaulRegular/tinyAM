root <- file.path('analysis', 'comp_assessments')
source(file.path(root, 'scripts', 'database',
  '498_screen_noaa_china_rockfish_2015.R'))
register <- read.csv(file.path(root, 'discovery', 'candidates.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
ids <- c('NOAA:14410', 'NOAA:14430', 'NOAA:14429')
units_expected <- c('noaa_china_rockfish_northern_pacific_coast',
  'noaa_china_rockfish_central_pacific_coast',
  'noaa_china_rockfish_southern_pacific_coast')
candidate <- register[match(ids, register$candidate_id), , drop = FALSE]
units <- read.csv(file.path(root, 'discovery', 'inventory_units.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
progress <- read.csv(file.path(root, 'discovery', 'inventory_progress.csv'),
  stringsAsFactors = FALSE, check.names = FALSE)
noaa_progress <- progress[progress$organization == 'NOAA', , drop = FALSE]
noaa <- aggregate(noaa_progress[c('inventory_rows', 'new', 'identity_unreviewed',
  'eligible', 'eligible_new', 'excluded', 'eligibility_unresolved',
  'screening_complete', 'screening_not_started')],
  list(organization = noaa_progress$organization), FUN = sum)
nw <- noaa_progress[noaa_progress$source_region == 'NWFSC / SWFSC', , drop = FALSE]
unit <- units[units$candidate_unit_id %in% units_expected, , drop = FALSE]
review <- paste(readLines(file.path(root, 'discovery', 'source_reviews',
  'noaa_china_rockfish_us_pacific_2015.md'), warn = FALSE), collapse = ' ')
stocks <- read.csv(file.path(root, 'database', 'stocks.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
assessments <- read.csv(file.path(root, 'database', 'assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
inputs <- read.csv(file.path(root, 'database', 'inputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
outputs <- read.csv(file.path(root, 'database', 'outputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)

stopifnot(nrow(candidate) == 3L, !anyDuplicated(candidate$candidate_unit_id),
  setequal(candidate$candidate_unit_id, units_expected),
  all(candidate$match_status == 'new'), all(candidate$screening_status == 'screened'),
  all(candidate$eligibility == 'eligible'), all(candidate$import_status == 'not_started'),
  all(grepl('2015 assessment', candidate$notes, fixed = TRUE)),
  all(grepl('not a tinyAM-readiness decision', candidate$decision_reason, fixed = TRUE)),
  grepl('three separate assessment units', review, fixed = TRUE),
  grepl('2019 action was a catch-only projection', review, fixed = TRUE),
  grepl('source_cache/discovery_noaa_china_rockfish_2019', review, fixed = TRUE),
  nrow(unit) == 3L, all(unit$inventory_entries == 1L),
  nrow(noaa) == 1L, noaa$inventory_rows == 432L, noaa$new == 282L,
  noaa$identity_unreviewed == 113L, noaa$eligible == 80L,
  noaa$eligible_new == 66L, noaa$excluded == 190L,
  noaa$eligibility_unresolved == 162L, noaa$screening_complete == 321L,
  noaa$screening_not_started == 111L,
  nrow(nw) == 1L, nw$new == 92L, nw$identity_unreviewed == 38L,
  nw$eligible == 32L, nw$eligible_new == 31L, nw$screening_complete == 100L,
  nw$screening_not_started == 38L,
  sum(progress$inventory_rows) == 905L, sum(progress$new) == 557L,
  sum(progress$identity_unreviewed) == 211L, sum(progress$eligible) == 161L,
  sum(progress$eligible_new) == 91L, sum(progress$excluded) == 522L,
  sum(progress$eligibility_unresolved) == 222L,
  sum(progress$screening_complete) == 794L,
  sum(progress$screening_not_started) == 111L,
  nrow(units) == 568L, sum(units$inventory_entries) == 604L,
  sum(!is.na(units$candidate_unit_id) & nzchar(units$candidate_unit_id)) == 468L,
  !any(stocks$scientific_name == 'Sebastes nebulosus', na.rm = TRUE),
  !any(assessments$assessment_id %in% units_expected, na.rm = TRUE),
  !any(inputs$assessment_id %in% units_expected, na.rm = TRUE),
  !any(outputs$assessment_id %in% units_expected, na.rm = TRUE))

cat('NOAA China rockfish 2015 area screens, source provenance, distinct-unit identities and no-import checks passed.\n')
