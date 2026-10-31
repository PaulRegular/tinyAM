source('analysis/comp_assessments/scripts/database/461_screen_noaa_iccat_white_marlin_2025.R')
source('analysis/comp_assessments/scripts/database/195_summarize_discovery.R')
source('analysis/comp_assessments/scripts/database/196_summarize_inventory_units.R')
root <- file.path('analysis', 'comp_assessments')
candidate <- read.csv(file.path(root, 'discovery', 'candidates.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
row <- candidate[candidate$candidate_id == 'NOAA:10955', , drop = FALSE]
stocks <- read.csv(file.path(root, 'database', 'stocks.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
assessments <- read.csv(file.path(root, 'database', 'assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
inputs <- read.csv(file.path(root, 'database', 'inputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
outputs <- read.csv(file.path(root, 'database', 'outputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
manifest <- read.csv(file.path(root, 'source_cache',
  'discovery_noaa_iccat_white_marlin_2025', 'sources.csv'),
  stringsAsFactors = FALSE)
progress <- read.csv(file.path(root, 'discovery', 'inventory_progress.csv'),
  stringsAsFactors = FALSE)
units <- read.csv(file.path(root, 'discovery', 'inventory_units.csv'),
  stringsAsFactors = FALSE)

stopifnot(nrow(row) == 1L,
  row$candidate_unit_id == 'iccat_atlantic_white_marlin',
  row$match_status == 'new',
  row$screening_status == 'screened',
  row$eligibility == 'ineligible',
  row$import_status == 'not_started',
  grepl('JABBA', row$decision_reason, fixed = TRUE),
  grepl('age-composition', row$catch_age_evidence, fixed = TRUE),
  grepl('aggregate fleet CPUE', row$survey_age_evidence, fixed = TRUE),
  grepl('2025_WHM_SA_ENG.pdf', row$detailed_source_url, fixed = TRUE),
  all(file.exists(file.path(root, 'source_cache',
    'discovery_noaa_iccat_white_marlin_2025', manifest$file))),
  !any(stocks$authority_stock_id == '10955', na.rm = TRUE),
  !any(assessments$assessment_id == 'iccat_atlantic_white_marlin', na.rm = TRUE),
  !any(inputs$assessment_id == 'iccat_atlantic_white_marlin', na.rm = TRUE),
  !any(outputs$assessment_id == 'iccat_atlantic_white_marlin', na.rm = TRUE),
  sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 434L,
  sum(progress$identity_unreviewed) == 336L,
  sum(progress$excluded) == 420L,
  sum(progress$eligibility_unresolved) == 344L,
  sum(progress$screening_complete) == 669L,
  sum(progress$screening_not_started) == 236L,
  nrow(units) == 458L,
  sum(!is.na(units$candidate_unit_id) & nzchar(units$candidate_unit_id)) == 358L)
cat('NOAA:10955 ICCAT white marlin source-screening checks passed.\n')
