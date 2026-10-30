register <- read.csv('analysis/comp_assessments/discovery/candidates.csv',
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
ids <- c('DFO:Maritimes:Eel (Large)', 'DFO:Maritimes:Elvers')
rows <- register[match(ids, register$candidate_id), , drop = FALSE]
stocks <- read.csv('analysis/comp_assessments/database/stocks.csv', stringsAsFactors = FALSE)
assessments <- read.csv('analysis/comp_assessments/database/assessments.csv', stringsAsFactors = FALSE)
inputs <- read.csv('analysis/comp_assessments/database/inputs.csv', stringsAsFactors = FALSE)
outputs <- read.csv('analysis/comp_assessments/database/outputs.csv', stringsAsFactors = FALSE)
unit_summary <- read.csv('analysis/comp_assessments/discovery/inventory_units.csv',
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
unit <- unit_summary[!is.na(unit_summary$candidate_unit_id) &
  unit_summary$candidate_unit_id == 'dfo_american_eel_maritimes', , drop = FALSE]
cache <- 'analysis/comp_assessments/source_cache/discovery_dfo_american_eel_maritimes_2026'
manifest <- read.csv(file.path(cache, 'sources.csv'), stringsAsFactors = FALSE)
files <- file.path(cache, manifest$file)

stopifnot(nrow(rows) == 2L,
  all(rows$scientific_name == 'Anguilla rostrata'),
  all(rows$candidate_unit_id == 'dfo_american_eel_maritimes'),
  all(rows$match_status == 'new'),
  all(rows$screening_status == 'screened'),
  all(rows$eligibility == 'ineligible'),
  all(rows$import_status == 'not_started'),
  all(vapply(rows[c('catch_age_evidence', 'survey_age_evidence', 'weights_evidence',
    'maturity_evidence', 'mortality_evidence', 'timing_evidence', 'output_evidence')],
    function(x) all(!is.na(x) & nzchar(x)), logical(1L))),
  nrow(unit) == 1L, unit$inventory_entries == 2L,
  unit$organizations == 'DFO', unit$screening_status == 'screened',
  unit$eligibility == 'ineligible',
  !any(stocks$scientific_name == 'Anguilla rostrata'),
  !any(assessments$assessment_id == 'dfo_american_eel_maritimes'),
  !any(inputs$assessment_id == 'dfo_american_eel_maritimes'),
  !any(outputs$assessment_id == 'dfo_american_eel_maritimes'),
  nrow(manifest) == 3L, all(file.exists(files)),
  all(file.info(files)$size == manifest$bytes),
  all(grepl('^[[:xdigit:]]{64}$', manifest$sha256)),
  any(grepl('41234996.pdf', manifest$url, fixed = TRUE)),
  file.exists('analysis/comp_assessments/discovery/source_reviews/dfo_american_eel_maritimes_2026.md'))
cat('DFO Maritimes American eel identity, source review and no-import checks passed.\n')
