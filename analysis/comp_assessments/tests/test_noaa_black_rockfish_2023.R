root <- file.path('analysis', 'comp_assessments')
register <- read.csv(file.path(root, 'discovery', 'candidates.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
ids <- c('NOAA:166727', 'NOAA:10402', 'NOAA:11302')
rows <- register[match(ids, register$candidate_id), , drop = FALSE]
unit_id <- 'noaa_black_rockfish_pacific_coast_2023'
unit <- read.csv(file.path(root, 'discovery', 'inventory_units.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
unit <- unit[unit$unit_key == paste0('candidate:', unit_id), , drop = FALSE]
review <- paste(readLines(file.path(root, 'discovery', 'source_reviews',
  'noaa_black_rockfish_pacific_coast_2023.md'), warn = FALSE), collapse = '\n')
stocks <- read.csv(file.path(root, 'database', 'stocks.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
assessments <- read.csv(file.path(root, 'database', 'assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
inputs <- read.csv(file.path(root, 'database', 'inputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
outputs <- read.csv(file.path(root, 'database', 'outputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))

stopifnot(
  nrow(rows) == 3L,
  rows$candidate_unit_id[1L] == unit_id,
  rows$official_stock_id[1L] == '166727',
  rows$match_status[1L] == 'new',
  rows$screening_status[1L] == 'screened',
  rows$eligibility[1L] == 'eligible',
  rows$import_status[1L] == 'not_started',
  rows$scientific_name[1L] == 'Sebastes melanops',
  grepl('Washington, Oregon, northern California, central California',
    rows$assessment_area[1L], fixed = TRUE),
  grepl('numeric N-at-age and biomass-at-age tables',
    rows$output_evidence[1L], fixed = TRUE),
  all(rows$match_status[-1L] == 'ambiguous'),
  all(rows$screening_status[-1L] == 'screened'),
  all(rows$eligibility[-1L] == 'ineligible'),
  all(is.na(rows$candidate_unit_id[-1L]) | !nzchar(rows$candidate_unit_id[-1L])),
  nrow(unit) == 1L,
  unit$inventory_entries == 1L,
  unit$organizations == 'NOAA',
  grepl('assert that they are one biological population', review, fixed = TRUE),
  grepl('not all supplied as direct observation', review, fixed = TRUE),
  grepl('Do not pool regional', review, fixed = TRUE),
  !any(stocks$authority_stock_id == '166727', na.rm = TRUE),
  !any(stocks$stock_id == unit_id, na.rm = TRUE),
  !any(assessments$assessment_id == unit_id, na.rm = TRUE),
  !any(inputs$assessment_id == unit_id, na.rm = TRUE),
  !any(outputs$assessment_id == unit_id, na.rm = TRUE)
)

cat('NOAA Pacific Coast black rockfish screen checks passed.\n')