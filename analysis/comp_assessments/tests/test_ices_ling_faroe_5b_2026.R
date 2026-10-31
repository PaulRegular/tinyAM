root <- file.path('analysis', 'comp_assessments')
old_dry_run <- Sys.getenv('TAM_SOURCE_DRY_RUN', unset = NA_character_)
Sys.setenv(TAM_SOURCE_DRY_RUN = '1')
screen <- new.env(parent = globalenv())
source(file.path(root, 'scripts/database', '454_screen_ices_ling_5b_2026.R'), local = screen)
if (is.na(old_dry_run)) Sys.unsetenv('TAM_SOURCE_DRY_RUN') else
  Sys.setenv(TAM_SOURCE_DRY_RUN, old_dry_run)

register <- read.csv(file.path(root, 'discovery/candidates.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'), check.names = FALSE)
candidate <- register[register$candidate_id == 'ICES:lin.27.5b', , drop = FALSE]
review <- screen$review
cache <- file.path(root, 'source_cache/discovery_ices_ling_5b_2026')
manifest <- read.csv(file.path(cache, 'sources.csv'), stringsAsFactors = FALSE)
stocks <- read.csv(file.path(root, 'database/stocks.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
assessments <- read.csv(file.path(root, 'database/assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
inputs <- read.csv(file.path(root, 'database/inputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
outputs <- read.csv(file.path(root, 'database/outputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
units <- read.csv(file.path(root, 'discovery/inventory_units.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
unit <- units[!is.na(units$candidate_unit_id) &
  units$candidate_unit_id == 'ices_ling_faroe_5b', , drop = FALSE]
source_review <- paste(readLines(file.path(root,
  'discovery/source_reviews/ices_ling_faroe_5b_2026.md'), warn = FALSE), collapse = '\n')

expected <- c(
  WGDEEP_2026_v1.pdf = '4dece836b45b56c03a34d2ddf7ee1dd88eda7cdf2ba8d6717f86949888da89cd',
  WGDEEP_2026_ling_5b_section.txt = 'f5f68201da8f1874c901e73ba591dec0cf6bca4fe4ac167c26a9d749f7cc0e89',
  WKBARFAR_2021_report.pdf = '94a98c08f64ca05ea6a5d3db32bba4d29a0cc2e863580619cefd8b181d3b7949',
  ICES_advice_7205.html = 'dd4fc9f818ef445a311e196aa81e6a220efc2727fdf3f0ef89cdcca3b59b6187',
  ICES_standard_graphs_source_23735.html = 'c00314737346f7ca637ccec1bc29375158e08c0544bddb5966b1758882bff6c1',
  ICES_SAG_23735.xml = '992bdf680900755292bf527f5615848947e7c8e71887d42590fdfb917073a5dc'
)
manifest_path <- function(file) file.path(cache, manifest$local_path[match(file, manifest$file)])
stopifnot(
  nrow(candidate) == 1L,
  candidate$official_stock_id == 'lin.27.5b',
  candidate$scientific_name == 'Molva molva',
  candidate$assessment_area == 'Faroes',
  candidate$inventory_assessment_year == 2026L,
  candidate$candidate_unit_id == 'ices_ling_faroe_5b',
  candidate$match_status == 'new',
  candidate$screening_status == 'screened',
  candidate$eligibility == 'eligible',
  candidate$import_status == 'not_started',
  grepl('ages 3–12+', candidate$catch_age_evidence, fixed = TRUE),
  grepl('ages 4–11', candidate$survey_age_evidence, fixed = TRUE),
  grepl('ages 3–11', candidate$survey_age_evidence, fixed = TRUE),
  grepl('2021–2023', candidate$survey_age_evidence, fixed = TRUE),
  grepl('does not publish numeric weights', candidate$survey_age_evidence, fixed = TRUE),
  grepl('does not label the weight unit', candidate$weights_evidence, fixed = TRUE),
  grepl('0.15', candidate$mortality_evidence, fixed = TRUE),
  grepl('February/March', candidate$timing_evidence, fixed = TRUE),
  grepl('August', candidate$timing_evidence, fixed = TRUE),
  grepl('F-at-age', candidate$output_evidence, fixed = TRUE),
  nrow(manifest) == length(expected),
  !anyDuplicated(manifest$file),
  setequal(manifest$file, names(expected)),
  all(manifest$retrieved_on == '2026-10-31'),
  all(file.exists(vapply(manifest$file, manifest_path, character(1L)))),
  all(file.info(vapply(manifest$file, manifest_path, character(1L)))$size == manifest$bytes),
  identical(tolower(manifest$sha256[match(names(expected), manifest$file)]), unname(expected)),
  nrow(unit) == 1L,
  unit$inventory_entries == 1L,
  unit$match_status == 'new',
  unit$screening_status == 'screened',
  unit$eligibility == 'eligible',
  grepl('eligible unit for source curation', source_review, ignore.case = TRUE),
  grepl('verify the survey index scale', source_review, fixed = TRUE),
  !any(tolower(stocks$authority_stock_id) == 'lin.27.5b', na.rm = TRUE),
  !any(assessments$stock_id == 'ices_ling_faroe_5b'),
  !any(inputs$assessment_id == 'ices_ling_faroe_5b'),
  !any(outputs$assessment_id == 'ices_ling_faroe_5b'),
  any(register$candidate_id == 'ICES:lin.27.5a' &
    register$candidate_unit_id == 'ices_ling_5a')
)

cat('ICES 2026 Faroes ling identity, source evidence, provenance, and no-import checks passed.\n')
