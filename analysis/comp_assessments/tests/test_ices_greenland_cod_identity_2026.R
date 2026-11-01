root <- file.path('analysis', 'comp_assessments')
old_dry_run <- Sys.getenv('TAM_SOURCE_DRY_RUN', unset = NA_character_)
Sys.setenv(TAM_SOURCE_DRY_RUN = '1')
review <- new.env(parent = globalenv())
source(file.path(root,
  'scripts/database/490_reconcile_ices_greenland_cod_identities_2026.R'), local = review)
if (is.na(old_dry_run)) Sys.unsetenv('TAM_SOURCE_DRY_RUN') else
  Sys.setenv(TAM_SOURCE_DRY_RUN = old_dry_run)

candidates <- read.csv(file.path(root, 'discovery/candidates.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
units <- read.csv(file.path(root, 'discovery/inventory_units.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
stocks <- read.csv(file.path(root, 'database/stocks.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
assessments <- read.csv(file.path(root, 'database/assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
inputs <- read.csv(file.path(root, 'database/inputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
outputs <- read.csv(file.path(root, 'database/outputs.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'), check.names = FALSE)
ids <- c('ICES:cod.21.1.osc', 'ICES:cod.21.1.isc', 'ICES:cod.21.27.1.14')
greenland <- candidates[match(ids, candidates$candidate_id), , drop = FALSE]
unit_ids <- c('ices_cod_west_greenland_offshore', 'ices_cod_west_greenland_inshore',
  'ices_cod_east_greenland_iceland_offshore')

stopifnot(nrow(greenland) == 3L,
  identical(greenland$candidate_unit_id, unit_ids),
  all(greenland$match_status == 'new'),
  all(greenland$screening_status == 'not_started'),
  all(greenland$eligibility == 'unresolved'),
  all(greenland$review_date == '2026-11-01'),
  all(grepl('Identity reviewed', greenland$decision_reason, fixed = TRUE)),
  all(grepl('eligibility remains unresolved', greenland$notes, fixed = TRUE)),
  all(unit_ids %in% units$candidate_unit_id),
  all(units$inventory_entries[match(unit_ids, units$candidate_unit_id)] == 1L),
  !any(stocks$authority_stock_id %in% c('cod.21.1.osc', 'cod.21.1.isc',
    'cod.21.27.1.14')),
  !any(assessments$stock_id %in% unit_ids),
  !any(inputs$assessment_id %in% unit_ids),
  !any(outputs$assessment_id %in% unit_ids),
  all(file.exists(file.path(root, 'source_cache', 'discovery_ices_greenland_cod_2026',
    c('cod_21_1_osc_advice_2026.html', 'cod_21_1_isc_advice_2026.html',
      'cod_21_27_1_14_advice_2026.html', 'sources.csv')))))

cat('ICES Greenland cod identity, component separation, and no-import checks passed.\n')
