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

lookup <- function(id) register[register$candidate_id == id, , drop = FALSE]
stopifnot(
  lookup('NOAA:10523')$matched_stock_id == 'joint_hake_coastal',
  lookup('DFO:Pacific:Yellowtail Rockfish - British Columbia coastwide')$eligibility == 'unresolved',
  lookup('DFO:Pacific:Yellowtail Rockfish - British Columbia coastwide')$import_status == 'not_started',
  lookup('DFO:Pacific:Pacific Hake - Offshore')$matched_stock_id == 'joint_hake_coastal',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 5Zjm')$matched_stock_id == 'ck_nefsc_trac_cert_morhua_gb',
  lookup('DFO:Pacific:Canary Rockfish')$matched_stock_id == 'dfo_canary_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 3CD-WCVI')$matched_stock_id == 'dfo_pop_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 5ABC-QCS')$matched_stock_id == 'dfo_pop_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 5DE-HS/DE/WHG')$matched_stock_id == 'dfo_pop_bc',
  lookup('NOAA:10014')$matched_stock_id == 'nefsc_witch_flounder_unit',
  lookup('NOAA:10026')$matched_stock_id == 'afsc_alaska_plaice_bsai',
  lookup('NOAA:10031')$matched_stock_id == 'afsc_yellowfin_sole_ebs',
  lookup('NOAA:10191')$matched_stock_id == 'nefsc_summer_flounder',
  lookup('NOAA:11050')$matched_stock_id == 'afsc_northern_rock_sole_ebs',
  lookup('NOAA:12746')$matched_stock_id == 'afsc_cod_ebs',
  lookup('NOAA:10510')$matched_stock_id == 'afsc_pollock_ebs',
  lookup('DFO:Gulf:Atlantic Herring - NAFO 4T (Spring Spawner)')$matched_stock_id == 'dfo_herring_4tvn_spring',
  lookup('DFO:Gulf:Atlantic Herring - NAFO 4T (Fall Spawner)')$matched_stock_id == 'ck_dfo_gsl_harengus_4t_fall',
  lookup('DFO:Gulf:American Plaice - NAFO 4T')$eligibility == 'unresolved',
  lookup('DFO:Pacific:Dungeness Crab')$eligibility == 'ineligible',
  lookup('DFO:Maritimes:Grey Seal')$eligibility == 'ineligible'
)
stopifnot(
  lookup('ICES:nep.fu.14')$eligibility == 'ineligible',
  lookup('ICES:pra.27.1-2')$eligibility == 'ineligible',
  lookup('DFO:Arctic:Lake Trout - Great Slave Lake')$eligibility == 'ineligible',
  lookup('NOAA:10523')$screening_status == 'not_started'
)
stopifnot(
  lookup('NOAA:10809')$eligibility == 'ineligible',
  lookup('NOAA:10977')$eligibility == 'ineligible',
  lookup('NOAA:11763')$eligibility == 'unresolved',
  lookup('NOAA:11763')$screening_status == 'not_started'
)
progress <- read.csv('analysis/comp_assessments/discovery/inventory_progress.csv', stringsAsFactors = FALSE)
stopifnot(sum(progress$inventory_rows) == nrow(register),
  all(progress$inventory_rows == progress$matched + progress$new + progress$ambiguous + progress$identity_unreviewed),
  all(progress$inventory_rows == progress$eligible + progress$excluded + progress$eligibility_unresolved),
  all(progress$inventory_rows == progress$screening_not_started + progress$screening_in_progress + progress$screening_complete))
for (org in unique(register$organization)) {
  stopifnot(sum(progress$inventory_rows[progress$organization == org]) == sum(register$organization == org),
    sum(progress$excluded[progress$organization == org]) == sum(register$organization == org & register$eligibility == 'ineligible'),
    sum(progress$matched[progress$organization == org]) == sum(register$organization == org & register$match_status == 'matched'))
}
cat('Discovery register structure, identity and eligibility checks passed.\n')
