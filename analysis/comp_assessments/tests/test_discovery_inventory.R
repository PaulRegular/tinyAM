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
  lookup('DFO:Pacific:Yellowtail Rockfish - British Columbia coastwide')$screening_status == 'screened',
  lookup('DFO:Pacific:Yellowtail Rockfish - British Columbia coastwide')$import_status == 'not_started',
  lookup('DFO:Newfoundland and Labrador:Atlantic Cod - NAFO 2J3KL')$screening_status == 'screened',
  lookup('DFO:Newfoundland and Labrador:Atlantic Cod - NAFO 2J3KL')$eligibility == 'eligible',
  lookup('DFO:Newfoundland and Labrador:Atlantic Cod - NAFO 2J3KL')$import_status == 'imported',
  lookup('DFO:Quebec:Atlantic Cod - NAFO 3Pn4RS')$screening_status == 'screened',
  lookup('DFO:Quebec:Atlantic Cod - NAFO 3Pn4RS')$eligibility == 'eligible',
  lookup('DFO:Quebec:Atlantic Cod - NAFO 3Pn4RS')$import_status == 'imported',
  lookup('DFO:Pacific:Pacific Hake - Offshore')$matched_stock_id == 'joint_hake_coastal',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 5Zjm')$matched_stock_id == 'ck_nefsc_trac_cert_morhua_gb',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 5Zjm')$screening_status == 'screened',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 5Zjm')$eligibility == 'unresolved',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 5Zjm')$import_status == 'imported',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 4X5Y')$screening_status == 'screened',
  lookup('DFO:Maritimes:Atlantic Cod - NAFO 4X5Y')$eligibility == 'unresolved',
  grepl('run-matched matrices', lookup('DFO:Maritimes:Atlantic Cod - NAFO 4X5Y')$decision_reason, fixed = TRUE),
  lookup('DFO:Pacific:Canary Rockfish')$matched_stock_id == 'dfo_canary_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 3CD-WCVI')$matched_stock_id == 'dfo_pop_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 5ABC-QCS')$matched_stock_id == 'dfo_pop_bc',
  lookup('DFO:Pacific:Pacific Ocean Perch - PMFC 5DE-HS/DE/WHG')$matched_stock_id == 'dfo_pop_bc',
  lookup('NOAA:10014')$matched_stock_id == 'nefsc_witch_flounder_unit',
  lookup('NOAA:10026')$screening_status == 'screened',
  lookup('NOAA:10026')$eligibility == 'ineligible',
  lookup('NOAA:10026')$import_status == 'imported',
  lookup('NOAA:10026')$matched_stock_id == 'afsc_alaska_plaice_bsai',
  lookup('NOAA:10031')$matched_stock_id == 'afsc_yellowfin_sole_ebs',
  lookup('NOAA:10191')$matched_stock_id == 'nefsc_summer_flounder',
  lookup('NOAA:11050')$matched_stock_id == 'afsc_northern_rock_sole_ebs',
  lookup('NOAA:12746')$matched_stock_id == 'afsc_cod_ebs',
  lookup('NOAA:10510')$matched_stock_id == 'afsc_pollock_ebs',
  lookup('DFO:Gulf:Atlantic Herring - NAFO 4T (Spring Spawner)')$matched_stock_id == 'dfo_herring_4tvn_spring',
  lookup('DFO:Gulf:Atlantic Herring - NAFO 4T (Fall Spawner)')$matched_stock_id == 'ck_dfo_gsl_harengus_4t_fall',
  lookup('DFO:Gulf:American Plaice - NAFO 4T')$eligibility == 'eligible',
  lookup('DFO:Pacific:Dungeness Crab')$eligibility == 'ineligible',
  lookup('DFO:Maritimes:Grey Seal')$eligibility == 'ineligible'
)
stopifnot(
  lookup('DFO:Maritimes:Atlantic Canada Spiny Dogfish - NAFO 4VWX5')$candidate_unit_id == 'spiny_dogfish_northwest_atlantic',
  lookup('DFO:Maritimes:Atlantic Canada Spiny Dogfish - NAFO 4VWX5')$match_status == 'ambiguous',
  lookup('DFO:Maritimes:Atlantic Canada Spiny Dogfish - NAFO 4VWX5')$eligibility == 'ineligible',
  lookup('NOAA:10701')$candidate_unit_id == 'spiny_dogfish_northwest_atlantic',
  lookup('NOAA:10701')$match_status == 'new',
  lookup('NOAA:10701')$screening_status == 'screened',
  lookup('NOAA:10701')$eligibility == 'ineligible',
  grepl('age and growth data', lookup('NOAA:10701')$decision_reason, ignore.case = TRUE),
  file.exists('analysis/comp_assessments/discovery/source_reviews/noaa_dfo_spiny_dogfish_northwest_atlantic.md'),
  lookup('ICES:nop.27.6a')$candidate_unit_id == 'ices_norway_pout_6a',
  lookup('ICES:nop.27.6a')$match_status == 'new',
  lookup('ICES:nop.27.6a')$screening_status == 'screened',
  lookup('ICES:nop.27.6a')$eligibility == 'ineligible',
  lookup('ICES:nop.27.6a')$import_status == 'not_started',
  grepl('status only', lookup('ICES:nop.27.6a')$output_evidence, ignore.case = TRUE),
  lookup('ICES:nep.fu.14')$eligibility == 'ineligible',
  lookup('ICES:wit.27.3a47d')$match_status == 'new',
  lookup('ICES:wit.27.3a47d')$screening_status == 'screened',
  lookup('ICES:wit.27.3a47d')$eligibility == 'unresolved',
  lookup('ICES:sol.27.7fg')$match_status == 'matched',
  lookup('ICES:sol.27.7fg')$screening_status == 'screened',
  lookup('ICES:sol.27.7fg')$eligibility == 'eligible',
  lookup('ICES:sol.27.7fg')$import_status == 'imported',
  lookup('ICES:sol.27.7e')$match_status == 'matched',
  lookup('ICES:sol.27.7e')$screening_status == 'screened',
  lookup('ICES:sol.27.7e')$eligibility == 'eligible',
  lookup('ICES:sol.27.7e')$import_status == 'imported',
  lookup('ICES:cod.27.7a')$match_status == 'matched',
  lookup('ICES:cod.27.7a')$screening_status == 'screened',
  lookup('ICES:cod.27.7a')$eligibility == 'eligible',
  lookup('ICES:cod.27.7a')$import_status == 'imported',
  lookup('ICES:had.27.6b')$match_status == 'matched',
  lookup('ICES:had.27.6b')$screening_status == 'screened',
  lookup('ICES:had.27.6b')$eligibility == 'eligible',
  lookup('ICES:had.27.6b')$import_status == 'imported',
  lookup('ICES:whg.27.6a')$match_status == 'matched',
  lookup('ICES:whg.27.6a')$screening_status == 'screened',
  lookup('ICES:whg.27.6a')$eligibility == 'eligible',
  lookup('ICES:whg.27.6a')$import_status == 'imported',
  lookup('ICES:cod.27.7e-k')$match_status == 'matched',
  lookup('ICES:cod.27.7e-k')$screening_status == 'screened',
  lookup('ICES:cod.27.7e-k')$eligibility == 'eligible',
  lookup('ICES:cod.27.7e-k')$import_status == 'imported',
  lookup('ICES:had.27.7b-k')$match_status == 'matched',
  lookup('ICES:had.27.7b-k')$screening_status == 'screened',
  lookup('ICES:had.27.7b-k')$eligibility == 'eligible',
  lookup('ICES:had.27.7b-k')$import_status == 'imported',
  lookup('ICES:ple.27.7a')$match_status == 'matched',
  lookup('ICES:ple.27.7a')$screening_status == 'screened',
  lookup('ICES:ple.27.7a')$eligibility == 'eligible',
  lookup('ICES:ple.27.7a')$import_status == 'imported',
  grepl('1981-2025', lookup('ICES:ple.27.7a')$catch_age_evidence, fixed = TRUE),
  grepl('UK(E&W) BTS-Q3', lookup('ICES:ple.27.7a')$survey_age_evidence, fixed = TRUE),
  lookup('ICES:sol.27.7a')$match_status == 'new',
  lookup('ICES:sol.27.7a')$screening_status == 'screened',
  lookup('ICES:sol.27.7a')$eligibility == 'unresolved',
  lookup('ICES:sol.27.8ab')$screening_status == 'screened',
  lookup('ICES:sol.27.8ab')$eligibility == 'unresolved',
  grepl('catch scale', lookup('ICES:sol.27.8ab')$decision_reason, fixed = TRUE),
  grepl('production-run M', lookup('ICES:sol.27.8ab')$decision_reason, fixed = TRUE),
  lookup('ICES:meg.27.8c9a')$screening_status == 'screened',
  lookup('ICES:meg.27.8c9a')$eligibility == 'unresolved',
  all(vapply(c('lez.27.4a6a', 'lez.27.6b', 'ple.27.7e', 'ple.27.7fg',
    'ple.27.7h-k', 'had.27.7a', 'sol.27.7h-k', 'bss.27.6a7bj',
    'cod.27.6b', 'ple.27.7bc', 'sol.27.7bc'), function(code) {
      row <- lookup(paste0('ICES:', code))
      nrow(row) == 1L && row$screening_status == 'screened' && row$eligibility == 'ineligible'
    }, logical(1L))),
  lookup('ICES:sol.27.8c9a')$screening_status == 'screened',
  lookup('ICES:sol.27.8c9a')$eligibility == 'ineligible',
  grepl('latest WGBIE update is length-based', lookup('ICES:sol.27.8c9a')$decision_reason, fixed = TRUE),
  lookup('ICES:meg.27.8c9a')$import_status == 'not_started',
  lookup('ICES:mac.27.nea')$match_status == 'new',
  lookup('ICES:mac.27.nea')$screening_status == 'screened',
  lookup('ICES:mac.27.nea')$eligibility == 'unresolved',
  grepl('age-specific fitted outputs', lookup('ICES:mac.27.nea')$decision_reason, fixed = TRUE),
  lookup('ICES:pra.27.1-2')$eligibility == 'ineligible',
  lookup('DFO:Arctic:Lake Trout - Great Slave Lake')$eligibility == 'ineligible',
  lookup('NOAA:10523')$screening_status == 'screened'
)
stopifnot(
  lookup('NOAA:10809')$eligibility == 'ineligible',
  lookup('NOAA:10011')$screening_status == 'screened',
  lookup('NOAA:10011')$eligibility == 'eligible',
  lookup('NOAA:10011')$import_status == 'not_started',
  grepl('Catch-at-age can be reconstructed', lookup('NOAA:10011')$catch_age_evidence, fixed = TRUE),
  grepl('Triennial survey has no age samples', lookup('NOAA:10011')$survey_age_evidence, fixed = TRUE),
  grepl('2023 spawning output', lookup('NOAA:10011')$output_evidence, fixed = TRUE),
  lookup('NOAA:17357')$match_status == 'new',
  lookup('NOAA:17357')$screening_status == 'screened',
  lookup('NOAA:17357')$eligibility == 'eligible',
  lookup('NOAA:17357')$import_status == 'not_started',
  grepl('Catch-at-age is not direct', lookup('NOAA:17357')$catch_age_evidence, fixed = TRUE),
  grepl('age-specific index requires', lookup('NOAA:17357')$survey_age_evidence, fixed = TRUE),
  grepl('N-at-age', lookup('NOAA:17357')$output_evidence, fixed = TRUE),
  lookup('NOAA:10395')$match_status == 'ambiguous',
  lookup('NOAA:10395')$screening_status == 'screened',
  lookup('NOAA:10395')$eligibility == 'ineligible',
  grepl('NOAA:17357', lookup('NOAA:10395')$notes, fixed = TRUE),
  lookup('NOAA:10977')$eligibility == 'ineligible',
  lookup('NOAA:11763')$eligibility == 'unresolved',
  lookup('NOAA:11763')$screening_status == 'not_started'
)
matched_noaa_review <- c(
  'NOAA:10014','NOAA:10019','NOAA:10029','NOAA:10031','NOAA:10033','NOAA:10034',
  'NOAA:10191','NOAA:10286','NOAA:10506','NOAA:10510','NOAA:10511','NOAA:10513',
  'NOAA:10515','NOAA:10518','NOAA:10520','NOAA:10523','NOAA:11050','NOAA:11055',
  'NOAA:11662','NOAA:12745','NOAA:12746'
)
ne_noaa_review <- c('NOAA:10164','NOAA:10508','NOAA:10509','NOAA:10519')
reviewed_ne_noaa <- register[match(ne_noaa_review, register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(reviewed_ne_noaa$candidate_id),
  identical(reviewed_ne_noaa$matched_stock_id, c('nefsc_atlantic_mackerel',
    'ck_nefsc_saw_morhua_gom','nefsc_cod_georges_bank','ck_nefsc_garmiii_aeglefinus_5y')),
  all(reviewed_ne_noaa$match_status == 'matched'),
  all(reviewed_ne_noaa$screening_status == 'screened'),
  all(reviewed_ne_noaa$import_status == 'imported'),
  identical(reviewed_ne_noaa$eligibility, c('eligible','ineligible','unresolved','unresolved')),
  grepl('combined U.S.-Canadian', lookup('NOAA:10164')$catch_age_evidence, fixed = TRUE),
  grepl('successor stock units', lookup('NOAA:10508')$decision_reason, fixed = TRUE),
  grepl('10509 / assessment 13596', lookup('NOAA:10509')$decision_reason, fixed = TRUE),
  grepl('2026 assessment', lookup('NOAA:10519')$decision_reason, fixed = TRUE)
)
reviewed_noaa_cod_units <- register[match(c('NOAA:17280', 'NOAA:17281', 'NOAA:17282'), register$candidate_id), , drop = FALSE]
stopifnot(
  nrow(reviewed_noaa_cod_units) == 3L,
  all(reviewed_noaa_cod_units$match_status == 'new'),
  all(is.na(reviewed_noaa_cod_units$matched_stock_id)),
  all(reviewed_noaa_cod_units$screening_status == 'screened'),
  all(reviewed_noaa_cod_units$eligibility == 'unresolved'),
  all(reviewed_noaa_cod_units$import_status == 'not_started'),
  all(grepl('2024', reviewed_noaa_cod_units$detailed_source_url, fixed = TRUE)),
  grepl('commercial and recreational catch-age composition', lookup('NOAA:17280')$catch_age_evidence, fixed = TRUE),
  grepl('since 2013', lookup('NOAA:17281')$catch_age_evidence, fixed = TRUE),
  grepl('no ages from the recreational fishery', lookup('NOAA:17282')$catch_age_evidence, fixed = TRUE),
  all(grepl('not', reviewed_noaa_cod_units$output_evidence, fixed = TRUE)),
  lookup('NOAA:12805')$matched_stock_id == 'ck_nefsc_trac_cert_morhua_gb',
  lookup('NOAA:12805')$match_status == 'matched',
  lookup('NOAA:12805')$screening_status == 'screened',
  lookup('NOAA:12805')$eligibility == 'unresolved',
  grepl('no assessment model', lookup('NOAA:12805')$catch_age_evidence, fixed = TRUE),
  grepl('2023 TRAC document is an indicator-based status report', lookup('NOAA:12805')$decision_reason, fixed = TRUE)
)
reviewed_noaa <- register[match(matched_noaa_review, register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(reviewed_noaa$candidate_id),
  all(reviewed_noaa$match_status == 'matched'),
  all(reviewed_noaa$screening_status == 'screened'),
  all(reviewed_noaa$import_status == 'imported'),
  identical(reviewed_noaa$eligibility, c(
    'ineligible','unresolved','eligible','eligible','unresolved','unresolved',
    'eligible','eligible','ineligible','eligible','eligible','eligible','eligible',
    'eligible','unresolved','eligible','eligible','eligible','eligible','ineligible','ineligible')),
  grepl('length compositions', lookup('NOAA:10506')$catch_age_evidence, fixed = TRUE),
  grepl('current official', lookup('NOAA:10014')$decision_reason, fixed = TRUE),
  grepl('age-based-M production run', lookup('NOAA:10033')$decision_reason, fixed = TRUE),
  grepl('2018 official ASAP record', lookup('NOAA:10191')$decision_reason, fixed = TRUE),
  grepl('no unconditional age-specific index', lookup('NOAA:12745')$survey_age_evidence, fixed = TRUE),
  all(vapply(c('catch_age_evidence','survey_age_evidence','weights_evidence',
    'maturity_evidence','mortality_evidence','timing_evidence','output_evidence',
    'detailed_source_url'), function(column) {
      values <- reviewed_noaa[[column]][reviewed_noaa$eligibility == 'eligible']
      all(!is.na(values) & nzchar(values))
    }, logical(1L)))
)
progress <- read.csv('analysis/comp_assessments/discovery/inventory_progress.csv', stringsAsFactors = FALSE)
stopifnot(sum(progress$inventory_rows) == nrow(register),
  all(progress$inventory_rows == progress$matched + progress$new + progress$ambiguous + progress$identity_unreviewed),
  all(progress$inventory_rows == progress$eligible + progress$excluded + progress$eligibility_unresolved),
  all(progress$inventory_rows == progress$screening_not_started + progress$screening_in_progress + progress$screening_complete))
stopifnot(all(progress$eligible_new <= progress$new),
  all(progress$eligible_new <= progress$eligible))
dfo_progress <- progress[progress$organization == 'DFO', , drop = FALSE]
stopifnot(sum(dfo_progress$inventory_rows) == 199L,
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$ambiguous) == 5L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$eligibility_unresolved) == 69L)
unit3_redfish <- lookup('DFO:Maritimes:Acadian Redfish - Unit 3')
stopifnot(
  nrow(unit3_redfish) == 1L,
  unit3_redfish$scientific_name == 'Sebastes spp.',
  unit3_redfish$match_status == 'new',
  is.na(unit3_redfish$matched_stock_id),
  unit3_redfish$screening_status == 'screened',
  unit3_redfish$eligibility == 'ineligible',
  unit3_redfish$import_status == 'not_started',
  grepl('not catch-at-age', unit3_redfish$catch_age_evidence, fixed = TRUE),
  grepl('abundance-at-length', unit3_redfish$survey_age_evidence, fixed = TRUE),
  grepl('above 22 cm', unit3_redfish$maturity_evidence, fixed = TRUE),
  grepl('not separated in commercial or survey catches', unit3_redfish$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_unit3_redfish_2022.md'))
)
noaa_progress <- progress[progress$organization == 'NOAA', , drop = FALSE]
stopifnot(sum(noaa_progress$new) == 50L, sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$ambiguous) == 3L, sum(noaa_progress$identity_unreviewed) == 351L,
  sum(noaa_progress$eligible) == 37L, sum(noaa_progress$excluded) == 63L,
  sum(noaa_progress$eligibility_unresolved) == 332L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L)
halibut <- lookup('NOAA:10056')
stopifnot(
  nrow(halibut) == 1L,
  halibut$match_status == 'ambiguous',
  is.na(halibut$matched_stock_id),
  halibut$screening_status == 'screened',
  halibut$eligibility == 'ineligible',
  halibut$import_status == 'not_started',
  grepl('aggregate annual catch', halibut$catch_age_evidence, fixed = TRUE),
  grepl('not age-specific', halibut$survey_age_evidence, fixed = TRUE),
  grepl('does not use maturity-at-age', halibut$maturity_evidence, fixed = TRUE),
  grepl('does not estimate a natural-mortality process', halibut$mortality_evidence, fixed = TRUE),
  grepl('does not estimate N-', halibut$output_evidence, fixed = TRUE),
  grepl('DFO halibut unit', halibut$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'noaa_atlantic_halibut_2024.md'))
)
greenland_turbot <- lookup('NOAA:10054')
stopifnot(
  nrow(greenland_turbot) == 1L,
  greenland_turbot$match_status == 'new',
  is.na(greenland_turbot$matched_stock_id),
  greenland_turbot$screening_status == 'screened',
  greenland_turbot$eligibility == 'unresolved',
  greenland_turbot$import_status == 'not_started',
  grepl('aggregate catch', greenland_turbot$catch_age_evidence, fixed = TRUE),
  grepl('excluded from the likelihood', greenland_turbot$survey_age_evidence, fixed = TRUE),
  grepl('0.112 per year', greenland_turbot$mortality_evidence, fixed = TRUE),
  grepl('60 cm', greenland_turbot$maturity_evidence, fixed = TRUE),
  grepl('numbers-at-age', greenland_turbot$output_evidence, fixed = TRUE),
  grepl('preliminary and not treated', greenland_turbot$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'noaa_greenland_turbot_bsai_2024.md'))
)
for (org in unique(register$organization)) {
  stopifnot(sum(progress$inventory_rows[progress$organization == org]) == sum(register$organization == org),
    sum(progress$excluded[progress$organization == org]) == sum(register$organization == org & register$eligibility == 'ineligible'),
    sum(progress$matched[progress$organization == org]) == sum(register$organization == org & register$match_status == 'matched'))
}
cat('Discovery register structure, identity and eligibility checks passed.\n')

wgcse <- function(code) lookup(paste0('ICES:', code))
stopifnot(
  wgcse('anf.27.3a46')$match_status == 'new',
  wgcse('anf.27.3a46')$screening_status == 'screened',
  wgcse('anf.27.3a46')$eligibility == 'ineligible',
  grepl('length-based', wgcse('anf.27.3a46')$decision_reason, fixed = TRUE),
  wgcse('bss.27.4bc7ad-h')$match_status == 'new',
  wgcse('bss.27.4bc7ad-h')$screening_status == 'screened',
  wgcse('bss.27.4bc7ad-h')$eligibility == 'unresolved',
  grepl('conditional age-at-length', wgcse('bss.27.4bc7ad-h')$decision_reason, fixed = TRUE),
  wgcse('pol.27.67')$match_status == 'new',
  wgcse('pol.27.67')$screening_status == 'screened',
  wgcse('pol.27.67')$eligibility == 'eligible',
  wgcse('pol.27.67')$import_status == 'not_started',
  grepl('since 2010', wgcse('pol.27.67')$catch_age_evidence, fixed = TRUE),
  grepl('0.920', wgcse('pol.27.67')$mortality_evidence, fixed = TRUE),
  wgcse('whg.27.7a')$match_status == 'matched',
  wgcse('whg.27.7a')$matched_stock_id == 'ck_ices_wgcse_merlangus_cs7a',
  wgcse('whg.27.7a')$eligibility == 'eligible',
  wgcse('whg.27.7a')$import_status == 'imported',
  wgcse('whg.27.7b-ce-k')$match_status == 'new',
  wgcse('whg.27.7b-ce-k')$screening_status == 'screened',
  wgcse('whg.27.7b-ce-k')$eligibility == 'eligible',
  wgcse('whg.27.7b-ce-k')$import_status == 'not_started',
  grepl('0-7+', wgcse('whg.27.7b-ce-k')$catch_age_evidence, fixed = TRUE),
  grepl('ages 0-2', wgcse('whg.27.7b-ce-k')$survey_age_evidence, fixed = TRUE),
  is.na(wgcse('whg.27.7b-ce-k')$matched_stock_id)
)
cat('WGCSE age-based screening decisions passed.\n')

alaska_flatfish <- register[match(paste0('NOAA:', c('10002','10003','10004','10015','10016')), register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(alaska_flatfish$candidate_id),
  identical(alaska_flatfish$match_status, c('new','matched','new','new','new')),
  identical(alaska_flatfish$matched_stock_id, c(NA, 'afsc_arrowtooth_flounder_bsai', NA, NA, NA)),
  identical(alaska_flatfish$screening_status, rep('screened', 5L)),
  identical(alaska_flatfish$eligibility, c('ineligible','ineligible','ineligible','eligible','ineligible')),
  identical(alaska_flatfish$import_status, c('not_started','imported','not_started','not_started','not_started')),
  grepl('does not report a fishery catch-age', lookup('NOAA:10002')$catch_age_evidence, fixed = TRUE),
  grepl('no accepted fishery age-composition', lookup('NOAA:10003')$catch_age_evidence, fixed = TRUE),
  grepl('Stock SMART lists Reinhardtius stomias', lookup('NOAA:10003')$notes, fixed = TRUE),
  grepl('Stock SMART snapshot says 2024 and lists Reinhardtius stomias', lookup('NOAA:10004')$notes, fixed = TRUE),
  grepl('fishery length compositions', lookup('NOAA:10004')$catch_age_evidence, fixed = TRUE),
  grepl('October 2026 SSC report selects Model 25.0', lookup('NOAA:10004')$output_evidence, fixed = TRUE),
  grepl('age compositions', lookup('NOAA:10015')$catch_age_evidence, fixed = TRUE),
  grepl('flathead sole only', lookup('NOAA:10015')$catch_age_evidence, fixed = TRUE),
  grepl('ghosted', lookup('NOAA:10016')$survey_age_evidence, fixed = TRUE),
  all(file.exists(file.path('analysis','comp_assessments','discovery','source_reviews',
    c('noaa_kamchatka_bsai.md','noaa_arrowtooth_bsai_identity.md','noaa_arrowtooth_goa.md',
      'noaa_flathead_bsai.md','noaa_flathead_goa.md'))))
)
cat('Alaska flatfish identity and age-input screening checks passed.\n')

goa_pop <- lookup('NOAA:11056')
stopifnot(
  goa_pop$match_status == 'new', is.na(goa_pop$matched_stock_id),
  goa_pop$screening_status == 'screened', goa_pop$eligibility == 'eligible',
  goa_pop$import_status == 'not_started',
  grepl('1961-2023', goa_pop$catch_age_evidence, fixed = TRUE),
  grepl('17 sampled years', goa_pop$catch_age_evidence, fixed = TRUE),
  grepl('14 years', goa_pop$survey_age_evidence, fixed = TRUE),
  grepl('fishery and survey', goa_pop$weights_evidence, fixed = TRUE),
  grepl('one constant parameter', goa_pop$mortality_evidence, fixed = TRUE),
  grepl('not survey timing', goa_pop$timing_evidence, fixed = TRUE),
  grepl('numbers-at-age', goa_pop$output_evidence, fixed = TRUE),
  grepl('2023 operational update', goa_pop$notes, fixed = TRUE),
  file.exists(file.path('analysis','comp_assessments','discovery','source_reviews',
    'noaa_goa_pacific_ocean_perch.md'))
)
cat('GOA Pacific ocean perch identity and source-screen checks passed.\n')

northern_rockfish <- register[match(paste0('NOAA:', c('10416', '10417')), register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(northern_rockfish$candidate_id),
  identical(northern_rockfish$match_status, c('new', 'new')),
  all(is.na(northern_rockfish$matched_stock_id)),
  identical(northern_rockfish$screening_status, c('screened', 'screened')),
  identical(northern_rockfish$eligibility, c('eligible', 'eligible')),
  identical(northern_rockfish$import_status, c('not_started', 'not_started')),
  grepl('2023 full assessment', lookup('NOAA:10416')$catch_age_evidence, fixed = TRUE),
  grepl('2024 full assessment', lookup('NOAA:10417')$catch_age_evidence, fixed = TRUE),
  grepl('2026 chapter is a harvest projection', lookup('NOAA:10416')$notes, fixed = TRUE),
  grepl('December 2026', lookup('NOAA:10417')$notes, fixed = TRUE),
  grepl('documented reconstruction', lookup('NOAA:10416')$decision_reason, fixed = TRUE),
  grepl('documented reconstruction', lookup('NOAA:10417')$decision_reason, fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    c('noaa_northern_rockfish_bsai.md', 'noaa_northern_rockfish_goa.md'))))
)
cat('Northern rockfish identity and source-screening checks passed.\n')

arrowtooth <- lookup('NOAA:10005')
stopifnot(
  arrowtooth$match_status == 'new',
  arrowtooth$screening_status == 'screened',
  arrowtooth$eligibility == 'unresolved',
  arrowtooth$scientific_name == 'Atheresthes stomias',
  grepl('could not be recovered', arrowtooth$decision_reason, fixed = TRUE),
  grepl('Stock SMART ID 10005', arrowtooth$notes, fixed = TRUE),
  !identical(arrowtooth$matched_stock_id, 'afsc_arrowtooth_flounder_bsai')
)

dover <- lookup('NOAA:10022')
stopifnot(
  dover$match_status == 'new',
  dover$screening_status == 'screened',
  dover$eligibility == 'eligible',
  dover$import_status == 'not_started',
  grepl('year/sex/fleet age-composition', dover$catch_age_evidence, fixed = TRUE),
  grepl('conditional age-at-length', dover$survey_age_evidence, fixed = TRUE),
  grepl('numbers-at-age, catch-at-age, F-at-age', dover$output_evidence, fixed = TRUE),
  grepl('dover sole', dover$decision_reason, ignore.case = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'noaa_pacific_coast_dover_sole.md'))
)
cat('Dover sole identity and source-screening checks passed.\n')

widow <- lookup('NOAA:10568')
stopifnot(
  widow$match_status == 'new',
  widow$screening_status == 'screened',
  widow$eligibility == 'eligible',
  widow$import_status == 'not_started',
  grepl('three commercial fleets', widow$catch_age_evidence, fixed = TRUE),
  grepl('conditional age-at-length', widow$survey_age_evidence, fixed = TRUE),
  grepl('billions of eggs', widow$output_evidence, fixed = TRUE),
  grepl('version labels', widow$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'noaa_west_coast_widow_rockfish.md'))
)
cat('Widow rockfish identity and source-screening checks passed.\n')

dfo_sablefish <- lookup('DFO:Pacific:Sablefish')
stopifnot(
  dfo_sablefish$scientific_name == 'Anoplopoma fimbria',
  dfo_sablefish$assessment_area == 'British Columbia coastwide',
  dfo_sablefish$match_status == 'new',
  is.na(dfo_sablefish$matched_stock_id),
  dfo_sablefish$screening_status == 'screened',
  dfo_sablefish$eligibility == 'unresolved',
  dfo_sablefish$import_status == 'not_started',
  grepl('no public numeric commercial age-composition arrays', dfo_sablefish$catch_age_evidence, fixed = TRUE),
  grepl('2003-2021', dfo_sablefish$survey_age_evidence, fixed = TRUE),
  grepl('prior (mean 0.1 and standard deviation 0.1', dfo_sablefish$mortality_evidence, fixed = TRUE),
  grepl('0.75', dfo_sablefish$timing_evidence, fixed = TRUE),
  grepl('N-at-age and F-at-age output matrices', dfo_sablefish$output_evidence, fixed = TRUE),
  grepl('2025-26 advice', dfo_sablefish$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_bc_sablefish.md'))
)
cat('B.C. Sablefish identity and source-screening checks passed.\n')

dfo_bocaccio <- lookup('DFO:Pacific:Bocaccio Rockfish')
stopifnot(
  dfo_bocaccio$scientific_name == 'Sebastes paucispinis',
  dfo_bocaccio$assessment_area == 'British Columbia coastwide (PMFC 3CD+5ABCDE)',
  dfo_bocaccio$match_status == 'new',
  is.na(dfo_bocaccio$matched_stock_id),
  dfo_bocaccio$screening_status == 'screened',
  dfo_bocaccio$eligibility == 'unresolved',
  dfo_bocaccio$import_status == 'not_started',
  grepl('only the combined trawl fishery has age-frequency observations', dfo_bocaccio$catch_age_evidence, fixed = TRUE),
  grepl('2021 QCS', dfo_bocaccio$survey_age_evidence, fixed = TRUE),
  grepl('0.07, 0.08, and 0.09', dfo_bocaccio$mortality_evidence, fixed = TRUE),
  grepl('fraction 0.5', dfo_bocaccio$timing_evidence, fixed = TRUE),
  grepl('B2024 = 23,789 tonnes', dfo_bocaccio$output_evidence, fixed = TRUE),
  grepl('no age composition for the non-trawl catch stream', dfo_bocaccio$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_bc_bocaccio.md'))
)
cat('B.C. Bocaccio identity and source-screening checks passed.\n')

dfo_sog_herring <- lookup('DFO:Pacific:Pacific Herring - Strait of Georgia')
stopifnot(
  dfo_sog_herring$scientific_name == 'Clupea pallasii',
  dfo_sog_herring$matched_stock_id == 'dfo_herring_sog',
  dfo_sog_herring$match_status == 'matched',
  dfo_sog_herring$screening_status == 'screened',
  dfo_sog_herring$eligibility == 'ineligible',
  dfo_sog_herring$import_status == 'imported',
  grepl('ages 1-10+', dfo_sog_herring$catch_age_evidence, fixed = TRUE),
  grepl('All native survey age samples are missing (-1)', dfo_sog_herring$survey_age_evidence, fixed = TRUE),
  grepl('density dependent', dfo_sog_herring$mortality_evidence, fixed = TRUE),
  grepl('0.766', dfo_sog_herring$timing_evidence, fixed = TRUE),
  grepl('84.33375 thousand tonnes', dfo_sog_herring$output_evidence, fixed = TRUE),
  grepl('cannot fit this aggregate spawning-biomass observation', dfo_sog_herring$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_herring_sog.md'))
)
cat('Strait of Georgia herring identity and source-screening checks passed.\n')

dfo_pop <- register[register$candidate_id %in% c(
  'DFO:Pacific:Pacific Ocean Perch - PMFC 3CD-WCVI',
  'DFO:Pacific:Pacific Ocean Perch - PMFC 5ABC-QCS',
  'DFO:Pacific:Pacific Ocean Perch - PMFC 5DE-HS/DE/WHG'), , drop = FALSE]
stopifnot(
  nrow(dfo_pop) == 3L,
  all(dfo_pop$scientific_name == 'Sebastes alutus'),
  all(dfo_pop$matched_stock_id == 'dfo_pop_bc'),
  all(dfo_pop$match_status == 'matched'),
  all(dfo_pop$screening_status == 'screened'),
  all(dfo_pop$eligibility == 'eligible'),
  all(dfo_pop$import_status == 'imported'),
  all(grepl('three subarea fisheries', dfo_pop$catch_age_evidence, fixed = TRUE)),
  all(grepl('five surveys', dfo_pop$survey_age_evidence, fixed = TRUE)),
  all(grepl('5/12', dfo_pop$timing_evidence, fixed = TRUE)),
  all(grepl('2024 spawning biomass', dfo_pop$output_evidence, fixed = TRUE)),
  all(grepl('do not create duplicate', dfo_pop$notes, fixed = TRUE)),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews', 'dfo_pop_bc.md'))
)
cat('British Columbia Pacific ocean perch identity and source-screening checks passed.\n')

dfo_canary <- lookup('DFO:Pacific:Canary Rockfish')
stopifnot(
  dfo_canary$scientific_name == 'Sebastes pinniger',
  dfo_canary$matched_stock_id == 'dfo_canary_bc',
  dfo_canary$match_status == 'matched',
  dfo_canary$screening_status == 'screened',
  dfo_canary$eligibility == 'eligible',
  dfo_canary$import_status == 'imported',
  grepl('36 years (1977-2017)', dfo_canary$catch_age_evidence, fixed = TRUE),
  grepl('QCS, WCVI and NMFS have source age compositions', dfo_canary$survey_age_evidence, fixed = TRUE),
  grepl('5/12', dfo_canary$timing_evidence, fixed = TRUE),
  grepl('500 official posterior-subset', dfo_canary$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews', 'dfo_canary_bc.md'))
)
cat('B.C. Canary Rockfish identity and source-screening checks passed.\n')

dfo_herring <- register[match(c(
  'DFO:Gulf:Atlantic Herring - NAFO 4T (Spring Spawner)',
  'DFO:Gulf:Atlantic Herring - NAFO 4T (Fall Spawner)'), register$candidate_id), , drop = FALSE]
stopifnot(
  nrow(dfo_herring) == 2L,
  all(dfo_herring$scientific_name == 'Clupea harengus'),
  identical(dfo_herring$matched_stock_id, c(
    'dfo_herring_4tvn_spring', 'ck_dfo_gsl_harengus_4t_fall')),
  all(dfo_herring$match_status == 'matched'),
  all(dfo_herring$screening_status == 'screened'),
  all(dfo_herring$eligibility == 'eligible'),
  all(dfo_herring$import_status == 'imported'),
  grepl('thousands of fish', dfo_herring$survey_age_evidence[1L], fixed = TRUE),
  grepl('duplicat', dfo_herring$catch_age_evidence[2L]),
  grepl('2026', dfo_herring$decision_reason[1L], fixed = TRUE),
  grepl('2026 advisory', dfo_herring$decision_reason[2L], fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_herring_4tvn_spring.md')),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_herring_4tvn_fall.md'))
)
cat('Southern Gulf herring identity and source-screening checks passed.\n')

goa_dover <- lookup('NOAA:10021')
stopifnot(
  goa_dover$scientific_name == 'Microstomus pacificus',
  goa_dover$match_status == 'new', is.na(goa_dover$matched_stock_id),
  goa_dover$screening_status == 'screened',
  goa_dover$eligibility == 'ineligible',
  goa_dover$import_status == 'not_started',
  grepl('no fishery age composition', goa_dover$catch_age_evidence, fixed = TRUE),
  grepl('conditional age-at-length', goa_dover$survey_age_evidence, fixed = TRUE),
  grepl('survey month 6', goa_dover$timing_evidence, fixed = TRUE),
  grepl('numbers-at-age time series', goa_dover$output_evidence, fixed = TRUE),
  grepl('Model 19.3.1', goa_dover$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_goa_dover_sole.md'))
)
cat('NOAA GOA Dover sole source-screening checks passed.\n')

starry <- register[match(c('NOAA:10025', 'NOAA:14804', 'NOAA:14805'), register$candidate_id), , drop = FALSE]
stopifnot(
  nrow(starry) == 3L,
  starry$match_status[1L] == 'ambiguous',
  starry$eligibility[1L] == 'ineligible',
  all(starry$screening_status == 'screened'),
  all(starry$eligibility == 'ineligible'),
  all(starry$match_status[-1L] == 'new'),
  all(starry$inventory_assessment_year[-1L] == 2017L),
  all(grepl('DB-SRA', starry$decision_reason[-1L], fixed = TRUE)),
  all(grepl('2017 DB-SRA', starry$notes[-1L], fixed = TRUE)),
  all(grepl('catch numbers-at-age', starry$catch_age_evidence[-1L], fixed = TRUE)),
  all(grepl('age compositions', starry$survey_age_evidence[-1L], fixed = TRUE)),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_starry_flounder.md'))
)
cat('NOAA starry flounder identity and age-data screening checks passed.\n')

winter <- register[match(c('NOAA:10028', 'NOAA:10030'), register$candidate_id), , drop = FALSE]
stopifnot(
  nrow(winter) == 2L,
  all(winter$scientific_name == 'Pseudopleuronectes americanus'),
  all(winter$match_status == 'new'),
  all(winter$screening_status == 'screened'),
  all(winter$inventory_assessment_year == 2025L),
  winter$eligibility[winter$official_stock_id == '10028'] == 'ineligible',
  grepl('30+ cm', winter$survey_age_evidence[winter$official_stock_id == '10028'], fixed = TRUE),
  winter$eligibility[winter$official_stock_id == '10030'] == 'unresolved',
  grepl('ASAP', winter$decision_reason[winter$official_stock_id == '10030'], fixed = TRUE),
  grepl('catch-age', winter$decision_reason[winter$official_stock_id == '10030'], fixed = TRUE),
  grepl('6,211 thousand fish', winter$output_evidence[winter$official_stock_id == '10030'], fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_winter_flounder_gom_sne.md'))
)
cat('NOAA winter flounder identity and source-screening checks passed.\n')

sablefish <- lookup('NOAA:10466')
stopifnot(
  sablefish$scientific_name == 'Anoplopoma fimbria',
  sablefish$assessment_area == 'Pacific Coast',
  sablefish$inventory_assessment_year == 2025L,
  sablefish$match_status == 'new', is.na(sablefish$matched_stock_id),
  sablefish$screening_status == 'screened',
  sablefish$eligibility == 'eligible',
  sablefish$import_status == 'not_started',
  grepl('1983', sablefish$catch_age_evidence, fixed = TRUE),
  grepl('length-age key', sablefish$catch_age_evidence, fixed = TRUE),
  grepl('no directly observed age-specific abundance', sablefish$survey_age_evidence, fixed = TRUE),
  grepl('0-30', sablefish$weights_evidence, fixed = TRUE),
  grepl('0.088 per year', sablefish$mortality_evidence, fixed = TRUE),
  grepl('July', sablefish$timing_evidence, fixed = TRUE),
  grepl('59,331-89,486 mt', sablefish$output_evidence, fixed = TRUE),
  grepl('3.30.23.1', sablefish$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_west_coast_sablefish_2025.md'))
)
cat('NOAA West Coast Sablefish identity and source-screening checks passed.\n')

yelloweye <- lookup('NOAA:10426')
stopifnot(
  yelloweye$scientific_name == 'Sebastes ruberrimus',
  yelloweye$assessment_area == 'Pacific Coast',
  yelloweye$inventory_assessment_year == 2025L,
  yelloweye$match_status == 'new', is.na(yelloweye$matched_stock_id),
  yelloweye$screening_status == 'screened',
  yelloweye$eligibility == 'eligible',
  yelloweye$import_status == 'not_started',
  grepl('seven fishery fleets', yelloweye$catch_age_evidence, fixed = TRUE),
  grepl('conditional age-at-length', yelloweye$catch_age_evidence, fixed = TRUE),
  grepl('direct age-specific', yelloweye$survey_age_evidence, fixed = TRUE),
  grepl('growth curve', yelloweye$weights_evidence, fixed = TRUE),
  grepl('0.044 per year', yelloweye$mortality_evidence, fixed = TRUE),
  grepl('January', yelloweye$timing_evidence, fixed = TRUE),
  grepl('477.63 billion eggs', yelloweye$output_evidence, fixed = TRUE),
  grepl('two areas', yelloweye$notes, fixed = TRUE),
  grepl('3.30.23.2', yelloweye$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_west_coast_yelloweye_2025.md'))
)
cat('NOAA West Coast Yelloweye identity and source-screening checks passed.\n')

yellowtail <- lookup('NOAA:10393')
stopifnot(
  yellowtail$scientific_name == 'Sebastes flavidus',
  yellowtail$stock_name == 'Yellowtail rockfish - Northern Pacific Coast',
  yellowtail$inventory_assessment_year == 2025L,
  yellowtail$match_status == 'new', is.na(yellowtail$matched_stock_id),
  yellowtail$screening_status == 'screened',
  yellowtail$eligibility == 'eligible',
  yellowtail$import_status == 'not_started',
  grepl('three fisheries', yellowtail$catch_age_evidence, fixed = TRUE),
  grepl('age and length compositions', yellowtail$catch_age_evidence, fixed = TRUE),
  grepl('four survey indices', yellowtail$survey_age_evidence, fixed = TRUE),
  grepl('232 histological samples', yellowtail$maturity_evidence, fixed = TRUE),
  grepl('age 43', yellowtail$mortality_evidence, fixed = TRUE),
  grepl('July', yellowtail$timing_evidence, fixed = TRUE),
  grepl('3.30.23.1', yellowtail$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_yellowtail_north_2025.md'))
)
cat('NOAA northern Yellowtail Rockfish identity and source-screening checks passed.\n')

rougheye <- lookup('NOAA:14516')
stopifnot(
  identical(rougheye$scientific_name, 'Sebastes aleutianus; Sebastes melanostictus'),
  rougheye$assessment_area == 'U.S. West Coast',
  rougheye$inventory_assessment_year == 2025L,
  rougheye$match_status == 'new', is.na(rougheye$matched_stock_id),
  rougheye$screening_status == 'screened',
  rougheye$eligibility == 'eligible',
  rougheye$import_status == 'not_started',
  grepl('six fisheries', rougheye$catch_age_evidence, fixed = TRUE),
  grepl('conditional age-at-length', rougheye$catch_age_evidence, fixed = TRUE),
  grepl('Triennial, Alaska slope', rougheye$survey_age_evidence, fixed = TRUE),
  grepl('Male natural mortality is fixed at 0.036', rougheye$mortality_evidence, fixed = TRUE),
  grepl('3.30.23.2', rougheye$notes, fixed = TRUE),
  grepl('not be counted as two cases', rougheye$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_rougheye_blackspotted_2025.md'))
)
cat('NOAA Rougheye/Blackspotted identity and source-screening checks passed.\n')

quillback <- lookup('NOAA:15342')
stopifnot(
  quillback$scientific_name == 'Sebastes maliger',
  quillback$stock_name == 'Quillback rockfish - California',
  quillback$assessment_area == 'California',
  quillback$inventory_assessment_year == 2025L,
  quillback$match_status == 'new', is.na(quillback$matched_stock_id),
  quillback$screening_status == 'screened',
  quillback$eligibility == 'eligible',
  quillback$import_status == 'not_started',
  grepl('total removals in metric', quillback$catch_age_evidence, fixed = TRUE),
  grepl('2007-2024', quillback$catch_age_evidence, fixed = TRUE),
  grepl('CCFRP', quillback$survey_age_evidence, fixed = TRUE),
  grepl('0.068 per year', quillback$mortality_evidence, fixed = TRUE),
  grepl('July', quillback$timing_evidence, fixed = TRUE),
  grepl('3.30.23.1', quillback$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_quillback_california_2025.md'))
)
cat('NOAA California Quillback identity and source-screening checks passed.\n')

ccgom_yellowtail <- lookup('NOAA:10035')
noaa_progress <- progress[progress$organization == 'NOAA', , drop = FALSE]
stopifnot(
  nrow(ccgom_yellowtail) == 1L,
  ccgom_yellowtail$match_status == 'new',
  is.na(ccgom_yellowtail$matched_stock_id),
  ccgom_yellowtail$screening_status == 'screened',
  ccgom_yellowtail$eligibility == 'unresolved',
  ccgom_yellowtail$import_status == 'not_started',
  grepl('catch-at-age diagnostics', ccgom_yellowtail$catch_age_evidence, fixed = TRUE),
  grepl('aggregate survey indices', ccgom_yellowtail$survey_age_evidence, fixed = TRUE),
  grepl('January 1 numbers-at-age', ccgom_yellowtail$output_evidence, fixed = TRUE),
  grepl('best available science', ccgom_yellowtail$decision_reason, fixed = TRUE),
  grepl('numeric age inputs', ccgom_yellowtail$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'noaa_yellowtail_ccgom_2025.md')),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$ambiguous) == 3L,
  sum(noaa_progress$identity_unreviewed) == 351L,
  sum(noaa_progress$eligibility_unresolved) == 332L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L
)
cat('NOAA Cape Cod/Gulf of Maine yellowtail identity and source-screening checks passed.\n')
icelandic_cod <- lookup('ICES:cod.27.5a')
stopifnot(
  nrow(icelandic_cod) == 1L,
  icelandic_cod$match_status == 'matched',
  icelandic_cod$matched_stock_id == 'ck_ices_nwwg_morhua_ice5a',
  icelandic_cod$screening_status == 'screened',
  icelandic_cod$eligibility == 'eligible',
  icelandic_cod$import_status == 'imported',
  grepl('1955-2025', icelandic_cod$catch_age_evidence, fixed = TRUE),
  grepl('IS-SMB', icelandic_cod$survey_age_evidence, fixed = TRUE),
  grepl('1955-2026', icelandic_cod$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_cod_iceland_2026.md'))
)
cat('ICES cod 5.a identity and current-run screen passed.\n')
northsea_ices_ids <- c('ICES:her.27.3a47d', 'ICES:had.27.46a20',
  'ICES:ple.27.420', 'ICES:pok.27.3a46', 'ICES:whg.27.47d')
northsea_ices <- register[match(northsea_ices_ids, register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(northsea_ices$candidate_id),
  identical(northsea_ices$matched_stock_id, c('ices_herring_north_sea',
    'ices_haddock_north_sea', 'ices_plaice_north_sea', 'ices_saithe_north_sea',
    'ices_whiting_north_sea')),
  all(northsea_ices$match_status == 'matched'),
  all(northsea_ices$screening_status == 'screened'),
  all(northsea_ices$eligibility == 'eligible'),
  all(northsea_ices$import_status == 'imported'),
  grepl('winter-ring', northsea_ices$catch_age_evidence[1], fixed = TRUE),
  grepl('partial spawning-component', northsea_ices$survey_age_evidence[1], fixed = TRUE),
  grepl('Q1', northsea_ices$survey_age_evidence[2], fixed = TRUE),
  grepl('Five age-specific', northsea_ices$survey_age_evidence[3], fixed = TRUE),
  grepl('Q3 and Q4', northsea_ices$survey_age_evidence[4], fixed = TRUE),
  grepl('April NSwhiting_2026n', northsea_ices$catch_age_evidence[5], fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    c('ices_herring_north_sea.md', 'ices_haddock_north_sea.md',
      'ices_plaice_north_sea.md', 'ices_saithe_north_sea.md',
      'ices_whiting_north_sea_2026.md'))))
)
cat('Matched ICES North Sea identity and source-screening checks passed.\n')
dfo_matched_ids <- c(
  'DFO:Newfoundland and Labrador:Atlantic Cod - NAFO 3Ps',
  'DFO:Gulf:Atlantic Cod - NAFO 4TVn',
  'DFO:Maritimes:Haddock - NAFO 4X5Y',
  'DFO:Pacific:Pacific Hake - Offshore')
dfo_matched <- register[match(dfo_matched_ids, register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(dfo_matched$candidate_id),
  identical(dfo_matched$matched_stock_id, c('ck_dfo_morhua_3ps',
    'dfo_cod_4t4vn', 'ck_dfo_aeglefinus_4x5y', 'joint_hake_coastal')),
  all(dfo_matched$match_status == 'matched'),
  all(dfo_matched$screening_status == 'screened'),
  all(dfo_matched$eligibility == 'eligible'),
  all(dfo_matched$import_status == 'imported'),
  grepl('2025 HYBRID update', dfo_matched$decision_reason[1], fixed = TRUE),
  grepl('separately', dfo_matched$decision_reason[2], fixed = TRUE),
  grepl('summary-only', dfo_matched$decision_reason[3], fixed = TRUE),
  grepl('composition-based', dfo_matched$decision_reason[4], fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    c('dfo_cod_3ps_2024.md', 'dfo_cod_4t4vn.md',
      'dfo_haddock_4x5y.md', 'joint_hake_coastal.md'))))
)
cat('Matched DFO identity and source-screening checks passed.\n')
nwwg_match_ids <- c('ICES:pok.27.5a', 'ICES:had.27.5a', 'ICES:cod.27.5b1',
  'ICES:had.27.5b', 'ICES:pok.27.5b')
nwwg_matches <- register[match(nwwg_match_ids, register$candidate_id), , drop = FALSE]
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  !anyNA(nwwg_matches$candidate_id),
  identical(nwwg_matches$matched_stock_id, c('ck_ices_nwwg_virens_fa5a',
    'ices_haddock_iceland', 'ck_ices_nwwg_morhua_fapl5b1',
    'ck_ices_nwwg_aeglefinus_fa5b', 'ck_ices_nwwg_virens_fa5b')),
  all(nwwg_matches$match_status == 'matched'),
  all(nwwg_matches$screening_status == 'screened'),
  all(nwwg_matches$eligibility == 'eligible'),
  all(nwwg_matches$import_status == 'imported'),
  grepl('3-14+', nwwg_matches$catch_age_evidence[1], fixed = TRUE),
  grepl('IS-SMB', nwwg_matches$survey_age_evidence[2], fixed = TRUE),
  grepl('1959-2025', nwwg_matches$catch_age_evidence[3], fixed = TRUE),
  grepl('official 2025 SAM production run', nwwg_matches$decision_reason[4], fixed = TRUE),
  grepl('1961-2025', nwwg_matches$catch_age_evidence[5], fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    c('ices_saithe_iceland_2026.md', 'ices_haddock_iceland_2026.md',
      'ices_cod_faroe_plateau_2025.md', 'ices_haddock_faroe_2025.md',
      'ices_saithe_faroe_2025.md')))),
  sum(ices_progress$inventory_rows) == 274L,
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$new) == 18L,
  sum(ices_progress$ambiguous) == 2L,
  sum(ices_progress$identity_unreviewed) == 197L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched Iceland/Faroe ICES identity and source-screening checks passed.\n')
wgcse_distinct <- register[match(c('ICES:ple.27.7h-k', 'ICES:had.27.7a',
  'ICES:cod.27.6b', 'ICES:ple.27.7bc', 'ICES:sol.27.7bc'), register$candidate_id), , drop = FALSE]
wgcse_ambiguous <- register[match(c('ICES:lez.27.4a6a', 'ICES:lez.27.6b'),
  register$candidate_id), , drop = FALSE]
stopifnot(
  !anyNA(wgcse_distinct$candidate_id),
  all(wgcse_distinct$match_status == 'new'),
  all(is.na(wgcse_distinct$matched_stock_id)),
  all(wgcse_distinct$screening_status == 'screened'),
  all(wgcse_distinct$eligibility == 'ineligible'),
  all(wgcse_distinct$import_status == 'not_started'),
  all(grepl('cross-organization reconciliation', wgcse_distinct$notes, fixed = TRUE)),
  !anyNA(wgcse_ambiguous$candidate_id),
  all(wgcse_ambiguous$match_status == 'ambiguous'),
  all(is.na(wgcse_ambiguous$matched_stock_id)),
  all(wgcse_ambiguous$screening_status == 'screened'),
  all(wgcse_ambiguous$eligibility == 'ineligible'),
  all(grepl('Lepidorhombus spp.', wgcse_ambiguous$notes, fixed = TRUE)),
  sum(ices_progress$new) == 18L,
  sum(ices_progress$ambiguous) == 2L,
  sum(ices_progress$identity_unreviewed) == 197L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L
)
cat('ICES WGCSE identity reconciliation checks passed.\n')
dfo_plaice <- lookup('DFO:Gulf:American Plaice - NAFO 4T')
dfo_progress <- progress[progress$organization == 'DFO', , drop = FALSE]
stopifnot(
  nrow(dfo_plaice) == 1L,
  dfo_plaice$matched_stock_id == 'ck_dfo_platessoides_4t',
  dfo_plaice$match_status == 'matched',
  dfo_plaice$screening_status == 'screened',
  dfo_plaice$eligibility == 'eligible',
  dfo_plaice$import_status == 'imported',
  grepl('1976-2015', dfo_plaice$catch_age_evidence, fixed = TRUE),
  grepl('uncalibrated vessel', dfo_plaice$survey_age_evidence, fixed = TRUE),
  grepl('0.075', dfo_plaice$mortality_evidence, fixed = TRUE),
  grepl('2025 indicator update', dfo_plaice$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ck_dfo_platessoides_4t.md')),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L
)
cat('Matched DFO Gulf American plaice source-screening checks passed.\n')
channel_plaice <- lookup('ICES:ple.27.7d')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(channel_plaice) == 1L,
  channel_plaice$matched_stock_id == 'ck_ices_wgnssk_platessa_ns_7d',
  channel_plaice$match_status == 'matched',
  channel_plaice$screening_status == 'screened',
  channel_plaice$eligibility == 'eligible',
  channel_plaice$import_status == 'imported',
  grepl('1980-2025', channel_plaice$catch_age_evidence, fixed = TRUE),
  grepl('UK BTS', channel_plaice$survey_age_evidence, fixed = TRUE),
  grepl('0.875', channel_plaice$timing_evidence, fixed = TRUE),
  grepl('N-at-age, F-at-age', channel_plaice$output_evidence, fixed = TRUE),
  grepl('native fitted-object endpoint was unavailable', channel_plaice$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_plaice_english_channel_7d_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Division 7.d plaice source-screening checks passed.\n')
northeast_arctic_haddock <- lookup('ICES:had.27.1-2')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(northeast_arctic_haddock) == 1L,
  northeast_arctic_haddock$matched_stock_id == 'ck_ices_afwg_aeglefinus_nea_1_2',
  northeast_arctic_haddock$match_status == 'matched',
  northeast_arctic_haddock$screening_status == 'screened',
  northeast_arctic_haddock$eligibility == 'eligible',
  northeast_arctic_haddock$import_status == 'imported',
  grepl('836 residual catch-at-age', northeast_arctic_haddock$catch_age_evidence, fixed = TRUE),
  grepl('790 non-missing observations', northeast_arctic_haddock$survey_age_evidence, fixed = TRUE),
  grepl('0.133', northeast_arctic_haddock$timing_evidence, fixed = TRUE),
  grepl('fixed in the final fit', northeast_arctic_haddock$mortality_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'jrn_afwg_haddock_northeast_arctic_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Northeast Arctic haddock source-screening checks passed.\n')
norpout <- lookup('ICES:nop.27.3a4')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(norpout) == 1L,
  norpout$matched_stock_id == 'ices_norway_pout_north_sea',
  norpout$match_status == 'matched',
  norpout$screening_status == 'screened',
  norpout$eligibility == 'eligible',
  norpout$import_status == 'imported',
  grepl('1984-2025 Q3', norpout$catch_age_evidence, fixed = TRUE),
  grepl('424 survey rows', norpout$survey_age_evidence, fixed = TRUE),
  grepl('not estimated by SESAM', norpout$mortality_evidence, fixed = TRUE),
  grepl('not an estimated output', norpout$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_norway_pout_north_sea.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES North Sea Norway pout source-screening checks passed.\n')
nssh <- lookup('ICES:her.27.1-24a514a')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(nssh) == 1L,
  nssh$matched_stock_id == 'ices_herring_norwegian_spring_spawning',
  nssh$match_status == 'matched',
  nssh$screening_status == 'screened',
  nssh$eligibility == 'eligible',
  nssh$import_status == 'imported',
  grepl('1950-2025', nssh$catch_age_evidence, fixed = TRUE),
  grepl('NASF acoustic', nssh$survey_age_evidence, fixed = TRUE),
  grepl('M deviations', nssh$mortality_evidence, fixed = TRUE),
  grepl('through 2026', nssh$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_herring_norwegian_spring_spawning.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Norwegian spring-spawning herring source-screening checks passed.\n')
golden_redfish <- lookup('ICES:reg.27.1-2')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(golden_redfish) == 1L,
  golden_redfish$matched_stock_id == 'ck_ices_afwg_norvegicus_deep1_2',
  golden_redfish$match_status == 'matched',
  golden_redfish$screening_status == 'screened',
  golden_redfish$eligibility == 'ineligible',
  golden_redfish$import_status == 'imported',
  grepl('total commercial landings', golden_redfish$catch_age_evidence, fixed = TRUE),
  grepl('aggregate bottom-trawl biomass index', golden_redfish$survey_age_evidence, fixed = TRUE),
  grepl('replaced the historical age-length Gadget model', golden_redfish$decision_reason, fixed = TRUE),
  grepl('B/BMSY', golden_redfish$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_golden_redfish_northeast_arctic.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Northeast Arctic golden redfish screening checks passed.\n')
gulf_riga_herring <- lookup('ICES:her.27.28')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(gulf_riga_herring) == 1L,
  gulf_riga_herring$matched_stock_id == 'ck_ices_wgbfas_harengus_bs_28_1',
  gulf_riga_herring$match_status == 'matched',
  gulf_riga_herring$screening_status == 'screened',
  gulf_riga_herring$eligibility == 'eligible',
  gulf_riga_herring$import_status == 'imported',
  grepl('1977-2025', gulf_riga_herring$catch_age_evidence, fixed = TRUE),
  grepl('GRAHS acoustic', gulf_riga_herring$survey_age_evidence, fixed = TRUE),
  grepl('0.675', gulf_riga_herring$timing_evidence, fixed = TRUE),
  grepl('conditional N/F state uncertainty', gulf_riga_herring$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_herring_gulf_riga_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Gulf of Riga herring source-screening checks passed.\n')
sandeel_ids <- c('ICES:san.sa.1r', 'ICES:san.sa.2r', 'ICES:san.sa.3r')
sandeels <- register[match(sandeel_ids, register$candidate_id), , drop = FALSE]
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  !anyNA(sandeels$candidate_id),
  identical(sandeels$matched_stock_id, c('ices_sandeel_area_1r',
    'ices_sandeel_area_2r', 'ices_sandeel_area_3r')),
  all(sandeels$match_status == 'matched'),
  all(sandeels$screening_status == 'screened'),
  all(sandeels$eligibility == 'eligible'),
  all(sandeels$import_status == 'imported'),
  grepl('1983-2025', sandeels$catch_age_evidence[1], fixed = TRUE),
  grepl('RTM indices', sandeels$survey_age_evidence[1], fixed = TRUE),
  grepl('long-term multispecies averages', sandeels$mortality_evidence[2], fixed = TRUE),
  grepl('2025 source zeros are retained', sandeels$catch_age_evidence[3], fixed = TRUE),
  grepl('acoustic ages 1-4+', sandeels$survey_age_evidence[3], fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    c('ices_sandeel_area_1r.md', 'ices_sandeel_area_2r.md', 'ices_sandeel_area_3r.md')))),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Sandeel Areas 1r-3r source-screening checks passed.\n')
area4_sandeel <- lookup('ICES:san.sa.4')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(area4_sandeel) == 1L,
  area4_sandeel$matched_stock_id == 'ices_sandeel_area_4',
  area4_sandeel$match_status == 'matched',
  area4_sandeel$screening_status == 'screened',
  area4_sandeel$eligibility == 'eligible',
  area4_sandeel$import_status == 'imported',
  grepl('1993-2024', area4_sandeel$catch_age_evidence, fixed = TRUE),
  grepl('Dredge and Dredge2', area4_sandeel$survey_age_evidence, fixed = TRUE),
  grepl('separate 2026 category-5 advice update', area4_sandeel$output_evidence, fixed = TRUE),
  grepl('2026 summary event separate', area4_sandeel$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_sandeel_area_4.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Sandeel Area 4 source-screening checks passed.\n')
north_sea_sole <- lookup('ICES:sol.27.4')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(north_sea_sole) == 1L,
  north_sea_sole$matched_stock_id == 'ck_ices_wgnssk_solea_ns_4',
  north_sea_sole$match_status == 'matched',
  north_sea_sole$screening_status == 'screened',
  north_sea_sole$eligibility == 'unresolved',
  north_sea_sole$import_status == 'imported',
  grepl('1957-2025', north_sea_sole$catch_age_evidence, fixed = TRUE),
  grepl('aggregate BTS', north_sea_sole$survey_age_evidence, fixed = TRUE),
  grepl('not recovered', north_sea_sole$survey_age_evidence, fixed = TRUE),
  grepl('prior 2025 age-index inputs are not substituted', north_sea_sole$survey_age_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_sole_north_sea.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES North Sea sole source-screening checks passed.\n')
baltic_sole <- lookup('ICES:sol.27.20-24')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(baltic_sole) == 1L,
  baltic_sole$matched_stock_id == 'ck_ices_wgbfas_solea_ns_iiia_22_24',
  baltic_sole$match_status == 'matched',
  baltic_sole$screening_status == 'screened',
  baltic_sole$eligibility == 'eligible',
  baltic_sole$import_status == 'imported',
  grepl('1984-2025', baltic_sole$catch_age_evidence, fixed = TRUE),
  grepl('Age-specific Q4 survey indices', baltic_sole$survey_age_evidence, fixed = TRUE),
  grepl('exact timing and observation precision are not reported', baltic_sole$survey_age_evidence, fixed = TRUE),
  grepl('fixed at 0.1', baltic_sole$mortality_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_sole_20_24_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Baltic sole source-screening checks passed.\n')
blue_whiting <- lookup('ICES:whb.27.1-91214')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(blue_whiting) == 1L,
  blue_whiting$matched_stock_id == 'ices_bluewhiting_northeast_atlantic',
  blue_whiting$match_status == 'matched',
  blue_whiting$screening_status == 'screened',
  blue_whiting$eligibility == 'eligible',
  blue_whiting$import_status == 'imported',
  grepl('1981-2026', blue_whiting$catch_age_evidence, fixed = TRUE),
  grepl('ages 1-8', blue_whiting$survey_age_evidence, fixed = TRUE),
  grepl('0.245', blue_whiting$timing_evidence, fixed = TRUE),
  grepl('fixed at 0.2', blue_whiting$mortality_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_bluewhiting_northeast_atlantic_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES blue whiting source-screening checks passed.\n')
icelandic_summer_herring <- lookup('ICES:her.27.5a')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(icelandic_summer_herring) == 1L,
  icelandic_summer_herring$matched_stock_id == 'ck_ices_nwwg_harengus_ice_summer',
  icelandic_summer_herring$match_status == 'matched',
  icelandic_summer_herring$screening_status == 'screened',
  icelandic_summer_herring$eligibility == 'unresolved',
  icelandic_summer_herring$import_status == 'imported',
  grepl('1980-2023', icelandic_summer_herring$catch_age_evidence, fixed = TRUE),
  grepl('mapping to the production run is unresolved', icelandic_summer_herring$survey_age_evidence, fixed = TRUE),
  grepl('infection-related', icelandic_summer_herring$mortality_evidence, fixed = TRUE),
  grepl('not splice later raw catch/acoustic', icelandic_summer_herring$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ck_ices_nwwg_harengus_ice_summer.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched Icelandic summer-spawning herring source-screening checks passed.\n')
baltic_sprat <- lookup('ICES:spr.27.22-32')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(baltic_sprat) == 1L,
  baltic_sprat$matched_stock_id == 'ices_sprat_baltic',
  baltic_sprat$match_status == 'matched',
  baltic_sprat$screening_status == 'screened',
  baltic_sprat$eligibility == 'eligible',
  baltic_sprat$import_status == 'imported',
  grepl('1974-2025', baltic_sprat$catch_age_evidence, fixed = TRUE),
  grepl('Four age-specific tuning fleets', baltic_sprat$survey_age_evidence, fixed = TRUE),
  grepl('cod predation', baltic_sprat$mortality_evidence, fixed = TRUE),
  grepl('confidence levels are unspecified', baltic_sprat$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_sprat_baltic.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Baltic sprat source-screening checks passed.\n')
southern_horse_mackerel <- lookup('ICES:hom.27.9a')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(southern_horse_mackerel) == 1L,
  southern_horse_mackerel$matched_stock_id == 'ck_ices_wghansa_trachurus_sp9a',
  southern_horse_mackerel$match_status == 'matched',
  southern_horse_mackerel$screening_status == 'screened',
  southern_horse_mackerel$eligibility == 'unresolved',
  southern_horse_mackerel$import_status == 'imported',
  grepl('1992-2023', southern_horse_mackerel$catch_age_evidence, fixed = TRUE),
  grepl('commercial CPUE and aggregate IBTS biomass series were not recovered',
    southern_horse_mackerel$survey_age_evidence, fixed = TRUE),
  grepl('exact sampling time', southern_horse_mackerel$timing_evidence, fixed = TRUE),
  grepl('no numerical age-specific N/F outputs', southern_horse_mackerel$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_horse_mackerel_9a.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Southern horse mackerel source-screening checks passed.\n')
northeast_arctic_cod <- lookup('ICES:cod.27.1-2')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(northeast_arctic_cod) == 1L,
  northeast_arctic_cod$matched_stock_id == 'ices_cod_northeast_arctic',
  northeast_arctic_cod$match_status == 'matched',
  northeast_arctic_cod$screening_status == 'screened',
  northeast_arctic_cod$eligibility == 'eligible',
  northeast_arctic_cod$import_status == 'imported',
  grepl('1946-2025', northeast_arctic_cod$catch_age_evidence, fixed = TRUE),
  grepl('five age-specific index series', northeast_arctic_cod$survey_age_evidence, fixed = TRUE),
  grepl('externally iterated cannibalism', northeast_arctic_cod$mortality_evidence, fixed = TRUE),
  grepl('2026 JRN-AFWG', northeast_arctic_cod$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_cod_northeast_arctic.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES Northeast Arctic cod source-screening checks passed.\n')
sole_7d <- lookup('ICES:sol.27.7d')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(sole_7d) == 1L,
  sole_7d$matched_stock_id == 'ck_ices_wgnssk_solea_ns_7d',
  sole_7d$match_status == 'matched',
  sole_7d$screening_status == 'screened',
  sole_7d$eligibility == 'eligible',
  sole_7d$import_status == 'imported',
  grepl('1982-2025', sole_7d$catch_age_evidence, fixed = TRUE),
  grepl('UK BTS indices', sole_7d$survey_age_evidence, fixed = TRUE),
  grepl('fixed at 0.1', sole_7d$mortality_evidence, fixed = TRUE),
  grepl('0.02%', sole_7d$output_evidence, fixed = TRUE),
  grepl('non-converged tinyAM fit', sole_7d$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_sole_eastern_channel_7d_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched ICES 7.d sole source-screening checks passed.\n')
western_baltic_herring <- lookup('ICES:her.27.20-24')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(western_baltic_herring) == 1L,
  western_baltic_herring$matched_stock_id == 'ices_herring_western_baltic',
  western_baltic_herring$match_status == 'matched',
  western_baltic_herring$screening_status == 'screened',
  western_baltic_herring$eligibility == 'eligible',
  western_baltic_herring$import_status == 'imported',
  grepl('1991-2025', western_baltic_herring$catch_age_evidence, fixed = TRUE),
  grepl('372 finite survey observations', western_baltic_herring$survey_age_evidence, fixed = TRUE),
  grepl('0.168 and 0.25', western_baltic_herring$timing_evidence, fixed = TRUE),
  grepl('0-5+', western_baltic_herring$notes, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_herring_western_baltic_2026.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched western Baltic herring source-screening checks passed.\n')
beaked_redfish <- lookup('ICES:reb.27.1-2')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(beaked_redfish) == 1L,
  beaked_redfish$matched_stock_id == 'ck_ices_afwg_mentella_deep_1_2',
  beaked_redfish$match_status == 'matched',
  beaked_redfish$screening_status == 'screened',
  beaked_redfish$eligibility == 'unresolved',
  beaked_redfish$import_status == 'imported',
  grepl('1992-2025', beaked_redfish$catch_age_evidence, fixed = TRUE),
  grepl('year class', beaked_redfish$survey_age_evidence, fixed = TRUE),
  grepl('0.05 per year', beaked_redfish$mortality_evidence, fixed = TRUE),
  grepl('JRN-AFWG 2026', beaked_redfish$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ck_ices_afwg_mentella_deep_1_2.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Matched Northeast Arctic beaked redfish source-screening checks passed.\n')

celtic_herring <- lookup('ICES:her.27.irls')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(celtic_herring) == 1L,
  is.na(celtic_herring$matched_stock_id),
  celtic_herring$match_status == 'new',
  celtic_herring$screening_status == 'screened',
  celtic_herring$eligibility == 'eligible',
  celtic_herring$import_status == 'not_started',
  grepl('1958-2025', celtic_herring$catch_age_evidence, fixed = TRUE),
  grepl('ages 2-7', celtic_herring$survey_age_evidence, fixed = TRUE),
  grepl('October', celtic_herring$timing_evidence, fixed = TRUE),
  grepl('unquantified mixing', celtic_herring$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'ices_celtic_herring_2026.md')),
  sum(ices_progress$new) == 18L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('ICES Celtic Sea herring source-screening checks passed.\n')

irish_herring <- lookup('ICES:her.27.nirs')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(irish_herring) == 1L,
  is.na(irish_herring$matched_stock_id),
  irish_herring$match_status == 'new',
  irish_herring$screening_status == 'screened',
  irish_herring$eligibility == 'eligible',
  irish_herring$import_status == 'not_started',
  grepl('1980-2025', irish_herring$catch_age_evidence, fixed = TRUE),
  grepl('ages 1-8+', irish_herring$survey_age_evidence, fixed = TRUE),
  grepl('borrowed from North Sea', irish_herring$mortality_evidence, fixed = TRUE),
  grepl('2.7 to 4.9', irish_herring$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'ices_irish_sea_herring_2026.md')),
  sum(ices_progress$new) == 18L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('ICES Irish Sea herring source-screening checks passed.\n')

northern_shelf_cod <- lookup('ICES:cod.27.46a7d20')
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(northern_shelf_cod) == 1L,
  northern_shelf_cod$matched_stock_id == 'ices_cod_north_sea',
  northern_shelf_cod$match_status == 'matched',
  northern_shelf_cod$screening_status == 'screened',
  northern_shelf_cod$eligibility == 'unresolved',
  northern_shelf_cod$import_status == 'not_started',
  grepl('cn_WS.dat', northern_shelf_cod$catch_age_evidence, fixed = TRUE),
  grepl('codData210926.Rdata', northern_shelf_cod$survey_age_evidence, fixed = TRUE),
  grepl('26001-26003', northern_shelf_cod$output_evidence, fixed = TRUE),
  grepl('2026-10-09', northern_shelf_cod$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'ices_cod_north_sea.md')),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('ICES Northern Shelf cod 2026 source screening checks passed.\n')

reviewed_ices <- c(
  'ICES:ane.27.9aS', 'ICES:ane.27.9aW', 'ICES:cod.27.1-2.coastN',
  'ICES:cod.27.2.coastS', 'ICES:cod.27.22-24', 'ICES:her.27.25-2932',
  'ICES:her.27.3031', 'ICES:her.27.6aN', 'ICES:her.27.6aS7bc',
  'ICES:pil.27.8abd', 'ICES:ple.27.21-32', 'ICES:pok.27.1-2'
)
reviewed <- register[match(reviewed_ices, register$candidate_id), , drop = FALSE]
ices_progress <- progress[progress$organization == 'ICES', , drop = FALSE]
stopifnot(
  nrow(reviewed) == 12L,
  all(reviewed$screening_status == 'screened'),
  all(reviewed$import_status == 'imported'),
  identical(reviewed$eligibility, c(
    'unresolved', 'ineligible', 'eligible', 'ineligible', 'eligible',
    'eligible', 'eligible', 'ineligible', 'ineligible', 'unresolved',
    'eligible', 'eligible'
  )),
  grepl('not recovered', reviewed$survey_age_evidence[1], fixed = TRUE),
  grepl('category-3.2', reviewed$decision_reason[2], fixed = TRUE),
  grepl('1994-2025', reviewed$output_evidence[3], fixed = TRUE),
  grepl('SPiCT', reviewed$decision_reason[4], fixed = TRUE),
  grepl('central 90%', reviewed$output_evidence[5], fixed = TRUE),
  grepl('central 90%', reviewed$output_evidence[6], fixed = TRUE),
  grepl('0.75 as an explicit approximation', reviewed$timing_evidence[7], fixed = TRUE),
  grepl('survey-biomass trend', reviewed$decision_reason[8], fixed = TRUE),
  grepl('exploratory age-based', reviewed$decision_reason[9], fixed = TRUE),
  grepl('omits the production numerical input files', reviewed$decision_reason[10], fixed = TRUE),
  grepl('0.1507111', reviewed$timing_evidence[11], fixed = TRUE),
  grepl('sampling time of 0.8', reviewed$timing_evidence[12], fixed = TRUE),
  all(file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    c('ices_anchovy_gulf_cadiz.md', 'ices_anchovy_west_iberia_2026.md',
      'ices_cod_northern_coastal_2026.md', 'ices_cod_southern_coastal_2026.md',
      'ck_ices_wgbfas_morhua_bs.md',
      'ck_ices_wgbfas_harengus_bs_25_29_32.md',
      'ices_herring_gulf_bothnia_2026.md',
      'ck_ices_hawg_harengus_cs_6a_n.md',
      'ck_ices_hawg_harengus_cs6_a_s_7b_c.md',
      'ck_ices_wghansa_pilchardus_sp8abd.md',
      'ices_plaice_baltic_21_32_2026.md',
      'ck_ices_afwg_virens_nea_1_2.md')))),
  sum(ices_progress$matched) == 57L,
  sum(ices_progress$eligible) == 45L,
  sum(ices_progress$excluded) == 51L,
  sum(ices_progress$eligibility_unresolved) == 178L,
  sum(ices_progress$screening_complete) == 109L,
  sum(ices_progress$screening_not_started) == 165L
)
cat('Twelve previously matched ICES candidates passed source-screening checks.\n')

pacific_herring_ids <- c(
  'DFO:Pacific:Pacific Herring - Central Coast',
  'DFO:Pacific:Pacific Herring - Haida Gwaii',
  'DFO:Pacific:Pacific Herring - Prince Rupert District',
  'DFO:Pacific:Pacific Herring - West Coast of Vancouver Island'
)
pacific_herring <- register[match(pacific_herring_ids, register$candidate_id), , drop = FALSE]
stopifnot(
  nrow(pacific_herring) == 4L,
  identical(pacific_herring$official_stock_id, c('CC', 'HG', 'PRD', 'WCVI')),
  identical(pacific_herring$scientific_name, rep('Clupea pallasii', 4L)),
  identical(pacific_herring$match_status, rep('new', 4L)),
  identical(pacific_herring$screening_status, rep('screened', 4L)),
  identical(pacific_herring$eligibility, rep('ineligible', 4L)),
  identical(pacific_herring$import_status, rep('not_started', 4L)),
  all(grepl('aggregate', pacific_herring$survey_age_evidence, fixed = TRUE)),
  all(grepl('2026-014|2026-006', pacific_herring$detailed_source_url)),
  all(grepl('2026-10-29', pacific_herring$review_date, fixed = TRUE)),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_pacific_herring_areas_discovery.md')),
  sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L
)
cat('DFO Pacific herring regional identity and source-screening checks passed.\n')

gulf_yellowtail <- lookup('DFO:Gulf:Yellowtail Flounder - NAFO 4T')
stopifnot(
  nrow(gulf_yellowtail) == 1L,
  is.na(gulf_yellowtail$official_stock_id),
  gulf_yellowtail$scientific_name == 'Myzopsetta ferruginea',
  gulf_yellowtail$inventory_assessment_year == 2025L,
  gulf_yellowtail$assessment_area == 'NAFO Division 4T, southern Gulf of St. Lawrence',
  gulf_yellowtail$match_status == 'new',
  is.na(gulf_yellowtail$matched_stock_id),
  gulf_yellowtail$screening_status == 'screened',
  gulf_yellowtail$eligibility == 'eligible',
  gulf_yellowtail$import_status == 'not_started',
  grepl('1985-2022', gulf_yellowtail$catch_age_evidence, fixed = TRUE),
  grepl('1985-2025', gulf_yellowtail$survey_age_evidence, fixed = TRUE),
  grepl('age-aggregated', gulf_yellowtail$survey_age_evidence, fixed = TRUE),
  grepl('4-7, 8-9 and 10+', gulf_yellowtail$mortality_evidence, fixed = TRUE),
  grepl('8,269 t', gulf_yellowtail$output_evidence, fixed = TRUE),
  grepl('2026-054', gulf_yellowtail$detailed_source_url, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_yellowtail_flounder_4t_2025.md')),
  sum(progress$new) == 89L,
  sum(progress$eligible) == 99L,
  sum(progress$eligible_new) == 30L,
  sum(progress$screening_complete) == 369L
)
cat('DFO NAFO 4T yellowtail identity and source-screening checks passed.\n')

gulf_white_hake <- lookup('DFO:Gulf:White Hake - NAFO 4T')
stopifnot(
  nrow(gulf_white_hake) == 1L,
  is.na(gulf_white_hake$official_stock_id),
  gulf_white_hake$scientific_name == 'Urophycis tenuis',
  gulf_white_hake$inventory_assessment_year == 2022L,
  gulf_white_hake$assessment_area == 'NAFO Division 4T, southern Gulf of St. Lawrence',
  gulf_white_hake$match_status == 'new',
  is.na(gulf_white_hake$matched_stock_id),
  gulf_white_hake$screening_status == 'screened',
  gulf_white_hake$eligibility == 'eligible',
  gulf_white_hake$import_status == 'not_started',
  grepl('1978-2022', gulf_white_hake$catch_age_evidence, fixed = TRUE),
  grepl('0-2', gulf_white_hake$catch_age_evidence, fixed = TRUE),
  grepl('fish/tow', gulf_white_hake$survey_age_evidence, fixed = TRUE),
  grepl('September', gulf_white_hake$timing_evidence, fixed = TRUE),
  grepl('6.1 kt', gulf_white_hake$output_evidence, fixed = TRUE),
  grepl('28963929', gulf_white_hake$detailed_source_url, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_white_hake_4t_2022.md')),
  sum(progress$new) == 89L,
  sum(progress$eligible) == 99L,
  sum(progress$eligible_new) == 30L,
  sum(progress$screening_complete) == 369L
)
cat('DFO NAFO 4T white hake identity and source-screening checks passed.\n')
noaa_sardine <- lookup('NOAA:10573')
dfo_sardine <- register[register$organization == 'DFO' &
  grepl('^DFO:Pacific:Sardine', register$candidate_id), , drop = FALSE]
pacific_sardine <- rbind(noaa_sardine, dfo_sardine)
stopifnot(
  nrow(pacific_sardine) == 2L,
  all(pacific_sardine$candidate_unit_id == 'pacific_sardine_northern_subpopulation'),
  all(pacific_sardine$match_status == 'new'),
  all(pacific_sardine$screening_status == 'screened'),
  all(pacific_sardine$eligibility == 'unresolved'),
  all(pacific_sardine$import_status == 'not_started'),
  pacific_sardine$inventory_assessment_year[1L] == 2025L,
  grepl('three commercial fleets', lookup('NOAA:10573')$catch_age_evidence, fixed = TRUE),
  grepl('not published', lookup('NOAA:10573')$catch_age_evidence, fixed = TRUE),
  grepl('age-1+ biomass', lookup('NOAA:10573')$output_evidence, fixed = TRUE),
  grepl('uses NOAA assessment biomass', dfo_sardine$decision_reason, fixed = TRUE),
  !any(stocks$scientific_name == 'Sardinops sagax'),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'pacific_sardine_northern_subpopulation_2025.md'))
)

eastern_georges_haddock <- register[
  !is.na(register$candidate_unit_id) &
    register$candidate_unit_id == 'haddock_eastern_georges_bank', , drop = FALSE
]
stopifnot(
  nrow(eastern_georges_haddock) == 2L,
  identical(sort(eastern_georges_haddock$candidate_id), sort(c(
    'DFO:Maritimes:Haddock - NAFO 5Zjm', 'NOAA:12806'
  ))),
  all(eastern_georges_haddock$match_status == 'new'),
  all(is.na(eastern_georges_haddock$matched_stock_id)),
  all(eastern_georges_haddock$screening_status == 'screened'),
  all(eastern_georges_haddock$eligibility == 'unresolved'),
  lookup('DFO:Maritimes:Haddock - NAFO 5Zjm')$official_stock_id == '5Zjm',
  lookup('DFO:Maritimes:Haddock - NAFO 5Zjm')$inventory_assessment_year == 2024L,
  lookup('NOAA:12806')$inventory_assessment_year == 2023L,
  lookup('NOAA:10520')$matched_stock_id == 'nefsc_haddock_georges_bank',
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'dfo_noaa_eastern_georges_haddock.md')),
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L
)
cat('Eastern Georges Bank haddock inventory crosswalk checks passed.\n')

goa_rock_sole <- register[register$candidate_unit_id %in% c(
  'noaa_goa_northern_rock_sole', 'noaa_goa_southern_rock_sole'
), , drop = FALSE]
goanrs <- goa_rock_sole[goa_rock_sole$candidate_unit_id ==
  'noaa_goa_northern_rock_sole', , drop = FALSE]
goasrs <- goa_rock_sole[goa_rock_sole$candidate_unit_id ==
  'noaa_goa_southern_rock_sole', , drop = FALSE]
stopifnot(
  nrow(goa_rock_sole) == 6L,
  identical(sort(goanrs$candidate_id), sort(c('NOAA:10836', 'NOAA:15412', 'NOAA:15413'))),
  identical(sort(goasrs$candidate_id), sort(c('NOAA:10040', 'NOAA:15414', 'NOAA:15415'))),
  all(goa_rock_sole$match_status == 'new'),
  all(goa_rock_sole$screening_status == 'screened'),
  all(goa_rock_sole$eligibility == 'ineligible'),
  lookup('NOAA:10040')$inventory_assessment_year == 2020L,
  lookup('NOAA:15414')$inventory_assessment_year == 2024L,
  lookup('NOAA:10836')$inventory_assessment_year == 2020L,
  lookup('NOAA:15412')$inventory_assessment_year == 2024L,
  all(grepl('length compositions', goa_rock_sole$catch_age_evidence, fixed = TRUE)),
  all(grepl('not catch-at-age', goa_rock_sole$catch_age_evidence, fixed = TRUE)),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_goa_northern_rock_sole.md')),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_goa_southern_rock_sole.md')),
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L
)
cat('NOAA Gulf of Alaska rock sole identity and age-data screen checks passed.\n')

canary_pacific <- lookup('NOAA:10415')
stopifnot(
  nrow(canary_pacific) == 1L,
  canary_pacific$inventory_assessment_year == 2025L,
  is.na(canary_pacific$candidate_unit_id),
  canary_pacific$match_status == 'new',
  is.na(canary_pacific$matched_stock_id),
  canary_pacific$screening_status == 'screened',
  canary_pacific$eligibility == 'eligible',
  canary_pacific$import_status == 'not_started',
  grepl('catch totals', canary_pacific$catch_age_evidence, fixed = TRUE),
  grepl('age-composition', canary_pacific$survey_age_evidence, fixed = TRUE),
  grepl('month 7 (July)', canary_pacific$timing_evidence, fixed = TRUE),
  grepl('2023 full assessment archive', canary_pacific$decision_reason, fixed = TRUE),
  canary_pacific$detailed_source_url == 'https://repository.library.noaa.gov/view/noaa/72904',
  lookup('DFO:Pacific:Canary Rockfish')$matched_stock_id == 'dfo_canary_bc',
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_canary_rockfish_pacific_coast.md')),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(progress$inventory_rows) == 905L,
  sum(progress$eligible_new) == 30L,
  sum(progress$screening_complete) == 369L
)
cat('NOAA Pacific Coast Canary Rockfish identity and source-screening checks passed.\n')
english_sole <- lookup('NOAA:10046')
stopifnot(
  nrow(english_sole) == 1L,
  english_sole$inventory_assessment_year == 2013L,
  is.na(english_sole$candidate_unit_id),
  english_sole$match_status == 'new',
  is.na(english_sole$matched_stock_id),
  english_sole$screening_status == 'screened',
  english_sole$eligibility == 'ineligible',
  english_sole$import_status == 'not_started',
  grepl('catch-and-survey-index-only', english_sole$catch_age_evidence, fixed = TRUE),
  grepl('does not include age-specific survey observations', english_sole$survey_age_evidence, fixed = TRUE),
  grepl('2019 catch-only projection update', english_sole$output_evidence, fixed = TRUE),
  grepl('2005/2007', english_sole$decision_reason, fixed = TRUE),
  english_sole$detailed_source_url ==
    'https://www.pcouncil.org/documents/2015/01/data-moderate-stock-assessments-for-brown-china-copper-sharpchin-stripetail-and-yellowtail-rockfishes-and-english-and-rex-soles-in-2013-published-january-2015.pdf',
  file.exists(file.path('analysis', 'comp_assessments', 'discovery',
    'source_reviews', 'noaa_english_sole_pacific_coast.md')),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$excluded) == 63L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 89L,
  sum(progress$excluded) == 227L,
  sum(progress$screening_complete) == 369L
)
cat('NOAA Pacific Coast English Sole identity and source-screening checks passed.\n')
darkblotched <- lookup("NOAA:10562")
stopifnot(
  nrow(darkblotched) == 1L,
  darkblotched$inventory_assessment_year == 2025L,
  is.na(darkblotched$candidate_unit_id),
  darkblotched$match_status == "new",
  is.na(darkblotched$matched_stock_id),
  darkblotched$screening_status == "screened",
  darkblotched$eligibility == "eligible",
  darkblotched$import_status == "not_started",
  grepl("three fleets", darkblotched$catch_age_evidence, fixed = TRUE),
  grepl("three active biomass survey-index series", darkblotched$survey_age_evidence, fixed = TRUE),
  grepl("commented out", darkblotched$survey_age_evidence, fixed = TRUE),
  grepl("fixes female natural mortality and estimates male", darkblotched$mortality_evidence, fixed = TRUE),
  grepl("2017 update", darkblotched$decision_reason, fixed = TRUE),
  grepl("2025 Stock SMART inventory year refers to a catch-only projection", darkblotched$decision_reason, fixed = TRUE),
  darkblotched$detailed_source_url == "https://www.pcouncil.org/documents/2018/01/status-of-the-darkblotched-rockfish-resource-off-the-continental-u-s-pacific-coast-in-2017-update-of-2015-assessment-model.pdf",
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews", "noaa_darkblotched_rockfish_pacific_coast.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$eligible) == 37L,
  sum(noaa_progress$identity_unreviewed) == 351L,
  sum(noaa_progress$eligibility_unresolved) == 332L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L,
  sum(progress$new) == 89L,
  sum(progress$eligible) == 99L,
  sum(progress$eligible_new) == 30L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("NOAA Pacific Coast Darkblotched Rockfish identity and source-screening checks passed.\n")
goathorny <- lookup("NOAA:10457")
stopifnot(
  nrow(goathorny) == 1L,
  goathorny$inventory_assessment_year == 2024L,
  is.na(goathorny$candidate_unit_id),
  goathorny$match_status == "new",
  goathorny$screening_status == "screened",
  goathorny$eligibility == "ineligible",
  goathorny$import_status == "not_started",
  grepl("does not use catch-at-age", goathorny$catch_age_evidence, fixed = TRUE),
  grepl("regional bottom-trawl survey biomass estimates", goathorny$survey_age_evidence, fixed = TRUE),
  grepl("Natural mortality is assumed to be 0.03", goathorny$mortality_evidence, fixed = TRUE),
  grepl("59,459 t", goathorny$output_evidence, fixed = TRUE),
  grepl("age-structured assessment was suspended in 2003", goathorny$decision_reason, fixed = TRUE),
  grepl("thornyhead complex", goathorny$notes, fixed = TRUE),
  goathorny$detailed_source_url == "https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/GOAthorny.pdf",
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews", "noaa_goa_thornyhead_complex_2024.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$excluded) == 63L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(progress$new) == 89L,
  sum(progress$excluded) == 227L,
  sum(progress$screening_complete) == 369L
)
cat("NOAA Gulf of Alaska Thornyhead current-method screening checks passed.\n")
shortspine <- lookup("NOAA:10458")
stopifnot(
  nrow(shortspine) == 1L,
  shortspine$inventory_assessment_year == 2025L,
  is.na(shortspine$candidate_unit_id),
  shortspine$match_status == "new",
  shortspine$screening_status == "screened",
  shortspine$eligibility == "ineligible",
  shortspine$import_status == "not_started",
  grepl("does not use catch-at-age or age-composition data", shortspine$catch_age_evidence, fixed = TRUE),
  grepl("no age-composition input", shortspine$survey_age_evidence, fixed = TRUE),
  grepl("fixed at 0.04 per year for both sexes", shortspine$mortality_evidence, fixed = TRUE),
  grepl("July (month 7)", shortspine$timing_evidence, fixed = TRUE),
  grepl("derived age-specific N, catch and F outputs", shortspine$output_evidence, fixed = TRUE),
  grepl("2025 update is a catch-only projection", shortspine$decision_reason, fixed = TRUE),
  shortspine$detailed_source_url ==
    "https://www.pcouncil.org/documents/2024/02/status-of-shortspine-thornyhead-sebastolobus-alascanus-along-the-us-west-coast-in-2023.pdf; https://www.pcouncil.org/documents/shortspine-thornyhead-2023-model-files-zip-file/; https://www.pcouncil.org/documents/2025/08/g-7-attachment-3-catch-only-projection-for-shortspine-thornyhead.pdf/",
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_shortspine_thornyhead_pacific_coast.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$excluded) == 63L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L,
  sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("NOAA Pacific Coast Shortspine Thornyhead identity and source-screening checks passed.\n")
longspine <- lookup("NOAA:10460")
stopifnot(
  nrow(longspine) == 1L,
  longspine$inventory_assessment_year == 2019L,
  is.na(longspine$candidate_unit_id),
  longspine$match_status == "new",
  longspine$screening_status == "screened",
  longspine$eligibility == "ineligible",
  longspine$import_status == "not_started",
  grepl("no age data are available", longspine$catch_age_evidence, fixed = TRUE),
  grepl("no age-composition input", longspine$survey_age_evidence, fixed = TRUE),
  grepl("fixes M at 0.1113 per year", longspine$mortality_evidence, fixed = TRUE),
  grepl("2019 update is catch-only", longspine$decision_reason, fixed = TRUE),
  grepl("DFO inventory also lists Longspine Thornyhead", longspine$notes, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_longspine_thornyhead_pacific_coast.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$excluded) == 63L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L,
  sum(progress$inventory_rows) == 905L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
california_scorpionfish <- lookup("NOAA:10461")
stopifnot(
  nrow(california_scorpionfish) == 1L,
  california_scorpionfish$match_status == "new",
  is.na(california_scorpionfish$matched_stock_id),
  california_scorpionfish$screening_status == "screened",
  california_scorpionfish$eligibility == "eligible",
  california_scorpionfish$import_status == "not_started",
  grepl("CAAL data from the NWFSC trawl survey (2005-2016)", california_scorpionfish$catch_age_evidence, fixed = TRUE),
  grepl("no fishery age compositions", california_scorpionfish$catch_age_evidence, fixed = TRUE),
  grepl("ages 0-21, 2005-2016", california_scorpionfish$survey_age_evidence, fixed = TRUE),
  grepl("fixed at 0.235 per year", california_scorpionfish$mortality_evidence, fixed = TRUE),
  grepl("July (month 7; year fraction 0.5)", california_scorpionfish$timing_evidence, fixed = TRUE),
  grepl("F-at-age is derived as Z minus fixed M", california_scorpionfish$output_evidence, fixed = TRUE),
  grepl("2017-california-scorpionfish-model-files", california_scorpionfish$detailed_source_url, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_california_scorpionfish_southern_california_2017.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$eligible) == 37L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(progress$new) == 89L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligible) == 99L,
  sum(progress$screening_complete) == 369L
)
cat("NOAA California Scorpionfish identity and source-screening checks passed.\n")
cat("NOAA Pacific Coast Longspine Thornyhead identity and source-screening checks passed.\n")
kelp_greenling <- lookup("NOAA:10462")
stopifnot(
  nrow(kelp_greenling) == 1L,
  kelp_greenling$match_status == "new",
  is.na(kelp_greenling$matched_stock_id),
  kelp_greenling$screening_status == "screened",
  kelp_greenling$eligibility == "eligible",
  kelp_greenling$import_status == "not_started",
  grepl("commercial age compositions (2003-2013)", kelp_greenling$catch_age_evidence, fixed = TRUE),
  grepl("recreational ocean-boat compositions (2005-2013)", kelp_greenling$catch_age_evidence, fixed = TRUE),
  grepl("Three fishery-dependent aggregate CPUE indices", kelp_greenling$survey_age_evidence, fixed = TRUE),
  grepl("No fishery-independent survey", kelp_greenling$survey_age_evidence, fixed = TRUE),
  grepl("0.360 per year for females and 0.318 for males", kelp_greenling$mortality_evidence, fixed = TRUE),
  grepl("year fraction 0.5", kelp_greenling$timing_evidence, fixed = TRUE),
  grepl("316 mt", kelp_greenling$output_evidence, fixed = TRUE),
  grepl("age-grid discrepancy", kelp_greenling$decision_reason, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_oregon_kelp_greenling_2015.md")),
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$eligible) == 37L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligible) == 99L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("NOAA Oregon Kelp Greenling identity and source-screening checks passed.\n")
rex_sole <- register[!is.na(register$candidate_unit_id) &
  register$candidate_unit_id == "noaa_goa_rex_sole", , drop = FALSE]
rex_pacific <- lookup("NOAA:10059")
dfo_rebs <- lookup("DFO:Pacific:Rougheye Rockfish")
stopifnot(
  identical(sort(rex_sole$candidate_id), sort(c("NOAA:10058", "NOAA:15410", "NOAA:15411"))),
  nrow(rex_sole) == 3L,
  all(rex_sole$scientific_name == "Glyptocephalus zachirus"),
  all(rex_sole$match_status == "new"),
  all(rex_sole$screening_status == "screened"),
  all(rex_sole$eligibility == "eligible"),
  all(rex_sole$import_status == "not_started"),
  all(grepl("age compositions", rex_sole$catch_age_evidence, fixed = TRUE)),
  all(grepl("conditional age-at-length", rex_sole$survey_age_evidence, fixed = TRUE)),
  all(grepl("0.17 per year", rex_sole$mortality_evidence, fixed = TRUE)),
  all(grepl("Appendix F", rex_sole$output_evidence, fixed = TRUE)),
  all(grepl("2026", rex_sole$decision_reason, fixed = TRUE)),
  is.na(rex_sole$matched_stock_id),
  nrow(rex_pacific) == 1L,
  rex_pacific$candidate_unit_id == "noaa_pacific_rex_sole",
  rex_pacific$match_status == "new",
  rex_pacific$screening_status == "screened",
  rex_pacific$eligibility == "eligible",
  rex_pacific$import_status == "not_started",
  grepl("no fishery age-composition", rex_pacific$catch_age_evidence, fixed = TRUE),
  grepl("2007-2019", rex_pacific$survey_age_evidence, fixed = TRUE),
  grepl("numbers-at-age", rex_pacific$output_evidence, fixed = TRUE),
  !any(stocks$scientific_name == "Glyptocephalus zachirus"),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_goa_rex_sole_2026.md")),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_pacific_rex_sole_2023.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(progress$new) == 89L,
  sum(progress$eligible_new) == 30L,
  sum(progress$screening_complete) == 369L
)
cat("NOAA rex sole identity and source-screening checks passed.\n")
stopifnot(
  nrow(dfo_rebs) == 1L,
  dfo_rebs$match_status == "ambiguous",
  is.na(dfo_rebs$candidate_unit_id),
  is.na(dfo_rebs$official_stock_id),
  dfo_rebs$screening_status == "screened",
  dfo_rebs$eligibility == "unresolved",
  dfo_rebs$inventory_assessment_year == 2020,
  dfo_rebs$scientific_name == "Sebastes aleutianus; Sebastes melanostictus",
  grepl("weighted age proportions", dfo_rebs$catch_age_evidence, fixed = TRUE),
  grepl("WCHG", dfo_rebs$survey_age_evidence, fixed = TRUE),
  grepl("0.035, 0.045 and 0.055", dfo_rebs$mortality_evidence, fixed = TRUE),
  grepl("machine-readable run object", dfo_rebs$output_evidence, fixed = TRUE),
  grepl("two current BC units", dfo_rebs$decision_reason, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "dfo_rougheye_blackspotted_bc_2020.md")),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$ambiguous) == 5L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$ambiguous) == 10L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("DFO Rougheye/Blackspotted identity and age-data screening checks passed.\n")

anchovy <- lookup("NOAA:10574")
historical_iberian_anchovy <- stocks[stocks$stock_id == "ck_ices_wghansa_mordax_sp9a", , drop = FALSE]
stopifnot(
  nrow(anchovy) == 1L,
  anchovy$candidate_unit_id == "noaa_northern_anchovy_central_subpopulation",
  anchovy$match_status == "new",
  is.na(anchovy$matched_stock_id),
  anchovy$screening_status == "screened",
  anchovy$eligibility == "eligible",
  anchovy$import_status == "not_started",
  grepl("Southern Pacific Coast", anchovy$assessment_area, fixed = TRUE),
  grepl("not direct catch numbers-at-age", anchovy$catch_age_evidence, fixed = TRUE),
  grepl("lognormal errors", anchovy$survey_age_evidence, fixed = TRUE),
  grepl("0.726, 0.931, 0.985 and 0.997", anchovy$maturity_evidence, fixed = TRUE),
  grepl("0.414 per year", anchovy$mortality_evidence, fixed = TRUE),
  grepl("2015-2021", anchovy$timing_evidence, fixed = TRUE),
  grepl("Tables 11-14", anchovy$output_evidence, fixed = TRUE),
  grepl("No numerical data have been imported", anchovy$decision_reason, fixed = TRUE),
  nrow(historical_iberian_anchovy) == 1L,
  historical_iberian_anchovy$scientific_name == "Engraulis encrasicolus",
  !any(stocks$scientific_name == "Engraulis mordax"),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_northern_anchovy_central_subpopulation_2021.md")),
  sum(noaa_progress$new) == 50L,
  sum(noaa_progress$eligible_new) == 23L,
  sum(noaa_progress$eligible) == 37L,
  sum(noaa_progress$screening_complete) == 116L,
  sum(noaa_progress$screening_not_started) == 316L,
  sum(progress$new) == 89L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligible) == 99L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("NOAA northern anchovy identity and source-screening checks passed.\n")

red_porgy <- lookup("NOAA:10100")
stopifnot(
  nrow(red_porgy) == 1L,
  red_porgy$candidate_unit_id == "noaa_south_atlantic_red_porgy",
  red_porgy$match_status == "new",
  is.na(red_porgy$matched_stock_id),
  red_porgy$screening_status == "screened",
  red_porgy$eligibility == "eligible",
  red_porgy$import_status == "not_started",
  grepl("1,000 lb whole weight", red_porgy$catch_age_evidence, fixed = TRUE),
  grepl("pooled length composition only", red_porgy$catch_age_evidence, fixed = TRUE),
  grepl("1973-1998", red_porgy$survey_age_evidence, fixed = TRUE),
  grepl("Table 15", red_porgy$weights_evidence, fixed = TRUE),
  grepl("0.167", red_porgy$timing_evidence, fixed = TRUE),
  grepl("Table 16", red_porgy$output_evidence, fixed = TRUE),
  grepl("digitizing published composition figures", red_porgy$decision_reason, fixed = TRUE),
  !any(stocks$scientific_name == "Pagrus pagrus"),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "noaa_south_atlantic_red_porgy_2020.md"))
)
cat("NOAA South Atlantic red porgy source-screening checks passed.\n")

halibut <- lookup("DFO:Maritimes:Atlantic Halibut - NAFO 3NOPs4VWX5Z")
stopifnot(
  nrow(halibut) == 1L,
  halibut$candidate_unit_id == "dfo_atlantic_halibut_3nops4vwx5zc",
  halibut$scientific_name == "Hippoglossus hippoglossus",
  halibut$match_status == "new",
  is.na(halibut$matched_stock_id),
  halibut$screening_status == "screened",
  halibut$eligibility == "ineligible",
  halibut$import_status == "not_started",
  grepl("catch-at-length", halibut$catch_age_evidence, fixed = TRUE),
  grepl("A50=11.5", halibut$maturity_evidence, fixed = TRUE),
  grepl("time-varying natural mortality", halibut$mortality_evidence, fixed = TRUE),
  grepl("2025/003", halibut$notes, fixed = TRUE),
  grepl("one-to-one", halibut$notes, fixed = TRUE),
  !any(stocks$scientific_name == "Hippoglossus hippoglossus"),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "dfo_atlantic_halibut_3nops4vwx5zc_2022.md"))
)
cat("DFO Atlantic halibut unit and catch-at-length screen checks passed.\n")

halibut_4rst <- lookup("DFO:Quebec:Atlantic Halibut - NAFO 4RST")
stopifnot(
  nrow(halibut_4rst) == 1L,
  halibut_4rst$candidate_unit_id == "dfo_atlantic_halibut_4rst",
  halibut_4rst$scientific_name == "Hippoglossus hippoglossus",
  halibut_4rst$inventory_assessment_year == 2024,
  halibut_4rst$match_status == "new",
  halibut_4rst$screening_status == "screened",
  halibut_4rst$eligibility == "ineligible",
  halibut_4rst$import_status == "not_started",
  grepl("delay-difference model", halibut_4rst$catch_age_evidence, fixed = TRUE),
  grepl("no age-specific index", halibut_4rst$survey_age_evidence, fixed = TRUE),
  grepl("2,031 trawl-survey fish", halibut_4rst$weights_evidence, fixed = TRUE),
  grepl("not run-specific documentation", halibut_4rst$maturity_evidence, fixed = TRUE),
  grepl("natural-mortality series", halibut_4rst$mortality_evidence, fixed = TRUE),
  grepl("July-September", halibut_4rst$timing_evidence, fixed = TRUE),
  grepl("exploitable biomass", halibut_4rst$output_evidence, fixed = TRUE),
  grepl("2024 assessment is a delay-difference model", halibut_4rst$decision_reason, fixed = TRUE),
  grepl("assessment-date field lists 2023", halibut_4rst$notes, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "dfo_atlantic_halibut_4rst_2024.md")),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("DFO 4RST halibut source-screening checks passed.\n")

greenland_halibut_4rst <- lookup("DFO:Quebec:Greenland Halibut - NAFO 4RST")
stopifnot(
  nrow(greenland_halibut_4rst) == 1L,
  greenland_halibut_4rst$candidate_unit_id == "dfo_greenland_halibut_4rst",
  greenland_halibut_4rst$scientific_name == "Reinhardtius hippoglossoides",
  greenland_halibut_4rst$inventory_assessment_year == 2025,
  greenland_halibut_4rst$match_status == "new",
  is.na(greenland_halibut_4rst$matched_stock_id),
  greenland_halibut_4rst$screening_status == "screened",
  greenland_halibut_4rst$eligibility == "ineligible",
  greenland_halibut_4rst$import_status == "not_started",
  grepl("no catch-at-age", greenland_halibut_4rst$catch_age_evidence, fixed = TRUE),
  grepl("size class, not age", greenland_halibut_4rst$survey_age_evidence, fixed = TRUE),
  grepl("no run-specific weight-at-age", greenland_halibut_4rst$weights_evidence, fixed = TRUE),
  grepl("no age-specific maturity input", greenland_halibut_4rst$maturity_evidence, fixed = TRUE),
  grepl("does not estimate an age-structured natural-mortality", greenland_halibut_4rst$mortality_evidence, fixed = TRUE),
  grepl("August", greenland_halibut_4rst$timing_evidence, fixed = TRUE),
  grepl("7,320 t", greenland_halibut_4rst$output_evidence, fixed = TRUE),
  grepl("age-1 recruitment index", greenland_halibut_4rst$output_evidence, fixed = TRUE),
  grepl("2026/013", greenland_halibut_4rst$notes, fixed = TRUE),
  file.exists(file.path("analysis", "comp_assessments", "discovery", "source_reviews",
    "dfo_greenland_halibut_4rst_2025_update.md")),
  !any(stocks$scientific_name == "Reinhardtius hippoglossoides"),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat("DFO 4RST Greenland halibut source-screening checks passed.\n")

greenland_halibut_2_3klmno <- lookup('DFO:Newfoundland and Labrador:Greenland Halibut (Turbot) - NAFO 2 + 3KLMNO')
stopifnot(
  nrow(greenland_halibut_2_3klmno) == 1L,
  greenland_halibut_2_3klmno$candidate_unit_id == 'dfo_greenland_halibut_2_3klmno',
  greenland_halibut_2_3klmno$scientific_name == 'Reinhardtius hippoglossoides',
  greenland_halibut_2_3klmno$inventory_assessment_year == 2025,
  greenland_halibut_2_3klmno$match_status == 'new',
  is.na(greenland_halibut_2_3klmno$matched_stock_id),
  greenland_halibut_2_3klmno$screening_status == 'screened',
  greenland_halibut_2_3klmno$eligibility == 'eligible',
  greenland_halibut_2_3klmno$import_status == 'not_started',
  grepl('1975-2021', greenland_halibut_2_3klmno$catch_age_evidence, fixed = TRUE),
  grepl('mean numbers per tow', greenland_halibut_2_3klmno$survey_age_evidence, fixed = TRUE),
  grepl('Catch weights-at-age', greenland_halibut_2_3klmno$weights_evidence, fixed = TRUE),
  grepl('proportion-mature-at-age', greenland_halibut_2_3klmno$maturity_evidence, fixed = TRUE),
  grepl('M = 0.12', greenland_halibut_2_3klmno$mortality_evidence, fixed = TRUE),
  grepl('fraction of year', greenland_halibut_2_3klmno$timing_evidence, fixed = TRUE),
  grepl('14,729 t', greenland_halibut_2_3klmno$output_evidence, fixed = TRUE),
  grepl('not the 2025 production advice fit', greenland_halibut_2_3klmno$decision_reason, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_greenland_halibut_2_3klmno_2025.md')),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligible) == 99L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO 2+3KLMNO Greenland halibut source-screening checks passed.\n')
herring_4rsw_spring <- lookup('DFO:Newfoundland and Labrador:Atlantic Herring - NAFO 4R (Spring Spawner)')
herring_4rsw_fall <- lookup('DFO:Newfoundland and Labrador:Atlantic Herring - NAFO 4R (Fall Spawner)')
herring_4rsw <- rbind(herring_4rsw_spring, herring_4rsw_fall)
stopifnot(
  nrow(herring_4rsw) == 2L,
  identical(herring_4rsw$candidate_unit_id,
    c('dfo_herring_4rsw_spring', 'dfo_herring_4rsw_fall')),
  identical(herring_4rsw$scientific_name, rep('Clupea harengus', 2L)),
  all(herring_4rsw$assessment_area == 'NAFO divisions 4RSw (4R and 4Sw; west coast Newfoundland and Lower North Shore Quebec)'),
  all(herring_4rsw$inventory_assessment_year == 2024),
  all(herring_4rsw$match_status == 'new'),
  all(is.na(herring_4rsw$matched_stock_id)),
  all(herring_4rsw$screening_status == 'screened'),
  all(herring_4rsw$eligibility == 'unresolved'),
  all(herring_4rsw$import_status == 'not_started'),
  all(grepl('1966-2024', herring_4rsw$catch_age_evidence, fixed = TRUE)),
  all(grepl('2019-2024', herring_4rsw$survey_age_evidence, fixed = TRUE)),
  all(grepl('component-specific weight-at-age', herring_4rsw$weights_evidence, fixed = TRUE)),
  all(grepl('assumed constant', herring_4rsw$mortality_evidence, fixed = TRUE)),
  all(grepl('within-year timing', herring_4rsw$timing_evidence, fixed = TRUE)),
  all(grepl('age-2 recruitment', herring_4rsw$output_evidence, fixed = TRUE)),
  grepl('2025/065', herring_4rsw_spring$catch_age_evidence, fixed = TRUE),
  grepl('2025/065', herring_4rsw_fall$catch_age_evidence, fixed = TRUE),
  all(grepl('WHAM state-space', herring_4rsw$decision_reason, fixed = TRUE)),
  all(grepl('4RSw', herring_4rsw$notes, fixed = TRUE)),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_herring_4rsw_2024.md')),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO 4RSw spring/fall herring source-screening checks passed.\n')
yellowtail_3lno <- lookup('DFO:Newfoundland and Labrador:Yellowtail Flounder - NAFO 3LNO')
stopifnot(
  nrow(yellowtail_3lno) == 1L,
  yellowtail_3lno$candidate_unit_id == 'dfo_yellowtail_flounder_3lno',
  yellowtail_3lno$scientific_name == 'Limanda ferruginea',
  yellowtail_3lno$assessment_area == 'NAFO divisions 3LNO',
  yellowtail_3lno$inventory_assessment_year == 2025,
  yellowtail_3lno$match_status == 'new',
  is.na(yellowtail_3lno$matched_stock_id),
  yellowtail_3lno$screening_status == 'screened',
  yellowtail_3lno$eligibility == 'ineligible',
  yellowtail_3lno$import_status == 'not_started',
  grepl('No catch-at-age', yellowtail_3lno$catch_age_evidence, fixed = TRUE),
  grepl('No age-specific survey index', yellowtail_3lno$survey_age_evidence, fixed = TRUE),
  grepl('No weight-at-age', yellowtail_3lno$weights_evidence, fixed = TRUE),
  grepl('No maturity-at-age', yellowtail_3lno$maturity_evidence, fixed = TRUE),
  grepl('no age-specific natural mortality', yellowtail_3lno$mortality_evidence, fixed = TRUE),
  grepl('B/BMSY', yellowtail_3lno$output_evidence, fixed = TRUE),
  yellowtail_3lno$detailed_source_url == 'https://www.nafo.int/Portals/0/PDFs/sc/2025/scr25-026.pdf',
  yellowtail_3lno$review_date == '2026-10-29',
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_yellowtail_flounder_3lno_2025.md')),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO 3LNO yellowtail flounder source-screening checks passed.
')
witch_3no <- lookup('DFO:Newfoundland and Labrador:Witch Flounder - NAFO 3NO')
stopifnot(
  nrow(witch_3no) == 1L,
  witch_3no$candidate_unit_id == 'dfo_witch_flounder_3no',
  witch_3no$scientific_name == 'Glyptocephalus cynoglossus',
  witch_3no$assessment_area == 'NAFO divisions 3NO',
  witch_3no$inventory_assessment_year == 2024,
  witch_3no$match_status == 'new',
  is.na(witch_3no$matched_stock_id),
  witch_3no$screening_status == 'screened',
  witch_3no$eligibility == 'ineligible',
  witch_3no$import_status == 'not_started',
  grepl('No catch-at-age', witch_3no$catch_age_evidence, fixed = TRUE),
  grepl('No age-specific survey index', witch_3no$survey_age_evidence, fixed = TRUE),
  grepl('Age-specific weights are not inputs', witch_3no$weights_evidence, fixed = TRUE),
  grepl('Age-specific maturity is not an input', witch_3no$maturity_evidence, fixed = TRUE),
  grepl('No age-specific natural mortality', witch_3no$mortality_evidence, fixed = TRUE),
  grepl('B/BMSY', witch_3no$output_evidence, fixed = TRUE),
  witch_3no$detailed_source_url == 'https://www.nafo.int/Portals/0/PDFs/sc/2024/scr24-018.pdf',
  witch_3no$review_date == '2026-10-29',
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_witch_flounder_3no_2024.md')),
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$new) == 89L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO NAFO 3NO witch flounder source-screening checks passed.
')
yellowtail_5z <- lookup('DFO:Maritimes:Yellowtail Flounder - NAFO 5Z')
stopifnot(
  nrow(yellowtail_5z) == 1L,
  yellowtail_5z$scientific_name == 'Limanda ferruginea',
  yellowtail_5z$assessment_area == 'NAFO Division 5Z (Georges Bank)',
  yellowtail_5z$inventory_assessment_year == 2025,
  yellowtail_5z$match_status == 'matched',
  yellowtail_5z$matched_stock_id == 'ck_dfo_nefsc_ferruginea_gb',
  yellowtail_5z$screening_status == 'screened',
  yellowtail_5z$eligibility == 'unresolved',
  yellowtail_5z$import_status == 'imported',
  grepl('run-matched catch-at-age', yellowtail_5z$catch_age_evidence, fixed = TRUE),
  grepl('NEFSC spring/fall and DFO spring', yellowtail_5z$survey_age_evidence, fixed = TRUE),
  grepl('age-based estimated M', yellowtail_5z$mortality_evidence, fixed = TRUE),
  grepl('Aggregate catch, SSB', yellowtail_5z$output_evidence, fixed = TRUE),
  grepl('NOAA:10033', yellowtail_5z$notes, fixed = TRUE),
  grepl('Georges-Bank-Yellowtail-Flounder-Assessment-Update-2025.pdf', yellowtail_5z$detailed_source_url, fixed = TRUE),
  lookup('NOAA:10033')$matched_stock_id == yellowtail_5z$matched_stock_id,
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_yellowtail_5z_crosswalk.md')),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$matched) == 105L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO 5Z yellowtail crosswalk and source-screening checks passed.
')
nl_herring <- lookup('DFO:Newfoundland and Labrador:Atlantic Herring - NAFO 2J3KLPs')
assessments <- read.csv('analysis/comp_assessments/database/assessments.csv', stringsAsFactors = FALSE)
stopifnot(
  nrow(nl_herring) == 1L,
  nl_herring$official_stock_id == 'NAFO 2J3KLPs',
  nl_herring$scientific_name == 'Clupea harengus',
  nl_herring$inventory_assessment_year == 2025,
  is.na(nl_herring$candidate_unit_id),
  nl_herring$match_status == 'ambiguous',
  nl_herring$screening_status == 'screened',
  nl_herring$eligibility == 'ineligible',
  nl_herring$import_status == 'not_started',
  grepl('annual commercial catch numbers-at-age', nl_herring$catch_age_evidence, fixed = TRUE),
  grepl('age-4 recruitment index', nl_herring$survey_age_evidence, fixed = TRUE),
  grepl('not age-specific abundance indices', nl_herring$survey_age_evidence, fixed = TRUE),
  grepl('survey-level mean weights', nl_herring$weights_evidence, fixed = TRUE),
  grepl('maturity-at-age surface', nl_herring$maturity_evidence, fixed = TRUE),
  grepl('no population-level', nl_herring$mortality_evidence, fixed = TRUE),
  grepl('April 1 and July 31', nl_herring$timing_evidence, fixed = TRUE),
  grepl('no fitted', nl_herring$output_evidence, fixed = TRUE),
  grepl('41311899.pdf', nl_herring$detailed_source_url, fixed = TRUE),
  nl_herring$review_date == '2026-10-29',
  file.exists('analysis/comp_assessments/discovery/source_reviews/dfo_nl_herring_3klps_2025.md'),
  !any(grepl('2J3KLPs', assessments$assessment_id, fixed = TRUE)),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$new) == 21L,
  sum(dfo_progress$ambiguous) == 5L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$matched) == 105L,
  sum(progress$new) == 89L,
  sum(progress$ambiguous) == 10L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligible) == 99L,
  sum(progress$excluded) == 227L,
  sum(progress$eligible_new) == 30L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO Newfoundland 2J3KLPs herring source-screening checks passed.\n')

dfo_4vwx_herring_components <- register[register$candidate_id %in% c(
  'DFO:Maritimes:Atlantic Herring - Coastal Nova Scotia',
  'DFO:Maritimes:Atlantic Herring - Offshore Scotian Shelf',
  'DFO:Maritimes:Atlantic Herring - Southwest Nova Scotia / Bay of Fundy'
), , drop = FALSE]
sw_nb_migrants <- lookup(
  'DFO:Maritimes:Atlantic Herring - Southwest New Brunswick (migrant juveniles)'
)
stopifnot(
  nrow(dfo_4vwx_herring_components) == 3L,
  all(dfo_4vwx_herring_components$candidate_unit_id == 'dfo_herring_4vwx_management_unit'),
  all(dfo_4vwx_herring_components$matched_stock_id == 'ck_dfo_mar_harengus_4vwx'),
  all(dfo_4vwx_herring_components$match_status == 'matched'),
  all(dfo_4vwx_herring_components$screening_status == 'screened'),
  all(dfo_4vwx_herring_components$eligibility == 'unresolved'),
  all(dfo_4vwx_herring_components$import_status == 'imported'),
  nrow(sw_nb_migrants) == 1L,
  sw_nb_migrants$match_status == 'ambiguous',
  sw_nb_migrants$screening_status == 'screened',
  sw_nb_migrants$eligibility == 'ineligible',
  sw_nb_migrants$import_status == 'imported',
  grepl('not indicative of abundance', sw_nb_migrants$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_4vwx_herring_component_crosswalk_2025.md')),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$ambiguous) == 5L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$matched) == 105L,
  sum(progress$ambiguous) == 10L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO 4VWX herring component crosswalk checks passed.\n')

dfo_mackerel <- lookup('DFO:National Capital Region:Atlantic Mackerel')
stopifnot(
  nrow(dfo_mackerel) == 1L,
  dfo_mackerel$candidate_unit_id == 'dfo_northern_contingent_atlantic_mackerel',
  dfo_mackerel$scientific_name == 'Scomber scombrus',
  dfo_mackerel$matched_stock_id == 'ck_dfo_qc_scombrus_4t',
  dfo_mackerel$match_status == 'matched',
  dfo_mackerel$screening_status == 'screened',
  dfo_mackerel$eligibility == 'eligible',
  dfo_mackerel$import_status == 'imported',
  grepl('ages 1-10+', dfo_mackerel$catch_age_evidence, fixed = TRUE),
  grepl('aggregate total egg production index', dfo_mackerel$survey_age_evidence, fixed = TRUE),
  grepl('0.28 per year', dfo_mackerel$mortality_evidence, fixed = TRUE),
  grepl('N-at-age and F-at-age', dfo_mackerel$output_evidence, fixed = TRUE),
  lookup('NOAA:10164')$matched_stock_id == 'nefsc_atlantic_mackerel',
  file.exists(file.path('analysis', 'comp_assessments', 'source_reviews',
    'dfo_northern_mackerel_2025.md')),
  sum(dfo_progress$matched) == 20L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$eligible) == 17L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$matched) == 105L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$eligible) == 99L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO Atlantic mackerel inventory crosswalk checks passed.\n')

quebec_4s_herring <- lookup('DFO:Quebec:Atlantic Herring - Herring Fishing Area 15')
stopifnot(
  nrow(quebec_4s_herring) == 1L,
  quebec_4s_herring$candidate_unit_id == 'dfo_quebec_4s_herring',
  quebec_4s_herring$official_stock_id == paste(
    'Quebec North Shore (NAFO Division 4S) spring- and fall-spawning stocks'),
  quebec_4s_herring$scientific_name == 'Clupea harengus',
  quebec_4s_herring$inventory_assessment_year == 2024,
  is.na(quebec_4s_herring$matched_stock_id),
  quebec_4s_herring$match_status == 'ambiguous',
  quebec_4s_herring$screening_status == 'screened',
  quebec_4s_herring$eligibility == 'ineligible',
  quebec_4s_herring$import_status == 'not_started',
  grepl('through 2024', quebec_4s_herring$catch_age_evidence, fixed = TRUE),
  grepl('age-specific abundances', quebec_4s_herring$survey_age_evidence, fixed = TRUE),
  grepl('unavailable', quebec_4s_herring$mortality_evidence, fixed = TRUE),
  grepl('do not provide', quebec_4s_herring$output_evidence, fixed = TRUE),
  file.exists(file.path('analysis', 'comp_assessments', 'discovery', 'source_reviews',
    'dfo_quebec_herring_4s_2025.md')),
  sum(dfo_progress$ambiguous) == 5L,
  sum(dfo_progress$identity_unreviewed) == 153L,
  sum(dfo_progress$excluded) == 113L,
  sum(dfo_progress$eligibility_unresolved) == 69L,
  sum(dfo_progress$screening_complete) == 144L,
  sum(dfo_progress$screening_not_started) == 55L,
  sum(progress$ambiguous) == 10L,
  sum(progress$identity_unreviewed) == 701L,
  sum(progress$excluded) == 227L,
  sum(progress$eligibility_unresolved) == 579L,
  sum(progress$screening_complete) == 369L,
  sum(progress$screening_not_started) == 536L
)
cat('DFO Quebec North Shore herring screening checks passed.\n')
