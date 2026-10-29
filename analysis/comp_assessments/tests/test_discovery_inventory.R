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
  lookup('DFO:Gulf:American Plaice - NAFO 4T')$eligibility == 'unresolved',
  lookup('DFO:Pacific:Dungeness Crab')$eligibility == 'ineligible',
  lookup('DFO:Maritimes:Grey Seal')$eligibility == 'ineligible'
)
stopifnot(
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
noaa_progress <- progress[progress$organization == 'NOAA', , drop = FALSE]
stopifnot(sum(noaa_progress$new) == 9L, sum(noaa_progress$eligible_new) == 6L,
  sum(noaa_progress$identity_unreviewed) == 395L,
  sum(noaa_progress$eligible) == 20L, sum(noaa_progress$excluded) == 46L,
  sum(noaa_progress$eligibility_unresolved) == 366L,
  sum(noaa_progress$screening_complete) == 72L,
  sum(noaa_progress$screening_not_started) == 360L)
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
