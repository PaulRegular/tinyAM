pkgload::load_all('.', quiet = TRUE)
source('analysis/comp_assessments/R/run_assessment.R')
source('analysis/comp_assessments/scripts/translation/review_translations.R')
source('analysis/comp_assessments/tests/helper_stock.R')

db <- read_database()
current <- db$assessments$assessment_id[db$assessments$is_current %in% TRUE]
stopifnot(setequal(names(review_inventory), current), length(current) == 23L,
          all(names(review_candidates) %in% current))
for (ids in capability_evidence) {
  stopifnot(!anyDuplicated(ids), all(ids %in% current))
}
stopifnot(length(capability_evidence$correlated_F_increments) == 10L,
          length(capability_evidence$correlated_observations) == 9L,
          length(capability_evidence$latent_SD_groups) == 9L,
          'dfo_cod_3pn4rs_2025' %in% capability_evidence$latent_SD_groups,
          !'ices_cod_northeast_arctic_2026' %in% capability_evidence$correlated_F_increments,
          !'ices_bluewhiting_northeast_atlantic_2026' %in% capability_evidence$latent_SD_groups)

record <- read.csv('analysis/comp_assessments/results/translation_review.csv',
                   na.strings = c('', 'NA'))
stopifnot(nrow(record) == length(current) + sum(lengths(review_candidates)),
          setequal(record$assessment_id, current),
          !anyDuplicated(record[c('assessment_id', 'candidate')]),
          !any(record$decision == 'pending_review'),
          all(record$baseline_revision == review_revision),
          sum(record$decision == 'blocked') == 2L,
          all(is.na(record$recruitment_percent[record$recruitment_status %in% 'non_equivalent'])))

# The frozen model comparison must use the same observations, ages and years.
ids <- sub('[.]R$', '', list.files('analysis/comp_assessments/scripts/translation/stocks'))
for (id in ids) {
  original <- .review_recipe(id, db)
  now <- run_assessment(id, db, fit = FALSE)
  stopifnot(isTRUE(all.equal(original$obs, now$obs, check.environment = FALSE)))
  old <- do.call(tinyAM::prepare_tam, c(list(data = original$obs,
    years = original$years, ages = original$ages), original$settings))
  stock <- .test_stock(read_assessment(id, db))
  stopifnot(isTRUE(all.equal(old, stock$dat, check.environment = FALSE)))
  starts <- all.equal(original$start_par, stock$start_par, check.environment = FALSE)
  if (!isTRUE(starts)) stop(id, ': ', paste(starts, collapse = '; '))
}

fixture <- data.frame(
  assessment_id = 'fixture',
  component = c('recruitment', 'recruitment', 'recruitment', 'index', 'index', 'F', 'M', 'M'),
  setting = c('stock_recruitment', 'stock_recruitment', 'deviation_distribution',
              'observation_correlation', 'observation_correlation', 'process', 'process', 'process_sd'),
  value = c('Beverton-Holt', 'hockey-stick', 'AR1', 'Independent', 'AR1 across ages',
            'Random walk with age-correlated increments', 'AR1', 'Fixed 0.075'),
  source_reference = 'offline fixture', notes = ''
)
audit <- audit_assumptions('fixture', fixture)
stopifnot(identical(audit$tinyam_support,
  c('partially_supported', 'unsupported', 'partially_supported', 'supported',
    'unsupported', 'partially_supported', 'partially_supported', 'partially_supported')),
  grepl('supplied SD', audit$audit_notes[8]),
  !any(grepl('basic random walk|does not fit.*autocorrelated', audit$audit_notes)))
cat('Review inventory and current-capability audit tests passed.\n')
