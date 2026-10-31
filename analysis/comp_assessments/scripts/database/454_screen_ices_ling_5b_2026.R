root <- file.path('analysis', 'comp_assessments')
cache <- file.path(root, 'source_cache', 'discovery_ices_ling_5b_2026')
register_path <- file.path(root, 'discovery', 'candidates.csv')
candidate_id <- 'ICES:lin.27.5b'
unit_id <- 'ices_ling_faroe_5b'

sha256_file <- function(path) {
  output <- system2('certutil.exe',
    c('-hashfile', shQuote(normalizePath(path, winslash = '\\', mustWork = TRUE)), 'SHA256'),
    stdout = TRUE)
  hashes <- tolower(trimws(grep('^[[:space:]]*[[:xdigit:]]{64}[[:space:]]*$',
    output, value = TRUE)))
  if (length(hashes) != 1L) stop('Could not verify source checksum: ', path, call. = FALSE)
  hashes
}

expected <- c(
  WGDEEP_2026_v1.pdf = '4dece836b45b56c03a34d2ddf7ee1dd88eda7cdf2ba8d6717f86949888da89cd',
  WGDEEP_2026_ling_5b_section.txt = 'f5f68201da8f1874c901e73ba591dec0cf6bca4fe4ac167c26a9d749f7cc0e89',
  WKBARFAR_2021_report.pdf = '94a98c08f64ca05ea6a5d3db32bba4d29a0cc2e863580619cefd8b181d3b7949',
  ICES_advice_7205.html = 'dd4fc9f818ef445a311e196aa81e6a220efc2727fdf3f0ef89cdcca3b59b6187',
  ICES_standard_graphs_source_23735.html = 'c00314737346f7ca637ccec1bc29375158e08c0544bddb5966b1758882bff6c1',
  ICES_SAG_23735.xml = '992bdf680900755292bf527f5615848947e7c8e71887d42590fdfb917073a5dc'
)
manifest <- read.csv(file.path(cache, 'sources.csv'), stringsAsFactors = FALSE,
  check.names = FALSE)
if (nrow(manifest) != length(expected) || anyDuplicated(manifest$file) ||
    !setequal(manifest$file, names(expected)))
  stop('The Faroes ling source manifest is incomplete.', call. = FALSE)
for (file in names(expected)) {
  source_row <- match(file, manifest$file)
  path <- normalizePath(file.path(cache, manifest$local_path[source_row]),
    winslash = '/', mustWork = FALSE)
  if (!file.exists(path) || manifest$retrieved_on[source_row] != '2026-10-31' ||
      as.numeric(manifest$bytes[source_row]) != file.info(path)$size ||
      tolower(manifest$sha256[source_row]) != expected[[file]] ||
      sha256_file(path) != expected[[file]])
    stop('Cached Faroes ling source differs from its manifest: ', file, call. = FALSE)
}
read_source <- function(file) {
  row <- match(file, manifest$file)
  paste(readLines(file.path(cache, manifest$local_path[row]), warn = FALSE), collapse = ' ')
}
report <- read_source('WGDEEP_2026_ling_5b_section.txt')
advice <- read_source('ICES_advice_7205.html')
source_data <- read_source('ICES_standard_graphs_source_23735.html')
stock_xml <- read_source('ICES_SAG_23735.xml')

report_evidence <- c(
  'Catch-at-age data are available from the Faroese fishery in 5.b since 1996.',
  'Table 3.2.7. Ling in 5.b. Catch numbers at age (*1000) used in the assessment.',
  'Table 3.2.8. Ling in 5.b. Weighted mean weights at age used in the assessment.',
  'Table 3.2.9. Ling in 5.b. Spring survey input to the tuning series in the assessment.',
  'Table 3.2.10. Ling in 5.b. Summer survey input to tuning series in the assessment.',
  'Fixed proportion mature at age used in the assessment',
  'Natural mortality is set to 0.15 for all ages and years.',
  'spring survey takes place in Febru- ary/March',
  'summer survey is conducted in August',
  'final assessment lin.27.5b_wgdeep2026_final',
  'Table 3.2.11. Ling in 5.b. Estimated recruitment, spawning stock biomass (SSB), and average fishing mortality.',
  'Table 3.2.12. Ling in 5.b. Estimated fishing mortality at age.',
  'Table 3.2.13. Ling in 5.b. Estimated stock numbers at age.'
)
source_checks <- c(
  report = all(vapply(report_evidence, grepl, logical(1L), x = report, fixed = TRUE)),
  advice = all(vapply(c('lin.27.5b', 'Assessment Year', '2026', 'Stock Category', 'zero catch'),
    grepl, logical(1L), x = advice, fixed = TRUE)),
  standard_graphs = all(vapply(c('lin.27.5b', 'Assessment Year', '2026',
    'Unit weights', 'tonnes', 'Unit recruitment', 'thousands',
    'Stock size description', 'lblStockSizeUnits'),
    grepl, logical(1L), x = source_data, fixed = TRUE)),
  sag_xml = all(vapply(c('<AssessmentYear>2026</AssessmentYear>',
    '<StockCode>lin.27.5b</StockCode>', '<StockCategory>1</StockCategory>',
    '<ConfidenceIntervalDefinition>95%</ConfidenceIntervalDefinition>'),
    grepl, logical(1L), x = stock_xml, fixed = TRUE))
)
if (!all(source_checks)) {
  print(source_checks)
  print(report_evidence[!vapply(report_evidence, grepl, logical(1L), x = report, fixed = TRUE)])
  stop('The cached sources do not support the documented Faroes ling screen.', call. = FALSE)
}

candidates <- read.csv(register_path, stringsAsFactors = FALSE,
  na.strings = c('', 'NA'), check.names = FALSE)
stocks <- read.csv(file.path(root, 'database', 'stocks.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
assessments <- read.csv(file.path(root, 'database', 'assessments.csv'),
  stringsAsFactors = FALSE, na.strings = c('', 'NA'))
inputs <- read.csv(file.path(root, 'database', 'inputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
outputs <- read.csv(file.path(root, 'database', 'outputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
row <- match(candidate_id, candidates$candidate_id)
if (is.na(row) || sum(candidates$candidate_id == candidate_id, na.rm = TRUE) != 1L ||
    candidates$official_stock_id[row] != 'lin.27.5b' ||
    candidates$scientific_name[row] != 'Molva molva' ||
    candidates$inventory_assessment_year[row] != 2026L ||
    candidates$assessment_area[row] != 'Faroes' ||
    any(candidates$candidate_unit_id[-row] == unit_id, na.rm = TRUE) ||
    any(tolower(stocks$authority_stock_id) == 'lin.27.5b', na.rm = TRUE) ||
    any(stocks$stock_id == unit_id) ||
    any(unit_id %in% c(assessments$stock_id, assessments$assessment_id,
      inputs$assessment_id, outputs$assessment_id), na.rm = TRUE))
  stop('The Faroes ling identity conflicts with the register or canonical database.', call. = FALSE)

review <- list(
  candidate_unit_id = unit_id,
  match_status = 'new',
  screening_status = 'screened',
  eligibility = 'eligible',
  import_status = 'not_started',
  catch_age_evidence = paste('WGDEEP 2026 Table 3.2.7 publishes catch numbers-at-age for ages 3–12+ over 1996–2025, in thousands of fish. The report describes commercial Faroese catch sampling; 2025 longliner catch-at-age uses weight-grader data because no longliner length samples were available.'),
  survey_age_evidence = paste('WGDEEP 2026 Tables 3.2.9 and 3.2.10 publish numeric age-disaggregated spring-survey indices for ages 4–11 in 1998–2025 and summer-survey indices for ages 3–11 in 1996–2025, each with a survey-effort column. The 2026 report states that summer indices for 2021–2023 and spring indices for 2022–2024 were downweighted for the change in rockhopper gear, but does not publish numeric weights or label the absolute index unit or scaling; confirm these against the final SAM data before import.'),
  weights_evidence = paste('Table 3.2.8 publishes annual weighted mean catch weights-at-age for ages 3–12+ over 1996–2025. The table does not label the weight unit. The 2021 benchmark stock annex records the decision to use catch weights as stock weights instead of spring-survey weights; verify the native unit before import.'),
  maturity_evidence = paste('The 2026 report publishes fixed, time-invariant Faroese-survey maturity proportions for ages 3–12+: 0.00, 0.04, 0.19, 0.50, 0.79, 0.93, 0.98, 1.00, 1.00, and 1.00.'),
  mortality_evidence = 'Natural mortality is fixed at M = 0.15 for all modeled ages and years.',
  timing_evidence = paste('The Faroese spring survey is conducted in February/March and the summer survey in August. The reports do not state the exact model timing fractions; month-based sampling fractions are an approximation and should be checked against the final SAM data.'),
  output_evidence = paste('WGDEEP 2026 Tables 3.2.11–3.2.13 report annual recruitment at age 3, SSB, Fbar ages 6–10, total biomass, F-at-age, and stock numbers-at-age for 1996–2025. The aggregate series include 95% intervals; the age-specific surfaces are point estimates without age-specific uncertainty. ICES Standard Graphs identifies SSB in tonnes and recruitment in thousands.'),
  detailed_source_url = paste('https://doi.org/10.17895/ices.pub.32051352.v1; ',
    'https://asd.ices.dk/viewAdvice/7205; ',
    'https://standardgraphs.ices.dk/ViewSourceData.aspx?key=23735; ',
    'https://sg.ices.dk/SAG_API/StockXmlDownload/23735; ',
    'https://stockassessment.org/datadisk/stockassessment/userdirs/user3/lin.27.5b_wgdeep2026_final/; ',
    'https://doi.org/10.17895/ices.pub.7920'),
  review_date = '2026-10-31',
  decision_reason = paste('This is a distinct Faroes stock unit with a 2026 active Category 1 SAM assessment, separate from the Icelandic ling unit lin.27.5a. Current official report tables publish numeric catch-at-age, two age-disaggregated survey series, catch weights-at-age, fixed maturity and M, plus comparable SSB, recruitment and F outputs including N-at-age and F-at-age. It is eligible for source curation. Verify survey-index scale and effort treatment, downweight factors, and explicit weight units against the final SAM source before canonical import; do not infer or relabel units from table values.'),
  notes = paste('The 2026 WGDEEP report identifies the final run as lin.27.5b_wgdeep2026_final. It reports downweighting summer indices in 2021–2023 and spring indices in 2022–2024 after a rockhopper-gear change, but the numeric weights are not in the report. Exact survey fractional-year timing, index scale/effort treatment and weight units should be checked against the native source. The direct 2021 ICES stock-annex URL returned HTTP 401 during retrieval; annex material and ALK methods are included in the cached WKBARFAR 2021 report. Source files and checksums are under source_cache/discovery_ices_ling_5b_2026; the large 2026 report PDF is reused from source_cache/discovery_ices_gss_5b6a_2026. No canonical records were added.')
)

is_unreviewed <- candidates$match_status[row] == 'unreviewed' &&
  candidates$screening_status[row] == 'not_started' &&
  candidates$eligibility[row] == 'unresolved' &&
  (is.na(candidates$candidate_unit_id[row]) || !nzchar(candidates$candidate_unit_id[row]))
current <- lapply(names(review), function(column) as.character(candidates[[column]][row]))
names(current) <- names(review)
same_review <- all(vapply(names(review), function(column)
  identical(current[[column]], as.character(review[[column]])), logical(1L)))
if (!is_unreviewed && !same_review)
  stop('The Faroes ling candidate already has a conflicting review.', call. = FALSE)
if (!identical(Sys.getenv('TAM_SOURCE_DRY_RUN'), '1') && (is_unreviewed || !same_review)) {
  for (column in names(review)) candidates[[column]][row] <- review[[column]]
  source(file.path(root, 'scripts/database/candidate_register_io.R'))
  replace_candidate_register_row(register_path, row, candidates[row, , drop = FALSE],
    candidate_id)
}

cat('Screened ICES Faroes ling as a new, eligible unit for source curation.\n')
