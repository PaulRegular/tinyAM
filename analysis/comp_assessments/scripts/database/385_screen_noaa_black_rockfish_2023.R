root <- file.path('analysis', 'comp_assessments')
cache <- file.path(root, 'source_cache', 'discovery_noaa_black_rockfish_2023')
register_path <- file.path(root, 'discovery', 'candidates.csv')
candidate_id <- 'NOAA:166727'
unit_id <- 'noaa_black_rockfish_pacific_coast_2023'
source_files <- c(
  'Black_Rockfish_Washington_2023_assessment.pdf',
  'Black_Rockfish_Oregon_2023_assessment.pdf',
  'Black_Rockfish_California_2023_assessment.pdf',
  'Black_Rockfish_Washington_2023_model_files.zip',
  'Black_Rockfish_Oregon_2023_model_files.zip',
  'Black_Rockfish_California_2023_model_files.zip',
  'Black_Rockfish_2023_STAR_panel.pdf',
  'NOAA_StockSMART_changelog.html')
source_urls <- c(
  'https://www.pcouncil.org/documents/2024/02/status-of-black-rockfish-sebastes-melanops-in-2023-in-state-and-federal-waters-off-washington-state.pdf/',
  'https://www.pcouncil.org/documents/2024/02/status-of-black-rockfish-sebastes-melanops-in-2023-in-state-and-federal-waters-off-oregon.pdf/',
  'https://www.pcouncil.org/documents/g-2-situation-summary-adopt-stock-assessments/g-2-attachment-8-draft-executive-summary-of-status-of-black-rockfish-sebastes-melanops-in-u-s-waters-off-california-in-2023/',
  'https://www.pcouncil.org/documents/black-rockfish-wa-2023-model-files-zip-file/',
  'https://www.pcouncil.org/documents/black-rockfish-or-2023-model-files-zip-file/',
  'https://www.pcouncil.org/documents/black-rockfish-ca-2023-model-files-zip/',
  'https://www.pcouncil.org/documents/2023/11/stock-assessment-review-star-panel-report-for-black-rockfish.pdf/',
  'https://noaa-edab.github.io/stocksmart/news/index.html')
source_revisions <- c(
  'PFMC 2023 Washington assessment, published February 2024',
  'PFMC 2023 Oregon assessment, published February 2024',
  'PFMC G.2 supplemental revised Attachment 8, retrieved 2026-10-30',
  'PFMC Washington 2023 base model files ZIP',
  'PFMC Oregon 2023 base model files ZIP',
  'PFMC California 2023 North and Central base model files ZIP',
  'PFMC Black Rockfish 2023 STAR panel report',
  'NOAA-EDAB stocksmart changelog; 2023 assessment-ID record')
source_paths <- file.path(cache, source_files)

sha256_file <- function(path) {
  output <- system2('certutil.exe', c('-hashfile', shQuote(normalizePath(path)), 'SHA256'),
    stdout = TRUE)
  hashes <- trimws(grep('^[[:space:]]*[[:xdigit:]]{64}[[:space:]]*$', output, value = TRUE))
  if (length(hashes) != 1L) stop('Could not calculate SHA-256 for ', path, call. = FALSE)
  toupper(hashes)
}

if (any(!file.exists(source_paths)) || any(file.info(source_paths)$size < 1000L))
  stop('Cache the reviewed NOAA black rockfish sources before screening.', call. = FALSE)
sources <- read.csv(file.path(cache, 'sources.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'), fileEncoding = 'UTF-8-BOM')
source_rows <- match(source_files, basename(sources$file))
source_valid <- c(
  urls = identical(as.character(sources$url[source_rows]), source_urls),
  revisions = identical(as.character(sources$revision[source_rows]), source_revisions),
  hashes = identical(toupper(as.character(sources$sha256[source_rows])),
    unname(vapply(source_paths, sha256_file, character(1L)))),
  sizes = identical(as.numeric(sources$bytes[source_rows]), as.numeric(file.info(source_paths)$size)))
if (anyNA(source_rows) || !all(source_valid))
  stop('Cached NOAA black rockfish sources differ from the reviewed manifest: ',
    paste(names(source_valid)[!source_valid], collapse = ', '), call. = FALSE)

register <- read.csv(register_path, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
stocks <- read.csv(file.path(root, 'database', 'stocks.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
assessments <- read.csv(file.path(root, 'database', 'assessments.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
inputs <- read.csv(file.path(root, 'database', 'inputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))
outputs <- read.csv(file.path(root, 'database', 'outputs.csv'), stringsAsFactors = FALSE,
  na.strings = c('', 'NA'))

if (any(register$candidate_unit_id == unit_id, na.rm = TRUE) &&
    !any(register$candidate_id == candidate_id & register$candidate_unit_id == unit_id, na.rm = TRUE))
  stop('The NOAA black rockfish candidate unit is already assigned to another row.', call. = FALSE)
if (any(stocks$authority_stock_id == '166727', na.rm = TRUE) ||
    any(stocks$stock_id == unit_id, na.rm = TRUE) ||
    any(assessments$assessment_id == unit_id, na.rm = TRUE) ||
    any(inputs$assessment_id == unit_id, na.rm = TRUE) ||
    any(outputs$assessment_id == unit_id, na.rm = TRUE))
  stop('Screening must not duplicate or create canonical black rockfish records.', call. = FALSE)

source_url_text <- paste(source_urls, collapse = '; ')
old_reviews <- list(
  'NOAA:10402' = list(
    catch_age_evidence = 'Legacy 2007 northern Pacific Coast identifier; the 2023 current assessment suite does not map one-to-one to this row.',
    survey_age_evidence = 'Do not assign current regional survey inputs to this legacy identifier; see the 2023 assessment-suite review.',
    weights_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    maturity_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    mortality_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    timing_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    output_evidence = 'Current 2023 regional model outputs are registered once under NOAA:166727, not under this legacy ID.'),
  'NOAA:11302' = list(
    catch_age_evidence = 'Legacy 2007 southern Pacific Coast identifier; the 2023 current assessment suite does not map one-to-one to this row.',
    survey_age_evidence = 'Do not assign current regional survey inputs to this legacy identifier; see the 2023 assessment-suite review.',
    weights_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    maturity_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    mortality_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    timing_evidence = 'Not assigned to this historical row because its boundary does not map one-to-one to the 2023 assessment areas.',
    output_evidence = 'Current 2023 regional model outputs are registered once under NOAA:166727, not under this legacy ID.'))
for (id in names(old_reviews)) {
  row <- match(id, register$candidate_id)
  if (is.na(row) || sum(register$candidate_id == id) != 1L)
    stop('Expected one legacy black rockfish inventory row: ', id, call. = FALSE)
  if (register$scientific_name[row] != 'Sebastes melanops')
    stop('Unexpected scientific name for ', id, call. = FALSE)
  old_state <- register$match_status[row] == 'ambiguous' &&
    register$screening_status[row] == 'screened' &&
    register$eligibility[row] == 'ineligible' &&
    (is.na(register$candidate_unit_id[row]) || !nzchar(register$candidate_unit_id[row]))
  initial_state <- register$match_status[row] == 'unreviewed' &&
    register$screening_status[row] == 'not_started' &&
    register$eligibility[row] == 'unresolved' &&
    (is.na(register$candidate_unit_id[row]) || !nzchar(register$candidate_unit_id[row]))
  if (!old_state && !initial_state)
    stop('Legacy black rockfish row has a conflicting review state: ', id, call. = FALSE)
  review <- c(list(candidate_unit_id = NA_character_, match_status = 'ambiguous',
    screening_status = 'screened', eligibility = 'ineligible', import_status = 'not_started'),
    old_reviews[[id]], list(detailed_source_url = source_url_text, review_date = '2026-10-30',
      decision_reason = paste('This 2007 Stock SMART area identifier does not map one-to-one to the four regional model components in the current 2023 assessment suite. Keep it as an ambiguous historical identity, not as a separate current candidate. The suite is registered once under NOAA:166727.'),
      notes = paste('The old northern/southern labels are not a verified crosswalk to Washington, Oregon, northern California, or central California. See discovery/source_reviews/noaa_black_rockfish_pacific_coast_2023.md.')))
  for (name in names(review)) register[[name]][row] <- review[[name]]
}

current_row <- data.frame(
  candidate_id = candidate_id,
  candidate_unit_id = unit_id,
  organization = 'NOAA',
  official_stock_id = '166727',
  scientific_name = 'Sebastes melanops',
  stock_name = 'Black rockfish - U.S. Pacific Coast 2023 assessment suite',
  assessment_area = 'Washington, Oregon, northern California, central California',
  inventory_assessment_year = '2023',
  inventory_url = 'https://apps-st.fisheries.noaa.gov/stocksmart/?app=browse_by_stock&stockid=166727',
  retrieved_on = '2026-10-30',
  matched_stock_id = NA_character_,
  match_status = 'new',
  screening_status = 'screened',
  eligibility = 'eligible',
  import_status = 'not_started',
  catch_age_evidence = paste('The four Stock Synthesis model archives contain fleet catch data and fishery age/length compositions. Catch-at-age is not supplied as one direct matrix; reconstruct it from source catches, age/length compositions, and component-specific biological inputs.'),
  survey_age_evidence = paste('Regional abundance-index series and survey age/length compositions are present in the native model inputs. The inputs support documented reconstruction of age-indexed observations; do not treat the aggregate survey indices as direct age-specific series.'),
  weights_evidence = paste('Sex- and region-specific growth and weight relationships are represented in the native control/data files and generated age tables. Preserve the regional model differences.'),
  maturity_evidence = paste('The regional model files contain maturity and fecundity settings; California reports spawning output in eggs, not mature biomass.'),
  mortality_evidence = paste('Natural-mortality assumptions are specified by regional model and sex. Preserve those component-specific values and estimation settings.'),
  timing_evidence = paste('Fleet and survey timing settings are retained in each region-specific Stock Synthesis model configuration; extract timing separately by component.'),
  output_evidence = paste('The model archives contain sex-specific numeric N-at-age and biomass-at-age tables. Ages are 0-40 in Washington and Oregon and 0-50 in northern and central California; output histories begin in 1940, 1892, and 1875, respectively. Spawning output is in millions of eggs.'),
  detailed_source_url = source_url_text,
  review_date = '2026-10-30',
  decision_reason = paste('Eligible for source curation: all four current regional assessment components have native inputs and numeric age-specific population outputs. Reconstructing catch-at-age and age-indexed survey series from the source totals and compositions requires documented transformations. The regional components remain separate in source data; grouping them here counts one assessment case.'),
  notes = paste('NOAA Stock SMART identifier 166727 is shared across the 2023 regional assessment suite. The four components are Washington, Oregon, northern California, and central California; this grouping is an inventory-counting convention, not a claim of one biological population. Do not count the legacy 2007 NOAA:10402 and NOAA:11302 rows separately. Spawning output is egg-based, not SSB in tonnes. See discovery/source_reviews/noaa_black_rockfish_pacific_coast_2023.md and the local source manifest.'),
  stringsAsFactors = FALSE, check.names = FALSE)

current_idx <- which(register$candidate_id == candidate_id)
if (length(current_idx) > 1L)
  stop('Duplicate current black rockfish inventory ID: ', candidate_id, call. = FALSE)
if (length(current_idx) == 1L) {
  row <- current_idx
  if (register$official_stock_id[row] != '166727' ||
      register$scientific_name[row] != 'Sebastes melanops' ||
      register$candidate_unit_id[row] != unit_id ||
      !register$match_status[row] %in% c('unreviewed', 'new'))
    stop('Current NOAA black rockfish inventory row conflicts with the reviewed identity.', call. = FALSE)
  current_state <- register$match_status[row] == 'new' &&
    register$screening_status[row] == 'screened' && register$eligibility[row] == 'eligible'
  if (current_state) {
    if (!all(vapply(names(current_row), function(name)
      identical(as.character(register[[name]][row]), as.character(current_row[[name]])), logical(1L))))
      stop('Current NOAA black rockfish row differs from the reviewed screen.', call. = FALSE)
  } else {
    register[row, names(current_row)] <- current_row
  }
} else {
  register <- rbind(register, current_row)
}

csv_line <- function(row) {
  tmp <- tempfile(fileext = '.csv')
  write.table(row, tmp, sep = ',', row.names = FALSE, col.names = FALSE,
    quote = TRUE, na = '', qmethod = 'double')
  raw <- readBin(tmp, 'raw', n = file.info(tmp)$size)
  unlink(tmp)
  charToRaw(sub('\\r?\\n$', '', rawToChar(raw), perl = TRUE))
}
write_row <- function(path, index, row) {
  bytes <- readBin(path, 'raw', n = file.info(path)$size)
  line_feeds <- which(bytes == as.raw(10L))
  physical_line <- index + 1L
  if (length(line_feeds) < physical_line)
    stop('CSV row is out of range: ', index, call. = FALSE)
  first <- line_feeds[physical_line - 1L] + 1L
  last <- line_feeds[physical_line]
  ending <- if (last > first && bytes[last - 1L] == as.raw(13L))
    as.raw(c(13L, 10L)) else as.raw(10L)
  before <- bytes[seq_len(first - 1L)]
  after <- if (last < length(bytes)) bytes[(last + 1L):length(bytes)] else raw(0L)
  writeBin(c(before, csv_line(row), ending, after), path)
}
append_row <- function(path, row) {
  bytes <- readBin(path, 'raw', n = file.info(path)$size)
  line_feeds <- which(bytes == as.raw(10L))
  ending <- if (length(line_feeds) > 0L && line_feeds[1L] > 1L &&
      bytes[line_feeds[1L] - 1L] == as.raw(13L)) as.raw(c(13L, 10L)) else as.raw(10L)
  ends_with_newline <- length(bytes) > 0L && bytes[length(bytes)] == as.raw(10L)
  suffix <- if (ends_with_newline) raw(0L) else ending
  writeBin(c(bytes, suffix, csv_line(row), ending), path)
}

if (!identical(Sys.getenv('TAM_SOURCE_DRY_RUN'), '1')) {
  for (id in names(old_reviews)) {
    row <- match(id, register$candidate_id)
    if (!identical(as.character(read.csv(register_path, stringsAsFactors = FALSE,
        na.strings = c('', 'NA'))[[1L]][row]), id))
      stop('Register row ordering changed during screen.', call. = FALSE)
    existing <- read.csv(register_path, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
    if (!identical(as.character(existing$match_status[row]), 'ambiguous') ||
        !identical(as.character(existing$screening_status[row]), 'screened'))
      write_row(register_path, row, register[row, , drop = FALSE])
  }
  existing <- read.csv(register_path, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
  row <- match(candidate_id, existing$candidate_id)
  if (is.na(row)) append_row(register_path, current_row) else
    if (!identical(as.character(existing$eligibility[row]), 'eligible'))
      write_row(register_path, row, current_row)
}

cat('Screened NOAA Pacific Coast black rockfish as one eligible assessment case with four preserved regional components.\n')
