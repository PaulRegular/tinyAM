base <- "analysis/comp_assessments"
database <- file.path(base, "database")
stock_id <- "ck_nefsc_garmiii_ferruginea_sne_ma"
cache <- file.path(base, "source_cache", stock_id)
framework_id <- "nefsc_yellowtail_snema_2024_framework"
update_id <- "nefsc_yellowtail_snema_2025"
source_url <- "https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2024&species_id=13&stock_id=11&review_type_id=5"
update_url <- "https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2025&species_id=13&stock_id=5&review_type_id=2"
panel_url <- "https://d23h0vhsm26o6d.cloudfront.net/10b.-Panel-Summary-Report_Yellowtail-Flounder_508_final_12_12_2024.pdf"
adoption_url <- "https://d23h0vhsm26o6d.cloudfront.net/260423_Groundfish_FW72_SIR_FINAL-Submission.pdf"
read_pdf <- function(file, hash) {
  path <- file.path(cache, file)
  if (!file.exists(path) || digest::digest(readBin(path, "raw", n = file.info(path)$size),
    serialize = FALSE, algo = "sha256") != hash) cli::cli_abort("Cache and verify the reviewed source {file} before importing.")
  pdftools::pdf_text(path)
}
report <- read_pdf("Yellowtail_Flounder_RT_WG_Report_2024.pdf",
  "805ba9617e2231bb000be98ea7a0982c8943cc0a6661034c8c256bc987857878")
panel <- read_pdf("peer_review_2024.pdf",
  "be557aec7f5b8e5a823acbe6f8917a857ffd519a63b086a8b6ec5e995948e707")
update <- read_pdf("2025_YEL_SNEMA/2025_Yellowtail_Flounder_SNEMA_REPORT_UPDATED.pdf",
  "5421a283dd3c50bfbe04d83c28b89ed3bcf6d20cb1fc0a1fc79a933109201f71")
adoption <- read_pdf("NEFMC_FW72_2026.pdf",
  "d961314f6fd146694956a50e8f128277314fe1d103aaa8a24c1236ccd1c72c9c")
stopifnot(grepl("age specific M values", panel[22], fixed = TRUE),
  grepl("All other settings and", panel[22], fixed = TRUE),
  grepl("38 mt", adoption[27], fixed = TRUE), grepl("0.051", adoption[27], fixed = TRUE),
  grepl("removing the IID random effects", update[4], fixed = TRUE))

annual <- function(text, years, columns) {
  lines <- strsplit(paste(text, collapse = "\n"), "\n")[[1L]]
  lines <- trimws(lines[grepl("^\\s*(19|20)[0-9]{2}\\s", lines, perl = TRUE)])
  cells <- strsplit(gsub(" *\\*", "", lines), "[[:space:]]+", perl = TRUE)
  if (length(cells) != length(years) || any(lengths(cells) != columns + 1L)) cli::cli_abort("Unexpected source table dimensions.")
  x <- do.call(rbind, cells)
  stopifnot(identical(as.integer(x[, 1]), years))
  missing <- x %in% c("--", "\u2013")
  values <- suppressWarnings(as.numeric(x))
  stopifnot(all(is.finite(values[!missing])), all(is.na(values[missing])))
  matrix(values, nrow = nrow(x))
}
catch <- annual(report[88:89], 1973:2022, 6L)
spring_years <- setdiff(1968:2022, 2020L)
spring <- cbind(annual(report[187:188], spring_years, 6L),
  annual(report[189:190], spring_years, 6L)[, -1, drop = FALSE])
fall_years <- setdiff(1963:2022, c(2017L, 2020L))
fall <- annual(report[191:192], fall_years, 10L)
winter <- annual(report[193], 1992:2007, 8L)
aggregate <- annual(report[185:186], 1963:2022, 8L)
spring_difference <- round(rowSums(spring[, -1]) - aggregate[match(spring_years, aggregate[, 1]), 4], 6)
fall_difference <- round(rowSums(fall[, -1]) - aggregate[match(fall_years, aggregate[, 1]), 2], 6)
stopifnot(all(is.na(catch[catch[, 1] >= 2020, -1])),
  all(is.finite(catch[catch[, 1] < 2020, -1])),
  identical(spring_years[abs(spring_difference) > .0065], c(2012L, 2015L, 2017L)),
  identical(spring_difference[spring_years %in% c(2012L, 2015L, 2017L)], c(-.016, -.007, -.012)),
  identical(fall_years[abs(fall_difference) > .0055], 2013L),
  fall_difference[fall_years == 2013] == -.008)

blank <- function(file, id, n = 1L) {
  fields <- names(read.csv(file.path(database, file), nrows = 0, check.names = FALSE))
  x <- as.data.frame(matrix("", n, length(fields)), stringsAsFactors = FALSE)
  names(x) <- fields; x$assessment_id <- id
  x
}
number <- function(x) ifelse(is.na(x), "", format(x, digits = 15, trim = TRUE, scientific = FALSE))
input <- function(id, type, measure, basis, year, age = NA, value, unit,
                  fleet = "", survey = "", season = "", reference, notes) {
  n <- length(value); x <- blank("inputs.csv", id, n)
  x$type <- type; x$measure <- measure; x$basis <- basis
  x$year <- year; x$year_basis <- "calendar_year"; x$age <- number(rep_len(age, n))
  x$value <- number(value); x$unit <- unit; x$fleet <- fleet
  x$survey <- survey; x$season <- season; x$source_type <- "official_table"
  x$source_reference <- reference; x$notes <- notes
  x
}
age_input <- function(table, ages, type, unit, reference, notes, survey = "", season = "", fleet = "") {
  g <- expand.grid(year = table[, 1], age = ages, KEEP.OUT.ATTRS = FALSE)
  input(framework_id, type, "numbers_at_age", "numbers", g$year, g$age,
    as.vector(table[, -1]), unit, fleet, survey, season, reference, notes)
}
inputs <- age_input(catch, 1:6, "catch", "thousand fish",
  paste(source_url, "Yellowtail Flounder RT WG Report, Table 2.8.3, pp. 88-89; peer review pp. 17, 22"),
  "Combined landings and discards; age 6 is 6+. Source numbers support the age composition, not a separate independent catch-at-age likelihood. The source fits an aggregate fleet total plus composition. Printed dashes in 2020-2022 remain missing; no ages are inferred from total catches.", fleet = "Combined commercial landings and dead discards")
for (survey in c("spring", "fall", "winter")) {
  table <- get(survey); ages <- switch(survey, spring = 0:11, fall = 1:10, winter = 1:8)
  table_id <- switch(survey, spring = "3.6.2a-b, pp. 187-190", fall = "3.6.3, pp. 191-192", winter = "3.6.4, p. 193")
  x <- age_input(table, ages, "index", "fish per tow",
    paste(source_url, "Yellowtail Flounder RT WG Report, Tables", table_id, "; survey description pp. 142-148"),
    "Original stratified-mean numbers-at-age, retaining finer ages and pre-1973 survey years from the source table. These provide age composition alongside an aggregate index; they are not all independent absolute observations in WHAM. Exact fitted-cell selection, aggregate-index basis and model timing fractions remain unresolved. No plus-group pooling, zero replacement or fitted-q correction is applied.",
    survey = paste("NEFSC", survey, "bottom trawl survey"), season = survey)
  if (survey == "fall") {
    flagged <- x$year == 2013 & x$age %in% c(2:7, 9)
    x$notes[flagged] <- paste(x$notes[flagged], "Source asterisk flags an incomplete age-length key; retain the reported value.")
  }
  discrepant <- if (survey == "spring") c(2012L, 2015L, 2017L) else if (survey == "fall") 2013L else integer()
  x$notes[x$year %in% discrepant] <- paste(x$notes[x$year %in% discrepant],
    "The sum of printed age values differs from the Table 3.6.1 aggregate beyond the three-decimal rounding bound. Cause unresolved; neither representation is adjusted or rescaled.")
  inputs <- rbind(inputs, x)
}
for (season in c("fall", "spring")) {
  offset <- if (season == "fall") 0L else 2L
  survey <- paste("NEFSC", season, "bottom trawl survey")
  for (biomass in c(FALSE, TRUE)) {
    column <- 2L + offset + if (biomass) 4L else 0L
    kind <- if (biomass) "biomass" else "abundance"
    reference <- paste(source_url, "Yellowtail Flounder RT WG Report, Table 3.6.1, pp. 185-186;", season, kind, "and CV columns")
    inputs <- rbind(inputs,
      input(framework_id, "index", if (biomass) "total_biomass" else "total_numbers",
        if (biomass) "biomass" else "numbers", aggregate[, 1], value = aggregate[, column],
        unit = if (biomass) "kg per tow" else "fish per tow", survey = survey, season = season,
        reference = reference,
        notes = "Published aggregate stratified mean; number and biomass representations are not independent duplicate likelihood series. Exact production likelihood basis and index-CV treatment remain unresolved. Printed missing values stay missing. The source prints spring 2020 as zero without a CV after incomplete coverage; it is retained as printed, not certified as a fitted zero observation."),
      input(framework_id, "index", if (biomass) "biomass_index_cv" else "abundance_index_cv", "relative_scale", aggregate[, 1],
        value = aggregate[, column + 1L], unit = "ratio",
        survey = survey, season = season, reference = reference,
        notes = "Source CV ratio for the corresponding aggregate index, not an SE or age-composition weight. No conversion to a fixed likelihood SD is made; missing CVs remain missing."))
  }
}
update_lines <- strsplit(update[2], "\n")[[1L]]
series <- function(label) {
  line <- trimws(update_lines[grepl(paste0("^\\s*", label, "\\s+[0-9]"), update_lines, perl = TRUE)])
  stopifnot(length(line) == 1L)
  cells <- strsplit(trimws(sub(label, "", line, fixed = TRUE)), "[[:space:]]+")[[1L]]
  values <- as.numeric(gsub(",", "", cells, fixed = TRUE))
  stopifnot(length(values) == 10L, all(is.finite(values)))
  values
}
inputs <- rbind(inputs, input(update_id, "catch", "total_biomass", "biomass", 2015:2024,
  value = series("Total Catch"), unit = "tonnes", fleet = "Combined commercial landings and dead discards",
  reference = paste(update_url, "2025 updated assessment report, Table 1, p. 1"),
  notes = "Source combined fleet catch total, including landings and dead discards. Values below 0.5 tonnes print as zero in the component columns; no unrounded component totals are reconstructed. This update revises catches and biology; no 2024 framework ages are attached to it."))
output <- function(id, measure, value, unit, year = NA, age = NA, group = "", reference, notes) {
  n <- length(value); x <- blank("outputs.csv", id, n)
  x$type <- switch(measure, SSB = "biomass", recruitment = "recruitment", "mortality")
  x$measure <- measure; x$year <- number(rep_len(year, n)); x$age <- number(rep_len(age, n))
  x$age_group <- group; x$value <- number(value); x$unit <- unit; x$source_type <- "official_table"
  x$source_reference <- reference; x$notes <- notes
  x
}
outputs <- rbind(
  output(update_id, "SSB", series("Spawning Stock Biomass"), "tonnes", 2015:2024,
    reference = paste(update_url, "2025 updated assessment report, Table 1, p. 1;", adoption_url, "p. 27"),
    notes = "Historical native WHAM SSB. The public report retains its draft watermark; the final Framework 72 submission confirms the applied update and terminal SSB of 38 tonnes. No numerical historical intervals or SEs are supplied in this table; none are inferred from plotted ribbons."),
  output(update_id, "Ffull", series("FF ull"), "per year", 2015:2024, group = "5-6+",
    reference = paste(update_url, "2025 updated assessment report, Table 1, p. 1;", adoption_url, "p. 27"),
    notes = "Fully selected fishing mortality for ages 5-6+, as labelled in Table 1; not population-weighted Fbar or a numerical F-at-age surface. Final Framework 72 confirms terminal F of 0.051. Uncertainty is not supplied numerically."),
  output(update_id, "recruitment", series("Recruits"), "thousand fish", 2015:2024, age = 1,
    reference = paste(update_url, "2025 updated assessment report, Tables 1-2, p. 1"),
    notes = "Age-1 WHAM recruitment from the 2025 update. No numerical uncertainty supplied; no shifting to another cohort year or mixing with the 2024 constant-M candidate."))
for (i in 1:3) {
  x <- output(framework_id, c("Fmsy", "Bmsy", "msy")[i],
    c(.36, 191, 70)[i], c("per year", "tonnes", "tonnes per year")[i],
    reference = paste(panel_url, "Table 5.5, p. 24, age-based-M row"),
    notes = "Selected revised framework F40% proxy, corresponding SSBMSY proxy or equilibrium MSY proxy at F40%, respectively; not a historical annual estimate or an estimated true MSY reference point. Parenthetical bounds are retained as reported; their interval level and construction are not stated in this table. No SE inferred. Superseded constant-M row excluded.")
  x$type <- c("mortality", "biomass", "catch")[i]
  x$lwr <- c(.30, 11, 4)[i]; x$upr <- c(.43, 3345, 1157)[i]
  outputs <- rbind(outputs, x)
}
assumption <- function(id, component, setting, value, reference, notes = "") {
  x <- blank("assumptions.csv", id); x$component <- component; x$setting <- setting
  x$value <- value; x$source_reference <- reference; x$notes <- notes; x
}
assumptions <- rbind(
  assumption(framework_id, "assessment", "selected_run", "Revised Candidate Model with age-specific M; not the original constant-M candidate", paste(panel_url, "pp. 17-20, 22-24"), "Final panel selected the age-based-M revision, with other settings and input unchanged for that comparison. Original candidate point estimates are excluded."),
  assumption(framework_id, "population", "dimensions", "1973-2022; ages 1-6+; one combined commercial fleet; three NEFSC trawl surveys; combined sexes", paste(panel_url, "p. 17")),
  assumption(framework_id, "M", "model", "Fixed age-specific M; higher at young ages, scaled from lifetime mortality using length-at-age", paste(panel_url, "pp. 22-24; Cadrin working paper pp. 2-3"), "Exact supplied production vector is unrecovered. The superseded constant M=0.5 is not imported as the revised model's M surface."),
  assumption(framework_id, "N", "process", "2D AR1 older-cohort N process; recruitment decoupled from ages 2-6+", paste(source_url, "WG report pp. 294-295;", panel_url, "p. 17")),
  assumption(framework_id, "recruitment", "process", "Random about a mean, with AR1 temporal recruitment and a one-year-lagged GSI controlling effect; GSI modelled as AR1", paste(source_url, "WG report pp. 294-295;", panel_url, "p. 17"), "No Beverton-Holt curve in the selected SNEMA model. Exact native covariate observations, scaling, coupling and lognormal corrections need production configuration."),
  assumption(framework_id, "F", "selectivity", "Logistic fleet selectivity with IID random effects", paste(panel_url, "p. 17")),
  assumption(framework_id, "catch", "likelihood", "Combined total landings and dead discards, CV 0.05; logistic-normal AR1 miss0 age composition", paste(source_url, "WG report pp. 290-295;", panel_url, "p. 17"), "Zero composition cells are treated as missing within the composition likelihood, not recoded source observations. Exact likelihood weights and selected cells are unresolved."),
  assumption(framework_id, "index", "included_surveys", "NEFSC spring, fall and winter bottom trawl surveys; larval indices excluded", paste(panel_url, "p. 17")),
  assumption(framework_id, "q", "selectivity", "Spring and winter logistic; fall age-specific; no q or survey-selectivity random effects", paste(panel_url, "p. 17")),
  assumption(framework_id, "index", "likelihood", "Aggregate indices plus logistic-normal AR1 miss0 age compositions", paste(source_url, "WG report pp. 292-295;", panel_url, "p. 17"), "Published finer age-specific means are retained. Native aggregate likelihood basis, variance weights and exact fitted age cells are unresolved."),
  assumption(framework_id, "index", "sampling_time", "unknown", paste(source_url, "WG report pp. 142-148; Hodgdon working paper section 3.2"), "Seasonal labels are known, but exact model timing fractions were not recovered. No generic seasonal midpoint is inserted."),
  assumption(framework_id, "index", "published_total_discrepancies", "Age sums are below number indices by 0.016, 0.007 and 0.012 fish/tow in spring 2012, 2015 and 2017; by 0.008 in fall 2013", paste(source_url, "Tables 3.6.1-3.6.3, pp. 185-192"), "Beyond the bound from rounding each printed cell to three decimals. Both source representations are preserved; no explanation or correction is invented. Spring 2020 prints zero without a CV but has no age table row after incomplete coverage."),
  assumption(framework_id, "biology", "coverage_gap", "Annual model weights and exact fixed maturity vector not recovered", paste(source_url, "WG report pp. 307-309; Hodgdon working paper section 6"), "Plots and lifetime growth approximations do not recover annual model biology or the exact age-6+ weights."),
  assumption(framework_id, "population", "initial_state_and_bias", "unknown", paste(source_url, "WG report pp. 290-295"), "Native initial-state activation, priors, bias corrections and pinned WHAM version remain unresolved."),
  assumption(update_id, "assessment", "production_verification", "2025 applied WHAM management-track update through 2024", paste(adoption_url, "p. 27;", update_url, "updated report Table 1"), "Report carries a draft watermark. Final Council submission corroborates use of this update, SSB 38 t, F 0.051 and reference points; no native fitted object recovered."),
  assumption(update_id, "F", "selectivity_change", "Fleet IID selectivity random effects removed to address convergence", paste(update_url, "updated report, pp. 2-3"), "Do not carry the 2024 random-selectivity configuration into this update."),
  assumption(update_id, "population", "bias_adjustment", "WHAM lognormal adjustment switched off", paste(update_url, "updated report, p. 3"), "Exact option mapping needs the production files; this statement does not recover every likelihood bias convention."),
  assumption(update_id, "biology", "revised_inputs", "CAMS discards, maturity and GSI revised; no catch age composition in the terminal five years", paste(update_url, "updated report, pp. 2-3"), "2024 framework arrays must not be attached to this update. Complete revised observations and biology remain unavailable."),
  assumption(update_id, "M", "model", "Age-specific fixed M retained from revised research-track model; numerical vector unrecovered", paste(panel_url, "pp. 22-24;", update_url, "pp. 1-3")),
  assumption(update_id, "recruitment", "GSI", "Recruitment strongly informed by Gulf Stream Index", paste(update_url, "updated report, pp. 2-3"), "Future GSI assumptions and projection outputs are separate from fitted historical recruitment."))
assessment <- function(id, year, terminal, estimate_terminal, current, applied, notes) {
  x <- blank("assessments.csv", id); x$stock_id <- stock_id
  x$assessment_year <- year; x$terminal_year <- terminal; x$estimate_terminal_year <- estimate_terminal
  x$assessment_type <- if (current) "research_track_framework" else "management_track_summary"
  x$model_family <- "WHAM"; x$model_version <- "Exact production revision unrecovered"
  x$is_current <- toupper(as.character(current)); x$is_applied <- toupper(as.character(applied))
  x$framework_year <- 2024; x$assessment_url <- if (current) source_url else update_url
  x$framework_url <- panel_url; x$data_url <- if (current) source_url else ""
  x$assumptions_status <- x$inputs_status <- x$outputs_status <- "partial"; x$notes <- notes; x
}
assessments <- rbind(
  assessment(framework_id, 2024, 2022, 2022, TRUE, FALSE,
    "Latest recoverable detailed approved framework: revised age-specific-M model selected in November 2024 for subsequent management-track use. Source catch-at-age, finer trawl indices and aggregates recovered; exact selected M/biology/covariate arrays, likelihood settings and numerical fitted age surfaces missing. Original constant-M candidate outputs are not official revised-run values. This framework was not itself the 2024 management-track advice run."),
  assessment(update_id, 2025, 2024, 2028, FALSE, TRUE,
    "Applied 2025 update corroborated by final NEFMC Framework 72. Public updated report retains draft watermark; canonical values retain that provenance. Historical catch totals, SSB, age-1 R and fully selected F recovered. Native inputs, biology and fitted N/F-at-age unavailable; do not mix with 2024 framework. Projection estimates through 2028 are not imported as historical observations or population estimates."))
stocks <- read.csv(file.path(database, "stocks.csv"), colClasses = "character", na.strings = "", check.names = FALSE)
at <- which(stocks$stock_id == stock_id)
stopifnot(length(at) == 1L, stocks$charbonneau_id[at] == "NEFSC-GARMIII_ferruginea_SNE, MA")
stocks$authority_stock_id[at] <- "Southern New England/Mid-Atlantic yellowtail flounder"
stocks$common_name[at] <- "Southern New England/Mid-Atlantic yellowtail flounder"
stocks$area[at] <- "Southern New England/Mid-Atlantic"
stocks$region[at] <- "Southern New England and Mid-Atlantic continental shelf"
stocks$notes[at] <- "Charbonneau-Keith NEFSC-GARMIII_ferruginea_SNE, MA; original PDF label SNE MA Limada ferrugina.pdf retained here. The selected 2024 WHAM framework uses age-specific M, superseding the original constant-M candidate; the applied 2025 update has revised inputs and selectivity settings. Partial official source tables are imported; full production biology, covariate/configuration arrays and fitted age surfaces remain unrecovered."

write_rows <- function(file, x, id) {
  path <- file.path(database, file)
  contents <- rawToChar(readBin(path, "raw", n = file.info(path)$size))
  lines <- strsplit(contents, "\n", fixed = TRUE)[[1L]]
  selected <- which(startsWith(lines, paste0('"', id, '",')))
  stopifnot(identical(names(x), names(read.csv(path, nrows = 0))))
  tmp <- tempfile(tmpdir = database); on.exit(unlink(tmp), add = TRUE)
  write.table(x, tmp, sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE,
    na = "", fileEncoding = "UTF-8")
  replacement <- readLines(tmp, encoding = "UTF-8", warn = FALSE)
  if (length(selected)) {
    if (!identical(selected, seq.int(min(selected), max(selected)))) cli::cli_abort("Existing {file} rows are not contiguous; review before replacing.")
    replacement <- paste0(replacement, if (endsWith(lines[selected[1L]], "\r")) "\r" else "")
    lines <- c(head(lines, min(selected) - 1L), replacement, tail(lines, length(lines) - max(selected)))
  } else lines <- c(lines, replacement)
  bytes <- charToRaw(paste0(paste(lines, collapse = "\n"), if (endsWith(contents, "\n")) "\n" else ""))
  if (identical(bytes, charToRaw(contents))) return(invisible(NULL))
  writeBin(bytes, tmp)
  if (!file.copy(tmp, path, overwrite = TRUE)) cli::cli_abort("Could not write {file}.")
}
write_rows("stocks.csv", stocks[at, ], stock_id)
for (id in c(framework_id, update_id)) {
  for (file in c("assessments", "assumptions", "inputs", "outputs")) {
    rows <- get(file)
    write_rows(paste0(file, ".csv"), rows[rows$assessment_id == id, ], id)
  }
}
cat("Imported separate SNEMA yellowtail framework/update:", nrow(inputs), "inputs,",
  nrow(assumptions), "assumptions,", nrow(outputs), "outputs; all remain partial.\n")
