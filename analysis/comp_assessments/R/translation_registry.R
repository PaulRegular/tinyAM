assessment_translations <- function(tier = NULL) {
  registry <- read.csv(file.path("analysis", "comp_assessments", "translations.csv"),
    stringsAsFactors = FALSE)
  if (anyNA(registry[c("assessment_id", "tier")]) ||
      any(!registry$tier %in% c("curated", "preliminary")) ||
      anyDuplicated(paste(registry$assessment_id, registry$tier))) {
    cli::cli_abort("The translation registry needs unique assessment IDs within each tier.")
  }
  if (!is.null(tier)) {
    tier <- match.arg(tier, c("curated", "preliminary"))
    registry <- registry[registry$tier == tier, , drop = FALSE]
  }
  registry$script <- file.path("analysis", "comp_assessments", "scripts",
    "translation", registry$tier, paste0(registry$assessment_id, ".R"))
  registry
}

.assessment_translation <- function(assessment_id, tier = NULL, required = TRUE) {
  entries <- assessment_translations(tier)
  entries <- entries[entries$assessment_id == assessment_id, , drop = FALSE]
  if (is.null(tier) && any(entries$tier == "curated")) {
    entries <- entries[entries$tier == "curated", , drop = FALSE]
  }
  if (nrow(entries) != 1L) {
    if (!required && !nrow(entries)) return(NULL)
    cli::cli_abort("No unique {if (is.null(tier)) '' else tier} translation is registered for {.val {assessment_id}}.")
  }
  entries
}

.assessment_script <- function(assessment_id, tier = NULL) {
  .assessment_translation(assessment_id, tier)$script[[1]]
}

.assessment_results <- function(tier = "curated") {
  tier <- match.arg(tier, c("curated", "preliminary"))
  file.path(if (exists(".assessment_root", inherits = TRUE)) .assessment_root else
    file.path("analysis", "comp_assessments"), "results", tier)
}

.assessment_cache_key <- function(database = NULL) {
  if (is.null(database)) return(.saio_revision())
  key <- database$database_key
  if (is.null(key) && !identical(database$database_source, "working_tree")) key <- database$commit
  if (!is.character(key) || length(key) != 1L || is.na(key) ||
      !grepl("^([[:xdigit:]]{40}|working-[[:xdigit:]]{32})$", key)) {
    cli::cli_abort("A cache requires a recorded SAIO revision or working-table fingerprint; reread the source database first.")
  }
  key
}

.assessment_cache <- function(assessment_id, tier = NULL, database = NULL) {
  if (is.null(tier)) tier <- .assessment_translation(assessment_id)$tier[[1]]
  file.path(.assessment_cache_root(tier, database), assessment_id)
}

.assessment_cache_root <- function(tier = "curated", database = NULL) {
  file.path(.assessment_results(tier), "cache", .assessment_cache_key(database))
}

.assessment_cached_fit <- function(assessment_id, tier = NULL, database = NULL) {
  file <- file.path(.assessment_cache(assessment_id, tier, database), "fit.rds")
  if (!file.exists(file)) return(NULL)
  fit <- readRDS(file)
  if (!identical(fit$assessment_provenance$database_key, .assessment_cache_key(database))) {
    cli::cli_abort("Cached fit {.file {file}} lacks matching source provenance; regenerate it from the selected data revision.")
  }
  fit
}

.assessment_result_key <- function(result) {
  key <- result$diagnostics$database_key
  if (is.null(key)) key <- result$diagnostics$database_revision
  if (is.null(key)) return(.assessment_cache_key(result$source))
  .assessment_cache_key(list(database_key = key))
}

.assessment_cached_runs <- function(tier = "curated", database = NULL) {
  file <- file.path(.assessment_cache_root(tier, database), "assessment_runs.rds")
  if (!file.exists(file)) return(NULL)
  runs <- readRDS(file)
  if (length(runs) && any(vapply(runs, .assessment_result_key, character(1)) != .assessment_cache_key(database))) {
    cli::cli_abort("Cached assessment runs do not match the selected SAIO revision; retain them as historical results and run the current batch separately.")
  }
  runs
}

.assessment_check_results_revision <- function(folder, database = NULL) {
  file <- file.path(folder, "fit_diagnostics.csv")
  if (!file.exists(file)) cli::cli_abort("Run the selected assessment batch before plotting; its source-revision diagnostics are unavailable.")
  diagnostics <- read.csv(file, stringsAsFactors = FALSE)
  keys <- if ("database_key" %in% names(diagnostics)) diagnostics$database_key else diagnostics$database_revision
  if (!length(keys) || anyNA(keys) || any(keys != .assessment_cache_key(database))) {
    cli::cli_abort("These result tables use another or unknown SAIO revision. Regenerate the selected batch before plotting current results, or explicitly select its historical revision.")
  }
  invisible(TRUE)
}

index_assessment_dashboards <- function(tier = "curated", database = NULL) {
  selected <- .assessment_cache_key(database)
  folder <- .assessment_results(tier)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(folder, "[.]html$", recursive = TRUE, full.names = TRUE)
  relative <- substring(gsub("\\", "/", files, fixed = TRUE), nchar(gsub("\\", "/", folder, fixed = TRUE)) + 2L)
  keep <- (basename(files) == "dashboard.html" | grepl("(^|/)dashboards/[^/]+[.]html$", relative)) &
    !grepl("/(check|source|tinyAM[.]Rcheck|experiments)/", relative)
  relative <- sort(relative[keep])
  escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  links <- vapply(relative, function(file) {
    key <- if (grepl("^cache/([[:xdigit:]]{40}|working-[[:xdigit:]]{32})/", file)) {
      strsplit(file, "/", fixed = TRUE)[[1]][[2]]
    } else NA_character_
    status <- if (is.na(key)) "historical output; source revision unknown" else
      if (identical(key, selected)) paste0("selected source ", key) else
        paste0("historical source ", key)
    paste0('<li><a href="', escape(utils::URLencode(file)), '">',
      escape(file), '</a> — ', escape(status), '</li>')
  }, character(1))
  output <- file.path(folder, "index.html")
  writeLines(c('<!doctype html><html><head><meta charset="utf-8">',
    paste0('<title>', tier, ' assessment dashboards</title></head><body>'),
    paste0('<h1>', tier, ' assessment dashboards</h1>'),
    paste0('<p>Selected SAIO source: ', escape(selected),
      '. Source labels identify cached data revisions, not current-script validation.</p>'),
    '<p>Historical dashboards retain their original estimates and HTML. They are not results from the selected source; review them against their original data and model revisions.</p>',
    '<ul>', links, '</ul></body></html>'), output, useBytes = TRUE)
  invisible(output)
}
