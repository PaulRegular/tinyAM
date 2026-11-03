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

.assessment_cache <- function(assessment_id, tier = NULL) {
  if (is.null(tier)) tier <- .assessment_translation(assessment_id)$tier[[1]]
  file.path(.assessment_results(tier), "cache", assessment_id)
}

.assessment_cache_root <- function(tier = "curated") {
  file.path(.assessment_results(tier), "cache")
}

index_assessment_dashboards <- function(tier = "curated") {
  folder <- .assessment_results(tier)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(folder, "[.]html$", recursive = TRUE, full.names = TRUE)
  relative <- substring(gsub("\\", "/", files, fixed = TRUE), nchar(gsub("\\", "/", folder, fixed = TRUE)) + 2L)
  keep <- (basename(files) == "dashboard.html" | grepl("(^|/)dashboards/[^/]+[.]html$", relative)) &
    !grepl("/(check|source|tinyAM[.]Rcheck)/", relative)
  relative <- sort(relative[keep])
  escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  links <- vapply(relative, function(file) paste0('<li><a href="',
    escape(utils::URLencode(file)), '">', escape(file), '</a>',
    if (grepl("legacy", file)) " — legacy; revalidation required" else "", '</li>'), character(1))
  output <- file.path(folder, "index.html")
  writeLines(c('<!doctype html><html><head><meta charset="utf-8">',
    paste0('<title>', tier, ' assessment dashboards</title></head><body>'),
    paste0('<h1>', tier, ' assessment dashboards</h1>'),
    '<p>Migrated outputs retain their original estimates. Older outputs require review against the current script and data revision.</p>',
    '<ul>', links, '</ul></body></html>'), output, useBytes = TRUE)
  invisible(output)
}
