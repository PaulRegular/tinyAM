.database_head <- function(repo = ".") {
  value <- suppressWarnings(tryCatch(system2("git", c("-C", shQuote(repo),
    "rev-parse", "--verify", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()))
  if (length(value) != 1L || !grepl("^[[:xdigit:]]{40}$", value)) return(NA_character_)
  value[[1]]
}

.assessment_snapshot <- function(manifest = file.path("analysis", "comp_assessments",
                                  "saio_revisions.csv")) {
  read.csv(manifest, stringsAsFactors = FALSE)$revision[[1]]
}

# Read the manually refreshed snapshot. This function needs neither SAIO nor Git.
read_database <- function(cache_dir = file.path("analysis", "comp_assessments", "cache", "saio"),
                          manifest = file.path("analysis", "comp_assessments", "saio_revisions.csv")) {
  if (!file.exists(manifest)) cli::cli_abort("Run analysis/comp_assessments/scripts/update_database_cache.R before fitting: the SAIO snapshot manifest is missing.")
  revision <- .assessment_snapshot(manifest)
  folder <- file.path(cache_dir, revision)
  cached <- file.path(folder, "database.rds")
  if (!file.exists(cached)) cli::cli_abort(
    "The local SAIO snapshot is unavailable or incomplete. Run analysis/comp_assessments/scripts/update_database_cache.R before fitting.")
  tables <- readRDS(cached)
  tables$database_source <- "saio_snapshot"
  tables$commit <- tables$database_ref <- tables$database_key <- revision
  tables$database_selection <- "cached_snapshot"
  tables
}

read_assessment <- function(assessment_id, database = NULL) {
  if (length(assessment_id) != 1L || is.na(assessment_id) ||
    !nzchar(assessment_id)) {
    stop("assessment_id must be one non-empty value.", call. = FALSE)
  }
  if (is.null(database)) database <- read_database()
  required <- c("stocks", "assessments", "assumptions", "inputs", "outputs")
  if (!all(required %in% names(database))) {
    stop("database must contain stocks, assessments, assumptions, inputs, and outputs tables.",
      call. = FALSE
    )
  }
  assessment <- database$assessments[
    !is.na(database$assessments$assessment_id) &
      database$assessments$assessment_id == assessment_id, ,
    drop = FALSE
  ]
  if (nrow(assessment) != 1L) {
    stop("Database must contain exactly one assessment_id: ", assessment_id,
      call. = FALSE
    )
  }
  tables <- database
  stock_id <- assessment$stock_id[[1]]
  tables$stocks <- tables$stocks[
    tables$stocks$stock_id == stock_id, ,
    drop = FALSE
  ]
  for (table in c("assumptions", "inputs", "outputs")) {
    rows <- tables[[table]]$assessment_id == assessment_id
    rows[is.na(rows)] <- FALSE
    tables[[table]] <- tables[[table]][rows, , drop = FALSE]
  }
  tables$assessment <- assessment
  tables
}
