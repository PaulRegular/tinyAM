.assessment_database_files <- c(
  "stocks.csv", "assessments.csv", "assumptions.csv", "inputs.csv", "outputs.csv"
)

.assessment_database <- function(read_table, source, commit = NA_character_) {
  tables <- lapply(.assessment_database_files, read_table)
  names(tables) <- c("stocks", "assessments", "assumptions", "inputs", "outputs")
  tables$database_source <- source
  tables$commit <- commit
  tables
}

.database_git <- function(repo, args, stdout = TRUE) {
  system2("git", c("-C", shQuote(normalizePath(repo, winslash = "/", mustWork = TRUE)),
    args), stdout = stdout, stderr = TRUE)
}

.database_head <- function(repo = ".", revision = "HEAD") {
  value <- suppressWarnings(tryCatch(.database_git(repo,
    c("rev-parse", "--verify", paste0(revision, "^{commit}"))),
    error = function(e) character()))
  if (length(value) != 1L || !grepl("^[[:xdigit:]]{40}$", value)) return(NA_character_)
  value[[1]]
}

saio_path <- function() {
  configured <- getOption("tinyAM.saio_path", Sys.getenv("SAIO_PATH", ""))
  common <- suppressWarnings(tryCatch(system2("git",
    c("rev-parse", "--path-format=absolute", "--git-common-dir"),
    stdout = TRUE, stderr = FALSE), error = function(e) character()))
  candidates <- if (nzchar(configured)) configured else c(
    if (length(common) == 1L) file.path(dirname(dirname(common)), "saio"),
    file.path("..", "saio")
  )
  found <- candidates[dir.exists(file.path(candidates, "database"))]
  if (!length(found)) cli::cli_abort(c(
    "Cannot find the independent SAIO checkout.",
    "i" = "Set SAIO_PATH or {.code options(tinyAM.saio_path = '/path/to/saio')}."
  ))
  normalizePath(found[[1]], winslash = "/", mustWork = TRUE)
}

.saio_revision <- function(tier = "curated") {
  tier <- match.arg(tier, c("curated", "preliminary"))
  pins <- read.csv(file.path("analysis", "comp_assessments", "saio_revisions.csv"),
    stringsAsFactors = FALSE)
  revision <- pins$revision[pins$tier == tier]
  if (length(revision) != 1L || !grepl("^[[:xdigit:]]{40}$", revision)) {
    cli::cli_abort("saio_revisions.csv must contain one complete SAIO commit ID per tier.")
  }
  revision
}

.saio_read_table <- function(name, repo = saio_path(), revision = .saio_revision(tier),
                             tier = "curated") {
  pointer <- tempfile("saio_table_", fileext = ".csv")
  on.exit(unlink(pointer), add = TRUE)
  status <- .database_git(repo, c("show", paste0(revision, ":database/", name)),
    stdout = pointer)
  code <- if (is.numeric(status)) status else attr(status, "status")
  if ((!is.null(code) && code != 0L) || !file.exists(pointer)) {
    cli::cli_abort("Cannot read {.file {name}} at SAIO revision {.val {revision}}.")
  }
  first <- readLines(pointer, n = 1L, warn = FALSE)
  table_path <- pointer
  if (identical(first, "version https://git-lfs.github.com/spec/v1")) {
    oid <- sub("^oid sha256:", "", readLines(pointer, warn = FALSE)[[2]])
    if (!grepl("^[[:xdigit:]]{64}$", oid)) cli::cli_abort("Invalid SAIO LFS pointer.")
    common <- .database_git(repo,
      c("rev-parse", "--path-format=absolute", "--git-common-dir"))
    table_path <- file.path(common[[1]], "lfs", "objects", substr(oid, 1L, 2L),
      substr(oid, 3L, 4L), oid)
    if (!file.exists(table_path)) cli::cli_abort(c(
      "The pinned SAIO data object for {.file {name}} is unavailable locally.",
      "i" = "Fetch the pinned revision's Git LFS objects in the SAIO checkout."
    ))
  }
  read.csv(table_path, stringsAsFactors = FALSE,
    na.strings = c("", "NA"), check.names = FALSE)
}

read_database <- function(database_dir = NULL, tier = "curated") {
  if (is.null(database_dir)) return(read_committed_database(tier = tier))
  if (!dir.exists(database_dir)) cli::cli_abort("Database directory does not exist: {.file {database_dir}}.")
  .assessment_database(function(name) read.csv(file.path(database_dir, name),
    stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE),
    source = "working_tree", commit = .database_head(database_dir))
}

read_committed_database <- function(revision = .saio_revision(tier), repo = saio_path(),
                                    tier = "curated") {
  commit <- .database_head(repo, revision)
  if (is.na(commit)) cli::cli_abort("Cannot resolve the pinned SAIO revision {.val {revision}}.")
  .assessment_database(function(name) .saio_read_table(name, repo, commit),
    source = "saio_committed", commit = commit)
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

read_committed_assessment <- function(assessment_id, database = NULL) {
  if (is.null(database)) database <- read_committed_database()
  read_assessment(assessment_id, database)
}
