.assessment_database_files <- c(
  "stocks.csv", "assessments.csv", "assumptions.csv", "inputs.csv", "outputs.csv"
)

.assessment_database <- function(read_table, source, commit = NA_character_,
                                 reference = NA_character_, selection = source,
                                 key = commit) {
  tables <- lapply(.assessment_database_files, read_table)
  names(tables) <- c("stocks", "assessments", "assumptions", "inputs", "outputs")
  tables$database_source <- source
  tables$commit <- commit
  tables$database_ref <- reference
  tables$database_selection <- selection
  tables$database_key <- key
  tables
}

.database_git <- function(repo, args, stdout = TRUE) {
  system2("git", c("-C", shQuote(normalizePath(repo, winslash = "/", mustWork = TRUE)),
    shQuote(args)), stdout = stdout, stderr = TRUE)
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

.saio_snapshot_revision <- function(tier = "curated") {
  tier <- match.arg(tier, c("curated", "preliminary"))
  pins <- read.csv(file.path("analysis", "comp_assessments", "saio_revisions.csv"),
    stringsAsFactors = FALSE)
  revision <- pins$revision[pins$tier == tier]
  if (length(revision) != 1L || is.na(revision) || !grepl("^[[:xdigit:]]{40}$", revision)) {
    cli::cli_abort("saio_revisions.csv must contain one complete SAIO commit ID per tier.")
  }
  revision
}

.saio_revision <- function(tier = "curated", repo = saio_path()) {
  match.arg(tier, c("curated", "preliminary"))
  main <- .database_head(repo, "refs/heads/main")
  if (is.na(main)) cli::cli_abort(c(
    "The local SAIO main branch is unavailable.",
    "i" = "Create or update main in the SAIO checkout before reading current data, or supply an explicit revision."
  ))
  remote <- .database_head(repo, "refs/remotes/origin/main")
  if (!is.na(remote)) {
    ancestor <- suppressWarnings(.database_git(repo,
      c("merge-base", "--is-ancestor", remote, main)))
    code <- attr(ancestor, "status")
    if (!is.null(code) && code != 0L) cli::cli_abort(c(
      "SAIO main is behind or diverged from the fetched origin/main.",
      "i" = "Update the SAIO checkout with git pull --ff-only, resolving any divergence explicitly, before reading current data."
    ))
  }
  main
}

.saio_metadata <- function(revision = NULL, repo = saio_path(), tier = "curated") {
  match.arg(tier, c("curated", "preliminary"))
  latest <- is.null(revision)
  if (latest) revision <- .saio_revision(tier, repo)
  if (!is.character(revision) || length(revision) != 1L || is.na(revision) || !nzchar(revision)) {
    cli::cli_abort("SAIO revision must identify one committed revision.")
  }
  commit <- .database_head(repo, revision)
  if (is.na(commit)) cli::cli_abort("Cannot resolve SAIO revision {.val {revision}}.")
  list(database_source = "saio_committed", commit = commit,
    database_ref = if (latest) "refs/heads/main" else revision,
    database_selection = if (latest) "latest_main" else "explicit_revision",
    database_key = commit)
}

.saio_read_table <- function(name, repo = saio_path(), revision = .saio_revision(tier, repo),
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
      "The committed SAIO data object for {.file {name}} is unavailable locally.",
      "i" = "Fetch this revision's Git LFS objects in the SAIO checkout before reading it."
    ))
  }
  read.csv(table_path, stringsAsFactors = FALSE,
    na.strings = c("", "NA"), check.names = FALSE)
}

read_database <- function(database_dir = NULL, tier = "curated", revision = NULL,
                          repo = saio_path()) {
  if (is.null(database_dir)) return(read_committed_database(revision, repo, tier))
  if (!is.null(revision)) cli::cli_abort("Supply a committed revision or a working database directory, not both.")
  if (!dir.exists(database_dir)) cli::cli_abort("Database directory does not exist: {.file {database_dir}}.")
  files <- file.path(database_dir, .assessment_database_files)
  if (any(!file.exists(files))) cli::cli_abort("The working database requires all five source tables.")
  hashes <- tools::md5sum(files)
  if (anyNA(hashes)) cli::cli_abort("The working database requires all five source tables.")
  fingerprint <- tempfile("saio_working_")
  on.exit(unlink(fingerprint), add = TRUE)
  writeLines(paste(.assessment_database_files, unname(hashes)), fingerprint)
  key <- paste0("working-", unname(tools::md5sum(fingerprint)))
  .assessment_database(function(name) read.csv(file.path(database_dir, name),
    stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE),
    source = "working_tree", commit = .database_head(database_dir),
    reference = "working_tree", key = key)
}

read_committed_database <- function(revision = NULL, repo = saio_path(),
                                    tier = "curated") {
  metadata <- .saio_metadata(revision, repo, tier)
  .assessment_database(function(name) .saio_read_table(name, repo, metadata$commit),
    source = metadata$database_source, commit = metadata$commit,
    reference = metadata$database_ref, selection = metadata$database_selection,
    key = metadata$database_key)
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
