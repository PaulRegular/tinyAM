.database_head <- function(repo = ".") {
  value <- suppressWarnings(tryCatch(system2("git", c("-C", shQuote(repo),
    "rev-parse", "--verify", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()))
  if (length(value) != 1L || !grepl("^[[:xdigit:]]{40}$", value)) return(NA_character_)
  value[[1]]
}

.saio_git <- function(repo, args, stdout = TRUE) {
  suppressWarnings(system2("git", c("-C", shQuote(repo), shQuote(args)),
    stdout = stdout, stderr = FALSE))
}

.saio_repo <- function() {
  configured <- Sys.getenv("SAIO_REPO")
  if (nzchar(configured)) return(configured)
  # Resolve the sibling of the primary tinyAM checkout from linked worktrees too.
  common_dir <- .saio_git(".", c("rev-parse", "--path-format=absolute", "--git-common-dir"))
  if (length(common_dir) == 1L) file.path(dirname(dirname(common_dir)), "saio") else "../saio"
}

.saio_commit <- function(saio_repo = .saio_repo(), revision = "HEAD") {
  if (!dir.exists(saio_repo)) cli::cli_abort("Set {.envvar SAIO_REPO} to your local SAIO repository. Could not find {.path {saio_repo}}.")
  commit <- .saio_git(saio_repo, c("rev-parse", "--verify", "--end-of-options", paste0(revision, "^{commit}")))
  if (length(commit) != 1L || !grepl("^[[:xdigit:]]{40}$", commit)) {
    cli::cli_abort("Cannot resolve SAIO revision {.val {revision}} in {.path {saio_repo}}. Check the local Git repository and revision.")
  }
  commit
}

# Read one committed revision. Local Git LFS objects are the CSV files for SAIO.
# Keep only the most recently read database in memory for repeated lookups.
if (!exists(".saio_read_cache", inherits = FALSE)) .saio_read_cache <- new.env(parent = emptyenv())

read_database <- function(saio_repo = .saio_repo(), revision = "HEAD") {
  commit <- .saio_commit(saio_repo, revision)
  repo <- normalizePath(saio_repo, winslash = "/", mustWork = TRUE)
  if (identical(repo, .saio_read_cache$repo) && identical(commit, .saio_read_cache$data$commit)) {
    return(.saio_read_cache$data)
  }

  table_names <- c("stocks", "assessments", "assumptions", "inputs", "outputs")
  blob <- tempfile("saio-table-")
  on.exit(unlink(blob))
  media_dir <- NULL
  tables <- lapply(table_names, function(name) {
    file <- paste0("database/", name, ".csv")
    status <- .saio_git(repo, c("cat-file", "blob", paste0(commit, ":", file)), stdout = blob)
    if (status != 0L) cli::cli_abort("SAIO commit {.val {commit}} does not provide {.file {file}}.")
    header <- if (file.size(blob) < 1024) readLines(blob, n = 3L, warn = FALSE) else character()
    path <- blob
    if (identical(header[1L], "version https://git-lfs.github.com/spec/v1")) {
      if (is.null(media_dir)) {
        env <- .saio_git(repo, c("lfs", "env"))
        media_dir <- sub("^LocalMediaDir=", "", grep("^LocalMediaDir=", env, value = TRUE))
        if (length(media_dir) != 1L) cli::cli_abort("Install Git LFS to read the local SAIO data.")
      }
      oid <- sub("^oid sha256:", "", header[2L])
      size <- as.numeric(sub("^size ", "", header[3L]))
      path <- file.path(media_dir, substr(oid, 1L, 2L), substr(oid, 3L, 4L), oid)
      if (.Platform$OS.type == "windows" && nchar(path) >= 260L) {
        # The LFS hash directories can push an otherwise valid path past MAX_PATH.
        path <- chartr("/", "\\", path)
        if (grepl("^[A-Za-z]:", path)) {
          path <- paste0("\\\\?\\", path)
        } else if (startsWith(path, "\\\\") && !startsWith(path, "\\\\?\\")) {
          path <- paste0("\\\\?\\UNC\\", substring(path, 3L))
        }
      }
      if (!isTRUE(file.size(path) == size)) {
        cli::cli_abort(c("The local Git LFS object for {.file {file}} is missing or incomplete.",
          "i" = "Fetch the data from SAIO with {.code git lfs fetch origin {commit}} and retry."))
      }
    }
    read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE)
  })
  names(tables) <- table_names
  tables$database_source <- "saio"
  tables$commit <- tables$database_ref <- tables$database_key <- commit
  tables$database_selection <- "local_commit"
  .saio_read_cache$repo <- repo
  .saio_read_cache$data <- tables
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
