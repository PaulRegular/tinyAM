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

# Keep only the most recently read database in memory for repeated lookups.
if (!exists(".saio_read_cache", inherits = FALSE)) .saio_read_cache <- new.env(parent = emptyenv())

# Keep the legacy snapshot and cache key used by the original 28-stock comparison.
.saio_comparison_revision <- "c36cc22df29dc0efb9a96b7e2500caf4bec36c86"
.saio_tables <- c("stocks", "assessments", "assumptions", "inputs", "outputs")

.saio_release_source <- function(release = NULL, database_dir = NULL) {
  if (!requireNamespace("saio", quietly = TRUE)) {
    cli::cli_abort(c("Install the SAIO package to read the assessment database.",
      "i" = 'Run {.code pak::pak("PaulRegular/saio")} in R.'))
  }
  if (is.null(database_dir)) {
    database_dir <- Sys.getenv("SAIO_DATABASE_DIR")
    if (!nzchar(database_dir)) {
      database_dir <- file.path(tools::R_user_dir("saio", "data"), "release")
    }
  }
  if (!is.character(database_dir) || length(database_dir) != 1L || is.na(database_dir) ||
      !nzchar(database_dir)) {
    cli::cli_abort("{.arg database_dir} must be one non-empty directory path.")
  }
  database_dir <- normalizePath(database_dir, winslash = "/", mustWork = FALSE)
  path <- file.path(database_dir, "manifest.csv")
  if (!file.exists(path)) {
    cli::cli_abort(c("No SAIO release manifest found in {.path {database_dir}}.",
      "i" = "Run {.code saio::download_database()} or point {.envvar SAIO_DATABASE_DIR} to the downloaded files."))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  manifest <- tryCatch(read.csv(text = lines, colClasses = "character",
    check.names = FALSE), error = function(e) NULL)
  tag <- unique(manifest$release)
  if (length(tag) != 1L || is.na(tag) || !nzchar(tag)) {
    cli::cli_abort("Invalid release tag in SAIO manifest {.file {path}}.")
  }
  if (!is.null(release) && !identical(release, tag)) {
    cli::cli_abort(c("The selected SAIO folder contains release {.val {tag}}.",
      "i" = "Point {.arg database_dir} to the downloaded copy of {.val {release}}."))
  }

  # Keep the original Git manifest key for saved fits, including on Windows.
  bytes <- charToRaw(paste0(paste(lines, collapse = "\n"), "\n"))
  blob <- c(charToRaw(paste("blob", length(bytes))), as.raw(0L), bytes)
  key <- unclass(as.character(openssl::sha1(blob)))
  manifest_commits <- c(
    "4cc958a3e2f7b51cfaf432298898e0d468590866" = "a34e616df69ad40f34be8b536531d67514b614f4",
    "00dd3b63d9934f01425ad1fa678f751a539c827b" = "14e86b223f04c376486b610a625efd2f35a48130"
  )
  list(database_dir = database_dir, release = tag, key = key,
    commit = unname(manifest_commits[key]), manifest = manifest)
}

.read_saio_release <- function(source, assessment_id = NULL, tables = .saio_tables) {
  if (!is.character(tables) || !length(tables) || anyDuplicated(tables) ||
      !all(tables %in% .saio_tables)) {
    cli::cli_abort("Choose from these SAIO tables: {.val {.saio_tables}}.")
  }
  files <- file.path(source$database_dir, c("manifest.csv", paste0(.saio_tables, ".csv.gz")))
  stamps <- c("size", "isdir", "mtime", "ctime")
  selection <- list(key = source$key, files = file.info(files)[, stamps, drop = FALSE])
  if (identical(selection, .saio_read_cache$selection)) {
    database <- .saio_read_cache$data
    if (!is.null(assessment_id)) {
      database <- .subset_saio_assessment(assessment_id, database)
    }
    return(database[c(tables, setdiff(names(database), .saio_tables))])
  }
  full <- is.null(assessment_id) && identical(tables, .saio_tables)
  if (full) .saio_read_cache$selection <- .saio_read_cache$data <- NULL
  database <- saio::read_database(database_dir = source$database_dir,
    assessment_id = assessment_id, tables = tables)
  if (!identical(database$database_release, source$release) ||
      !identical(selection$files, file.info(files)[, stamps, drop = FALSE])) {
    cli::cli_abort("The SAIO release files changed while being read. Verify the download and retry.")
  }
  if (is.null(assessment_id)) {
    expected <- source$manifest$rows[match(tables, source$manifest$table)]
    actual <- vapply(database[tables], nrow, integer(1L))
    wrong <- tables[actual != as.numeric(expected)]
    if (length(wrong)) {
      cli::cli_abort("SAIO table row counts differ from the release manifest: {.val {wrong}}.")
    }
  }
  database$database_source <- "saio"
  database$commit <- source$commit
  database$database_ref <- source$release
  database$database_key <- source$key
  database$database_selection <- "local_release"
  if (full) {
    .saio_read_cache$selection <- selection
    .saio_read_cache$data <- database
  }
  database
}

read_database <- function(saio_repo = NULL, revision = NULL, release = NULL,
                          database_dir = NULL, assessment_id = NULL, tables = .saio_tables) {
  if (!is.null(assessment_id)) .saio_assessment_id(assessment_id)
  if (!is.null(revision)) {
    if (!is.null(release) || !is.null(database_dir) || !is.null(assessment_id) ||
        !identical(tables, .saio_tables)) {
      cli::cli_abort("Read a historical SAIO revision in full, then select the assessment from it.")
    }
    if (is.null(saio_repo)) saio_repo <- .saio_repo()
    return(.read_saio_revision(saio_repo, revision))
  }
  if (!is.null(saio_repo)) {
    cli::cli_abort("Use {.arg saio_repo} with a historical {.arg revision}, or {.arg database_dir} for a downloaded release.")
  }
  .read_saio_release(.saio_release_source(release, database_dir), assessment_id, tables)
}

# Read committed Git blobs or their local Git LFS objects for a legacy revision.
.read_saio_revision <- function(saio_repo, revision) {
  commit <- .saio_commit(saio_repo, revision)
  repo <- normalizePath(saio_repo, winslash = "/", mustWork = TRUE)
  selection <- list(repo = repo, commit = commit)
  if (identical(selection, .saio_read_cache$selection)) return(.saio_read_cache$data)
  .saio_read_cache$selection <- .saio_read_cache$data <- NULL

  blob <- tempfile("saio-table-")
  on.exit(unlink(blob))
  media_dir <- NULL
  tables <- lapply(.saio_tables, function(name) {
    file <- paste0("database/", name, ".csv")
    status <- .saio_git(repo, c("cat-file", "blob", paste0(commit, ":", file)), stdout = blob)
    if (status != 0L) {
      cli::cli_abort(c("SAIO commit {.val {commit}} does not provide legacy table {.file {file}}.",
        "i" = "For current SAIO data, run {.code saio::download_database()} and read it without {.arg revision}."))
    }
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
  names(tables) <- .saio_tables
  tables$database_source <- "saio"
  tables$commit <- tables$database_ref <- tables$database_key <- commit
  tables$database_selection <- "local_commit"
  .saio_read_cache$selection <- selection
  .saio_read_cache$data <- tables
  tables
}

.saio_assessment_id <- function(assessment_id) {
  if (!is.character(assessment_id) || length(assessment_id) != 1L ||
      is.na(assessment_id) || !nzchar(assessment_id)) {
    cli::cli_abort("{.arg assessment_id} must be one non-empty assessment ID.")
  }
}

.subset_saio_assessment <- function(assessment_id, database) {
  if (!all(.saio_tables %in% names(database))) {
    cli::cli_abort("{.arg database} must contain all five SAIO tables.")
  }
  database$assessments <- database$assessments[
    which(database$assessments$assessment_id == assessment_id), , drop = FALSE]
  if (nrow(database$assessments) != 1L) {
    cli::cli_abort("The database must contain exactly one assessment with ID {.val {assessment_id}}.")
  }
  stock_id <- database$assessments$stock_id[[1L]]
  database$stocks <- database$stocks[
    which(database$stocks$stock_id == stock_id), , drop = FALSE]
  for (table in c("assumptions", "inputs", "outputs")) {
    database[[table]] <- database[[table]][
      which(database[[table]]$assessment_id == assessment_id), , drop = FALSE]
  }
  for (table in .saio_tables) rownames(database[[table]]) <- NULL
  database
}

read_assessment <- function(assessment_id, database = NULL) {
  .saio_assessment_id(assessment_id)
  tables <- if (is.null(database)) read_database(assessment_id = assessment_id) else
    .subset_saio_assessment(assessment_id, database)
  # A recipe using numeric bins still receives numbers when all labels are numeric.
  bins <- intersect(c("length_bin", "weight_bin"), names(tables$inputs))
  tables$inputs[bins] <- lapply(tables$inputs[bins], utils::type.convert, as.is = TRUE)
  age <- tables$inputs$age
  if (is.double(age) && all(is.na(age) | age == trunc(age))) {
    tables$inputs$age <- as.integer(age)
  }
  tables$assessment <- tables$assessments
  tables
}
