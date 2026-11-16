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

.saio_release_source <- function(saio_repo = .saio_repo(), release = NULL, database_dir = NULL) {
  head <- .saio_commit(saio_repo)
  repo <- normalizePath(saio_repo, winslash = "/", mustWork = TRUE)
  paths <- .saio_git(repo, c("ls-tree", "-r", "--name-only", head, "--", "database/releases"))
  paths <- paths[grepl("^database/releases/[^/]+[.]csv$", paths)]
  # SAIO names releases with a sortable ISO date prefix.
  available <- sort(sub("[.]csv$", "", basename(paths)), decreasing = TRUE, method = "radix")
  if (!length(available)) {
    cli::cli_abort(c("No committed SAIO release manifests were found in {.path {repo}}.",
      "i" = "Update the SAIO checkout, or select a legacy Git snapshot with {.arg revision}."))
  }
  if (is.null(release)) release <- available[[1L]]
  if (!is.character(release) || length(release) != 1L || is.na(release) ||
      !release %in% available) {
    cli::cli_abort(c("Select a committed SAIO release by its name without {.file .csv}.",
      "i" = "Available releases: {.val {available}}."))
  }
  path <- paste0("database/releases/", release, ".csv")
  ref <- paste0(head, ":", path)
  key <- .saio_git(repo, c("rev-parse", "--verify", "--end-of-options", ref))
  commit <- .saio_git(repo, c("log", "-1", "--format=%H", head, "--", path))
  if (length(key) != 1L || !grepl("^[[:xdigit:]]{40}$", key) ||
      length(commit) != 1L || !grepl("^[[:xdigit:]]{40}$", commit)) {
    cli::cli_abort("Could not verify the Git history for SAIO manifest {.file {path}}.")
  }
  csv <- .saio_git(repo, c("cat-file", "blob", key))
  manifest <- tryCatch(read.csv(text = csv, colClasses = "character",
    na.strings = NULL, check.names = FALSE), error = function(e) NULL)
  fields <- c("release", "schema_version", "table", "asset", "rows", "csv_bytes",
    "csv_sha256", "asset_bytes", "asset_sha256")
  if (!identical(names(manifest), fields) || nrow(manifest) != length(.saio_tables)) {
    cli::cli_abort("Invalid or unsupported SAIO release manifest {.file {path}}.")
  }
  numbers <- unlist(manifest[c("rows", "csv_bytes", "asset_bytes")], use.names = FALSE)
  hashes <- unlist(manifest[c("csv_sha256", "asset_sha256")], use.names = FALSE)
  valid <- setequal(manifest$table, .saio_tables) && isTRUE(all(manifest$release == release)) &&
    isTRUE(all(manifest$schema_version == "1")) &&
    isTRUE(all(manifest$asset == paste0(manifest$table, ".csv.gz"))) &&
    all(grepl("^[0-9]+$", numbers)) && all(is.finite(suppressWarnings(as.numeric(numbers)))) &&
    all(grepl("^[[:xdigit:]]{64}$", hashes))
  if (!valid) cli::cli_abort("Invalid or unsupported SAIO release manifest {.file {path}}.")
  manifest <- manifest[match(.saio_tables, manifest$table), , drop = FALSE]

  if (is.null(database_dir)) {
    database_dir <- Sys.getenv("SAIO_DATABASE_DIR", unset = file.path(repo, "database"))
    if (!nzchar(database_dir)) database_dir <- file.path(repo, "database")
  }
  if (!is.character(database_dir) || length(database_dir) != 1L || is.na(database_dir) ||
      !nzchar(database_dir)) {
    cli::cli_abort("{.arg database_dir} must be one non-empty directory path.")
  }
  list(repo = repo, release = release, key = key, commit = commit,
    database_dir = normalizePath(database_dir, winslash = "/", mustWork = FALSE),
    manifest = manifest)
}

.saio_sha256 <- function(path) {
  stream <- file(path, "rb")
  on.exit(close(stream))
  unclass(as.character(openssl::sha256(stream)))
}

.read_saio_csv <- function(path, table) {
  # SAIO identifiers and labels can look numeric, so keep text fields as text.
  types <- readr::cols(
    .default = readr::col_character(),
    assessment_year = readr::col_integer(), terminal_year = readr::col_integer(),
    estimate_terminal_year = readr::col_integer(), framework_year = readr::col_integer(),
    is_current = readr::col_logical(), is_applied = readr::col_logical(),
    year = readr::col_integer(), age = readr::col_integer(),
    value = readr::col_double(), sampling_time = readr::col_double(),
    length_bin_lower = readr::col_double(),
    length_bin_upper = readr::col_double(),
    sample_size = readr::col_double(), se = readr::col_double(),
    lwr = readr::col_double(), upr = readr::col_double()
  )
  if (table == "assumptions") types$cols$value <- NULL
  # Ageing-error SDs in SAIO inputs can refer to half-year midpoints.
  if (table == "inputs") types$cols$age <- readr::col_double()
  columns <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  types$cols <- types$cols[intersect(names(types$cols), columns)]
  data <- suppressWarnings(readr::with_edition(2L, readr::read_csv(
    path, col_types = types, name_repair = "minimal", trim_ws = FALSE,
    num_threads = 1L, progress = FALSE, lazy = FALSE)))
  problems <- readr::problems(data)
  # Some source exports omit trailing, optional input fields. The reader
  # fills those with NA, as in the SAIO validator's base-R reader.
  optional <- c("observation_id", "length_bin", "length_bin_lower", "length_bin_upper",
    "sample_size", "age_error", "partition", "weight_bin")
  if (table == "inputs" && nrow(problems)) {
    short <- grepl("^[0-9]+ columns?$", problems$actual) &
      problems$expected == paste(ncol(data), "columns")
    fields <- as.integer(sub(" columns?$", "", problems$actual[short]))
    trailing <- vapply(fields, function(n) {
      n < length(columns) && all(columns[seq.int(n + 1L, length(columns))] %in% optional)
    }, logical(1L))
    if (any(trailing)) {
      accepted <- which(short)[trailing]
      # readr can retain a CR in the final present text field of a short
      # CRLF record. Its row numbers include the header.
      for (field in unique(fields[trailing])) {
        rows <- problems$row[accepted[fields[trailing] == field]] - 1L
        if (!is.character(data[[field]])) next
        value <- sub("\r$", "", data[[field]][rows])
        data[[field]][rows] <- replace(value, !is.na(value) & value == "", NA_character_)
      }
      problems <- problems[-accepted, , drop = FALSE]
    }
  }
  if (nrow(problems)) {
    first <- problems[1L, ]
    cli::cli_abort(c("Cannot parse SAIO table {.file {path}}.",
      "x" = "At row {first$row}, column {first$col}: expected {first$expected}, found {.val {first$actual}}."))
  }
  as.data.frame(data)
}

.read_saio_release <- function(source) {
  manifest <- source$manifest
  release <- source$release
  files <- file.path(source$database_dir, paste0(manifest$table, ".csv"))
  file_stats <- file.info(files)[, c("size", "isdir", "mtime", "ctime"), drop = FALSE]
  bad <- which(is.na(file_stats$size) | file_stats$isdir |
    file_stats$size != as.numeric(manifest$csv_bytes))
  if (length(bad)) {
    cli::cli_abort(c("Missing or size-mismatched tables for SAIO release {.val {release}}.",
      "x" = "{.file {files[bad]}}",
      "i" = "Restore the selected SAIO release or set {.envvar SAIO_DATABASE_DIR} to its restored tables."))
  }
  selection <- list(repo = source$repo, commit = source$commit, key = source$key,
    files = file_stats)
  if (identical(selection, .saio_read_cache$selection)) return(.saio_read_cache$data)
  .saio_read_cache$selection <- .saio_read_cache$data <- NULL
  if (!requireNamespace("openssl", quietly = TRUE) || !requireNamespace("readr", quietly = TRUE)) {
    cli::cli_abort("Install the R packages {.pkg openssl} and {.pkg readr} to verify and read SAIO releases.")
  }
  valid <- vapply(files, .saio_sha256, character(1L)) == tolower(manifest$csv_sha256)
  if (!all(valid)) {
    cli::cli_abort(c("SHA-256 checksum mismatch for SAIO release {.val {release}}.",
      "x" = "{.file {files[!valid]}}",
      "i" = "Restore the selected SAIO release or use its verified database directory."))
  }
  tables <- stats::setNames(lapply(seq_along(files), function(i) {
    path <- files[[i]]
    table <- .read_saio_csv(path, manifest$table[i])
    if (nrow(table) != as.numeric(manifest$rows[i])) {
      cli::cli_abort("Row count in {.file {path}} does not match SAIO release {.val {release}}.")
    }
    required <- switch(manifest$table[i], stocks = "stock_id",
      assessments = c("assessment_id", "stock_id"), "assessment_id")
    if (anyDuplicated(names(table)) || !all(required %in% names(table))) {
      cli::cli_abort("Table {.file {path}} does not have the expected SAIO columns.")
    }
    table
  }), manifest$table)
  if (!identical(file_stats, file.info(files)[, names(file_stats), drop = FALSE])) {
    cli::cli_abort("SAIO release files changed while being read. Finish restoring the release and retry.")
  }
  tables$database_source <- "saio"
  tables$commit <- source$commit
  tables$database_release <- tables$database_ref <- release
  tables$database_key <- source$key
  tables$database_selection <- "local_release"
  .saio_read_cache$selection <- selection
  .saio_read_cache$data <- tables
  tables
}

read_database <- function(saio_repo = .saio_repo(), revision = NULL, release = NULL,
                          database_dir = NULL) {
  if (!is.null(revision)) {
    if (!is.null(release) || !is.null(database_dir)) {
      cli::cli_abort("Use {.arg revision} for a legacy SAIO snapshot, or use {.arg release} and {.arg database_dir} for restored release tables.")
    }
    return(.read_saio_revision(saio_repo, revision))
  }
  .read_saio_release(.saio_release_source(saio_repo, release, database_dir))
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
        "i" = "For current SAIO data, restore a release and read it without {.arg revision}."))
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
  # SAIO preserves source bin labels (including intervals such as "25-40 cm").
  # A recipe using numeric bins still receives numbers when all its labels are numeric.
  bins <- intersect(c("length_bin", "weight_bin"), names(tables$inputs))
  tables$inputs[bins] <- lapply(tables$inputs[bins], utils::type.convert, as.is = TRUE)
  age <- tables$inputs$age
  if (is.double(age) && all(is.na(age) | age == trunc(age))) {
    tables$inputs$age <- as.integer(age)
  }
  tables$assessment <- assessment
  tables
}
