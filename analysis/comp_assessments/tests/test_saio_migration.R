source("analysis/comp_assessments/R/run_assessment.R")
registry <- assessment_translations()
stopifnot(all(file.exists(registry$script)),
  !anyDuplicated(paste(registry$assessment_id, registry$tier)),
  nrow(assessment_translations("curated")) == 41L)

local({
  original_wd <- getwd()
  original_database_dir <- Sys.getenv("SAIO_DATABASE_DIR", unset = NA_character_)
  on.exit({
    if (is.na(original_database_dir)) {
      Sys.unsetenv("SAIO_DATABASE_DIR")
    } else {
      Sys.setenv(SAIO_DATABASE_DIR = original_database_dir)
    }
  }, add = TRUE)
  Sys.unsetenv("SAIO_DATABASE_DIR")
  fixture <- tempfile("saio_repository_")
  dir.create(fixture)
  on.exit({setwd(original_wd); unlink(fixture, recursive = TRUE)}, add = TRUE)
  relative_root <- file.path("analysis", "comp_assessments")
  dir.create(file.path(fixture, relative_root, "R"), recursive = TRUE)
  file.copy(list.files(file.path(relative_root, "R"), full.names = TRUE),
    file.path(fixture, relative_root, "R"))
  setwd(fixture)
  dir.create("data_repo/database", recursive = TRUE)
  git <- function(...) {
    result <- system2("git", c("-C", shQuote("data_repo"), shQuote(c(...))),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"))
    result
  }
  stocks <- data.frame(stock_id = "source_stock", scientific_name = "Source species")
  assessments <- data.frame(assessment_id = c("shared", "extra"),
    stock_id = "source_stock", is_current = c(FALSE, TRUE))
  tables <- list(stocks = stocks, assessments = assessments,
    assumptions = data.frame(assessment_id = assessments$assessment_id,
      component = "population", setting = "age_range", value = "1-2",
      source_reference = "fixture", notes = "fixture"),
    inputs = data.frame(assessment_id = assessments$assessment_id),
    outputs = data.frame(assessment_id = assessments$assessment_id, source_marker = "Committed output"))
  for (name in names(tables)) write.csv(tables[[name]],
    file.path("data_repo/database", paste0(name, ".csv")), row.names = FALSE)
  git("init", "-b", "main")
  lfs_version <- suppressWarnings(system2("git", c("lfs", "version"),
    stdout = TRUE, stderr = TRUE))
  has_lfs <- is.null(attr(lfs_version, "status"))
  store_lfs_output <- function() {
    if (!has_lfs) return(invisible(NULL))
    # Generate a real local LFS object without changing Git filter configuration.
    path <- "data_repo/database/outputs.csv"
    pointer <- system2("git", c("-C", shQuote("data_repo"), "lfs", "clean", "--",
      "database/outputs.csv"), stdin = path, stdout = TRUE)
    stopifnot(is.null(attr(pointer, "status")),
      identical(pointer[1L], "version https://git-lfs.github.com/spec/v1"))
    writeLines(pointer, path, useBytes = TRUE)
  }
  store_lfs_output()
  git("add", "database")
  commit <- function(text) git("-c", "user.name=SourceTest", "-c", "user.email=source@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", text)
  commit("fixture")
  revision <- git("rev-parse", "HEAD")[[1]]
  if (has_lfs) stopifnot(identical(
    git("show", "HEAD:database/outputs.csv")[[1]],
    "version https://git-lfs.github.com/spec/v1"
  ))

  # The first read must use the commit even when the working tree is dirty.
  stocks$scientific_name <- "Uncommitted species"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  write.csv(transform(tables$outputs, source_marker = "Uncommitted output"),
    "data_repo/database/outputs.csv", row.names = FALSE)
  pinned <- read_database(saio_repo = "data_repo", revision = "HEAD")
  stopifnot(pinned$stocks$scientific_name == "Source species",
    all(pinned$outputs$source_marker == "Committed output"),
    all(names(tables) %in% names(pinned)),
    identical(pinned$commit, revision), identical(pinned$database_ref, revision),
    identical(pinned$database_key, revision),
    identical(pinned$database_source, "saio"),
    identical(pinned$database_selection, "local_commit"))

  stocks$scientific_name <- "New committed species"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  write.csv(transform(tables$outputs, source_marker = "New committed output"),
    "data_repo/database/outputs.csv", row.names = FALSE)
  store_lfs_output()
  git("add", "database")
  commit("new data")
  next_revision <- git("rev-parse", "HEAD")[[1]]
  stopifnot(file.rename("data_repo/database/outputs.csv", "data_repo/database/outputs.saved"))

  current <- read_database(saio_repo = "data_repo", revision = "HEAD")
  stopifnot(current$commit != revision, identical(current$commit, next_revision),
    current$stocks$scientific_name == "New committed species",
    all(current$outputs$source_marker == "New committed output"),
    identical(read_database(saio_repo = "data_repo", revision = revision), pinned),
    identical(read_database(saio_repo = "data_repo", revision = "HEAD~1"), pinned),
    !dir.exists(file.path(relative_root, "cache", "saio")))
  invalid <- tryCatch(read_database(saio_repo = "data_repo", revision = "missing-ref"),
    error = identity)
  stopifnot(inherits(invalid, "error"))
  entries <- data.frame(assessment_id = c("shared", "shared", "extra"),
    tier = c("curated", "preliminary", "preliminary"), source_commit = revision)
  write.csv(entries, file.path(relative_root, "translations.csv"), row.names = FALSE)
  for (tier in c("curated", "preliminary")) {
    folder <- file.path(relative_root, "scripts", tier)
    dir.create(folder, recursive = TRUE)
    for (id in entries$assessment_id[entries$tier == tier]) writeLines(c(
      paste0('obs <- list(marker = "', tier, '")'),
      'fit <- NULL', 'background <- "Fixture translation"'
    ), file.path(folder, paste0(id, ".R")))
  }
  source(file.path(relative_root, "R", "run_assessment.R"))
  preferred <- run_assessment("shared", pinned, fit = FALSE)
  alternative <- run_assessment("shared", pinned, fit = FALSE, tier = "preliminary")
  stopifnot(preferred$obs$marker == "curated",
    alternative$obs$marker == "preliminary",
    preferred$diagnostics$translation_tier == "curated",
    alternative$diagnostics$translation_tier == "preliminary",
    preferred$diagnostics$database_revision == revision,
    grepl(paste0("results/curated/cache/", revision, "/shared$"), .assessment_cache("shared", "curated", pinned)),
    grepl(paste0("results/preliminary/cache/", revision, "/shared$"), .assessment_cache("shared", "preliminary", pinned)))
  historical_cache <- .assessment_cache("shared", "curated", pinned)
  current_cache <- .assessment_cache("shared", "curated", current)
  legacy_cache <- file.path(.assessment_results("curated"), "cache", "shared")
  for (folder in c(historical_cache, legacy_cache)) dir.create(folder, recursive = TRUE)
  old_fit <- structure(list(marker = "historical", assessment_provenance = list(database_key = revision)), class = "tam_fit")
  saveRDS(old_fit, file.path(historical_cache, "fit.rds"))
  saveRDS(old_fit, file.path(legacy_cache, "fit.rds"))
  original_cache_hash <- tools::md5sum(file.path(c(historical_cache, legacy_cache), "fit.rds"))
  stopifnot(is.null(.assessment_cached_fit("shared", "curated", current)),
    .assessment_cached_fit("shared", "curated", pinned)$marker == "historical",
    current_cache != historical_cache,
    .assessment_cache_key(pinned) != .assessment_cache_key(current))
  dir.create(current_cache, recursive = TRUE)
  saveRDS(old_fit, file.path(current_cache, "fit.rds"))
  mismatched <- tryCatch(.assessment_cached_fit("shared", "curated", current), error = identity)
  stopifnot(inherits(mismatched, "error"),
    grepl("lacks matching source provenance", conditionMessage(mismatched), fixed = TRUE),
    identical(tools::md5sum(file.path(c(historical_cache, legacy_cache), "fit.rds")), original_cache_hash))
  curated <- run_assessments(database = pinned, fit = FALSE,
    progress = FALSE, save_results = TRUE)
  preliminary <- run_assessments(database = pinned, tier = "preliminary",
    fit = FALSE, progress = FALSE, save_results = TRUE)
  stopifnot(identical(names(curated), "shared"),
    setequal(names(preliminary), c("shared", "extra")),
    !"extra" %in% read.csv(file.path(.assessment_results("curated"),
      "fit_diagnostics.csv"))$assessment_id,
    setequal(read.csv(file.path(.assessment_results("preliminary"),
      "fit_diagnostics.csv"))$assessment_id, c("shared", "extra")))
  stopifnot(isTRUE(.assessment_check_results_revision(.assessment_results("curated"), pinned)))
  stale <- tryCatch(.assessment_check_results_revision(.assessment_results("curated"), current), error = identity)
  stopifnot(inherits(stale, "error"), grepl("another or unknown SAIO revision", conditionMessage(stale), fixed = TRUE))
  latest_result <- run_assessment("shared", current, fit = FALSE)
  mixed_revision <- tryCatch(.assessment_result_tables(list(preferred, latest_result),
    save_results = TRUE), error = identity)
  stopifnot(inherits(mixed_revision, "error"), grepl("mixed SAIO revisions", conditionMessage(mixed_revision), fixed = TRUE))
  .assessment_save_rds(list(shared = preferred), file.path(.assessment_cache_root("curated", pinned), "assessment_runs.rds"))
  stopifnot(.assessment_cached_runs("curated", pinned)$shared$diagnostics$database_revision == revision,
    is.null(.assessment_cached_runs("curated", current)))
  .assessment_save_rds(list(shared = preferred), file.path(.assessment_cache_root("curated", current), "assessment_runs.rds"))
  wrong_batch <- tryCatch(.assessment_cached_runs("curated", current), error = identity)
  stopifnot(inherits(wrong_batch, "error"),
    grepl("do not match the selected SAIO revision", conditionMessage(wrong_batch), fixed = TRUE))
  rejected <- tryCatch(run_assessments("extra", pinned,
    fit = FALSE, progress = FALSE), error = identity)
  mixed <- tryCatch(.assessment_result_tables(list(preferred, alternative),
    save_results = TRUE), error = identity)
  stopifnot(inherits(rejected, "error"), inherits(mixed, "error"))
  for (tier in c("curated", "preliminary")) {
    folder <- file.path(.assessment_results(tier), "dashboards")
    dir.create(folder, recursive = TRUE)
    writeLines("<html>Fixture</html>", file.path(folder, paste0(tier, "_only.html")))
    trial <- file.path(.assessment_results(tier), "cache", "shared", "experiments", "setup")
    dir.create(trial, recursive = TRUE)
    writeLines("<html>Retired setup</html>", file.path(trial, "dashboard.html"))
    for (database in list(pinned, current)) {
      cache <- .assessment_cache("shared", tier, database)
      dir.create(cache, recursive = TRUE, showWarnings = FALSE)
      writeLines("<html>Cached dashboard</html>", file.path(cache, "dashboard.html"))
    }
    .index_assessment_dashboards(tier, current)
  }
  curated_index <- readLines(file.path(.assessment_results("curated"), "index.html"))
  preliminary_index <- readLines(file.path(.assessment_results("preliminary"), "index.html"))
  stopifnot(any(grepl("curated_only.html", curated_index, fixed = TRUE)),
    !any(grepl("preliminary_only.html", curated_index, fixed = TRUE)),
    any(grepl("preliminary_only.html", preliminary_index, fixed = TRUE)),
    !any(grepl("curated_only.html", preliminary_index, fixed = TRUE)),
    any(grepl(paste0("selected source ", current$commit), curated_index, fixed = TRUE)),
    any(grepl(paste0("historical source ", pinned$commit), curated_index, fixed = TRUE)),
    any(grepl("historical output; source revision unknown", curated_index, fixed = TRUE)),
    !any(grepl("/experiments/", curated_index, fixed = TRUE)),
    !any(grepl("/experiments/", preliminary_index, fixed = TRUE)))

  # Released data can be read directly from compressed files outside Git.
  stopifnot(requireNamespace("saio", quietly = TRUE), requireNamespace("openssl", quietly = TRUE))
  stopifnot(identical(.saio_comparison_revision,
    "c36cc22df29dc0efb9a96b7e2500caf4bec36c86"))
  dir.create("data_repo/database/releases")
  release_ids <- c("002", "003")
  source_note <- 'Reported as "main series", then checked.\nReview appendix.'
  csv_bytes <- function(lines, eol) charToRaw(paste0(paste(lines, collapse = eol), eol))
  digest <- function(x) unclass(as.character(openssl::sha256(x)))
  write_gzip <- function(bytes, path) {
    connection <- gzfile(path, "wb")
    on.exit(close(connection))
    writeBin(bytes, connection)
  }
  write_release <- function(tag, destination, species) {
    dir.create(destination, showWarnings = FALSE)
    release_tables <- tables
    release_tables$stocks <- transform(release_tables$stocks, stock_id = "001", scientific_name = species)
    release_tables$assessments <- transform(release_tables$assessments,
      stock_id = "001", assessment_id = release_ids)
    release_tables$assumptions <- transform(release_tables$assumptions,
      assessment_id = release_ids, value = "03.00", notes = source_note)
    release_tables$inputs <- transform(release_tables$inputs,
      assessment_id = release_ids, type = "ageing_error_sd", age = c(0.5, 1.5),
      value = c(0.15, 0.25), length_bin = c("25-40 cm", "20"), weight_bin = NA_real_)
    release_tables$outputs$assessment_id <- release_ids
    lines <- lapply(release_tables, function(table) capture.output(write.csv(table, row.names = FALSE)))
    # A release may omit trailing optional fields in individual input rows.
    lines$inputs[[2L]] <- sub(",NA$", "", lines$inputs[[2L]])
    csv <- lapply(lines, csv_bytes, eol = "\r\n")
    asset_names <- paste0(names(release_tables), ".csv.gz")
    files <- file.path(destination, asset_names)
    for (i in seq_along(files)) write_gzip(csv[[i]], files[[i]])
    assets <- lapply(files, function(path) readBin(path, "raw", n = file.size(path)))
    manifest <- data.frame(release = tag, schema_version = 1, table = names(release_tables),
      asset = asset_names, rows = vapply(release_tables, nrow, integer(1)),
      csv_bytes = lengths(csv), csv_sha256 = vapply(csv, digest, character(1)),
      asset_bytes = lengths(assets), asset_sha256 = vapply(assets, digest, character(1)))
    manifest_lines <- capture.output(write.csv(manifest, row.names = FALSE))
    # Match the CRLF release asset and LF manifest blob used by published comparisons.
    writeBin(csv_bytes(manifest_lines, "\r\n"), file.path(destination, "manifest.csv"))
    path <- file.path("database", "releases", paste0(tag, ".csv"))
    writeBin(csv_bytes(manifest_lines, "\n"), file.path("data_repo", path))
    path
  }
  previous_tag <- "saio-2026-11-14-expansion"
  latest_tag <- "saio-2026-11-16-curation"
  previous_dir <- "previous data"
  latest_dir <- "current data"
  previous_path <- write_release(previous_tag, previous_dir, "Earlier released species")
  latest_path <- write_release(latest_tag, latest_dir, "Current released species")
  git("add", previous_path, latest_path)
  commit("released manifests")
  previous_blob <- git("rev-parse", paste0("HEAD:", previous_path))[[1L]]
  manifest_blob <- git("rev-parse", paste0("HEAD:", latest_path))[[1L]]

  Sys.setenv(SAIO_DATABASE_DIR = latest_dir)
  selected <- read_database(assessment_id = "002")
  latest <- read_database()
  stopifnot(identical(latest$stocks$scientific_name, "Current released species"),
    identical(latest$stocks$stock_id, "001"),
    identical(latest$assessments$stock_id, rep("001", 2L)),
    identical(latest$assessments$assessment_id, release_ids),
    identical(latest$inputs$assessment_id, release_ids),
    identical(latest$inputs$age, c(0.5, 1.5)),
    identical(latest$inputs$value, c(0.15, 0.25)),
    identical(latest$inputs$length_bin, c("25-40 cm", "20")),
    identical(read_assessment("002", latest)$inputs$length_bin, "25-40 cm"),
    identical(read_assessment("003", latest)$inputs$length_bin, 20L),
    identical(latest$assumptions$value, rep("03.00", 2L)),
    identical(gsub("\r\n", "\n", latest$assumptions$notes, fixed = TRUE), rep(source_note, 2L)),
    identical(latest$database_source, "saio"),
    identical(latest$database_selection, "local_release"),
    identical(latest$database_release, latest_tag),
    identical(latest$database_ref, latest_tag),
    identical(latest$commit, NA_character_),
    identical(latest$database_key, manifest_blob),
    identical(latest, read_database(release = latest_tag)),
    identical(.assessment_cache_key(), manifest_blob))
  stopifnot(identical(selected, read_database(assessment_id = "002")),
    identical(selected$assessments$assessment_id, "002"),
    identical(selected$stocks, latest$stocks),
    all(vapply(selected[c("assumptions", "inputs", "outputs")],
      function(table) identical(table$assessment_id, "002"), logical(1L))),
    identical(selected$database_key, manifest_blob))
  previous <- read_database(release = previous_tag, database_dir = previous_dir)
  stopifnot(identical(previous$stocks$scientific_name, "Earlier released species"),
    identical(previous$database_release, previous_tag),
    identical(previous$database_key, previous_blob),
    identical(pinned, read_database(saio_repo = "data_repo", revision = revision)),
    length(unique(vapply(list(latest, previous, pinned), .assessment_cache_key, character(1)))) == 3L)
  testthat::expect_error(read_database(release = previous_tag), previous_tag, fixed = TRUE)
  Sys.setenv(SAIO_DATABASE_DIR = previous_dir)
  stopifnot(identical(previous, read_database(release = previous_tag)),
    identical(latest, read_database(database_dir = latest_dir)))
  Sys.setenv(SAIO_DATABASE_DIR = latest_dir)

  # A fit saved with the old manifest key is still found for the same release.
  release_cache <- file.path(.assessment_results("curated"), "cache", manifest_blob, "002")
  dir.create(release_cache, recursive = TRUE)
  release_fit <- structure(list(marker = "released fit",
    assessment_provenance = list(database_key = manifest_blob)), class = "tam_fit")
  saveRDS(release_fit, file.path(release_cache, "fit.rds"))
  stopifnot(identical(.assessment_cached_fit("002", "curated", latest), release_fit),
    is.null(.assessment_cached_fit("002", "curated", previous)))

  # Plotting needs only the small catalogue tables, even if inputs are absent.
  released_inputs <- file.path(latest_dir, "inputs.csv.gz")
  stopifnot(file.rename(released_inputs, paste0(released_inputs, ".saved")))
  catalogue <- read_database(tables = c("stocks", "assessments"))
  stopifnot(identical(catalogue$stocks, latest$stocks),
    identical(catalogue$assessments, latest$assessments),
    identical(catalogue$database_key, manifest_blob),
    !any(c("assumptions", "inputs", "outputs") %in% names(catalogue)))
  testthat::expect_error(read_database(), "inputs\\.csv\\.gz")
  stopifnot(file.rename(paste0(released_inputs, ".saved"), released_inputs))

  # The parsed table still has to match the manifest's declared row count.
  manifest_path <- file.path(latest_dir, "manifest.csv")
  manifest_bytes <- readBin(manifest_path, "raw", n = file.size(manifest_path))
  wrong_rows <- read.csv(manifest_path)
  wrong_rows$rows[wrong_rows$table == "stocks"] <- 2L
  write.csv(wrong_rows, manifest_path, row.names = FALSE)
  testthat::expect_error(read_database(tables = "stocks"), "row counts differ.*stocks")
  writeBin(manifest_bytes, manifest_path)

  # Exercise missing assets after a successful read and same-size tampering.
  stopifnot(identical(latest, read_database()))
  released_stocks <- file.path(latest_dir, "stocks.csv.gz")
  stopifnot(file.rename(released_stocks, paste0(released_stocks, ".saved")))
  testthat::expect_error(read_database(), "stocks\\.csv\\.gz")
  stopifnot(file.rename(paste0(released_stocks, ".saved"), released_stocks))
  size <- file.size(released_stocks)
  stamp <- file.mtime(released_stocks)
  corrupted <- readBin(released_stocks, "raw", n = size)
  corrupted[[20L]] <- as.raw(bitwXor(as.integer(corrupted[[20L]]), 1L))
  writeBin(corrupted, released_stocks)
  Sys.setFileTime(released_stocks, stamp + 2)
  stopifnot(file.size(released_stocks) == size)
  testthat::expect_error(read_database(), "[Cc]hecksum mismatch.*stocks\\.csv\\.gz")
})
cat("SAIO release and legacy reads, tier dispatch, and cache provenance checks passed.\n")
