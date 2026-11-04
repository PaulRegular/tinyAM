# Manual step: commit SAIO tables and resolve Git LFS files before running.
update_database_cache <- function(saio_repo,
    cache_dir = file.path("analysis", "comp_assessments", "cache", "saio"),
    manifest = file.path("analysis", "comp_assessments", "saio_revisions.csv")) {
  files <- c("stocks.csv", "assessments.csv", "assumptions.csv", "inputs.csv", "outputs.csv")
  revision <- system2("git", c("-C", shQuote(saio_repo), "rev-parse", "HEAD"),
    stdout = TRUE)
  if (length(revision) != 1L || !is.null(attr(revision, "status"))) {
    cli::cli_abort("Cannot read SAIO HEAD. Set saio_repo to the local committed SAIO checkout.")
  }
  target <- file.path(cache_dir, revision)
  if (!dir.exists(target)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    staging <- tempfile(".snapshot-", tmpdir = cache_dir)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    if (!all(file.copy(file.path(saio_repo, "database", files), file.path(staging, files)))) {
      cli::cli_abort("Could not copy all five SAIO tables; the previous manifest remains selected.")
    }
    if (!file.rename(staging, target)) cli::cli_abort("Could not publish the SAIO snapshot.")
  }
  # Parse the copied tables once, using the same base-R reader as before.
  parsed <- file.path(target, "database.rds")
  if (!file.exists(parsed)) {
    tables <- lapply(file.path(target, files), read.csv, stringsAsFactors = FALSE,
      na.strings = c("", "NA"), check.names = FALSE)
    names(tables) <- sub("[.]csv$", "", files)
    temporary <- tempfile(".database-", tmpdir = target)
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(tables, temporary)
    if (!file.rename(temporary, parsed)) cli::cli_abort("Could not save the parsed SAIO cache.")
  }
  # Select the revision only after its five CSVs have been copied successfully.
  temporary <- tempfile(".saio-manifest-", tmpdir = dirname(manifest))
  on.exit(unlink(temporary), add = TRUE)
  write.csv(data.frame(revision = revision), temporary, row.names = FALSE)
  if (!file.rename(temporary, manifest)) cli::cli_abort("Could not replace the SAIO snapshot manifest.")
  cli::cli_inform("SAIO snapshot selected: {revision}. Both translation tiers use these local CSVs.")
  invisible(revision)
}
