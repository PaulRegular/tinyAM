source("analysis/comp_assessments/R/run_assessment.R")
registry <- assessment_translations()
stopifnot(all(file.exists(registry$script)),
  !anyDuplicated(paste(registry$assessment_id, registry$tier)),
  nrow(assessment_translations("curated")) == 40L)

local({
  original_wd <- getwd()
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
  pinned <- read_database(saio_repo = "data_repo")
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

  current <- read_database(saio_repo = "data_repo")
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
})
cat("Committed SAIO reads, tier dispatch and source-keyed cache checks passed.\n")
