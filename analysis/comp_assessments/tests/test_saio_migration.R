source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/update_database_cache.R")
registry <- assessment_translations()
stopifnot(all(file.exists(registry$script)),
  !anyDuplicated(paste(registry$assessment_id, registry$tier)),
  nrow(assessment_translations("curated")) == 38L)

local({
  original_wd <- getwd()
  fixture <- tempfile("saio_snapshot_")
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
    outputs = data.frame(assessment_id = assessments$assessment_id))
  for (name in names(tables)) write.csv(tables[[name]],
    file.path("data_repo/database", paste0(name, ".csv")), row.names = FALSE)
  git("init", "-b", "main")
  git("add", "database")
  commit <- function(text) git("-c", "user.name=SnapshotTest", "-c", "user.email=snapshot@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", text)
  commit("fixture")
  revision <- git("rev-parse", "HEAD")[[1]]
  update_database_cache("data_repo")
  pinned <- read_database()
  stopifnot(pinned$stocks$scientific_name == "Source species",
    pinned$commit == revision, pinned$database_selection == "cached_snapshot")
  # Reading remains local and stable when upstream data advance or Git is unavailable.
  stocks$scientific_name <- "New committed species"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  stopifnot(read_database()$stocks$scientific_name == "Source species")
  git("add", "database"); commit("new data")
  stopifnot(read_database()$commit == revision)
  update_database_cache("data_repo")
  current <- read_database()
  stopifnot(current$commit != revision, current$stocks$scientific_name == "New committed species",
    dir.exists(file.path(relative_root, "cache", "saio", revision)))
  # An incomplete source copy cannot change the selected manifest.
  writeLines("next source revision", "data_repo/note.txt")
  git("add", "note.txt"); commit("copy interruption fixture")
  file.rename("data_repo/database/outputs.csv", "data_repo/database/outputs.saved")
  selected <- readLines(file.path(relative_root, "saio_revisions.csv"))
  copy_error <- suppressWarnings(tryCatch(update_database_cache("data_repo"), error = identity))
  stopifnot(inherits(copy_error, "error"),
    identical(readLines(file.path(relative_root, "saio_revisions.csv")), selected),
    identical(read_database()$commit, current$commit))
  file.rename("data_repo/database/outputs.saved", "data_repo/database/outputs.csv")
  rename <- file.rename("data_repo", "hidden_source")
  original_head <- .database_head
  .database_head <- function(...) stop("Reader must not invoke Git")
  stopifnot(rename, identical(read_database()$commit, current$commit))
  .database_head <- original_head
  # A missing copy gives the manual-refresh instruction.
  path <- file.path(relative_root, "cache", "saio", current$commit, "database.rds")
  file.rename(path, paste0(path, ".saved"))
  bad <- tryCatch(read_database(), error = identity)
  stopifnot(inherits(bad, "error"), grepl("unavailable or incomplete", conditionMessage(bad)))
  file.rename(paste0(path, ".saved"), path)
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
cat("Manual SAIO snapshots, stable local reads, tier dispatch and source-keyed cache checks passed.\n")
