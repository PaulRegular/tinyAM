source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/run_assessment.R")

registry <- assessment_translations()
stopifnot(all(file.exists(registry$script)),
  !anyDuplicated(paste(registry$assessment_id, registry$tier)),
  nrow(assessment_translations("curated")) == 35L,
  !anyNA(.saio_read_table("stocks.csv")$stock_id),
  identical(.database_head(saio_path(), "refs/heads/main"), .saio_revision()),
  identical(.saio_revision("curated"), .saio_revision("preliminary")))

local({
  original_wd <- getwd()
  fixture <- tempfile("saio_migration_")
  dir.create(fixture)
  stopifnot(startsWith(normalizePath(fixture, winslash = "/"),
    normalizePath(tempdir(), winslash = "/")))
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
  assumptions <- data.frame(assessment_id = assessments$assessment_id,
    component = "population", setting = "age_range", value = "1-2",
    source_reference = "fixture", notes = "fixture")
  tables <- list(stocks = stocks, assessments = assessments,
    assumptions = assumptions,
    inputs = data.frame(assessment_id = assessments$assessment_id),
    outputs = data.frame(assessment_id = assessments$assessment_id))
  for (name in names(tables)) write.csv(tables[[name]],
    file.path("data_repo/database", paste0(name, ".csv")), row.names = FALSE)
  git("init", "-b", "main")
  git("add", "database")
  git("-c", "user.name=MigrationTest", "-c", "user.email=migration@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", "fixture")
  revision <- git("rev-parse", "HEAD")[[1]]
  stocks$scientific_name <- "Uncommitted replacement"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  source(file.path(relative_root, "R", "read_database.R"))
  original_options <- options(tinyAM.saio_path = normalizePath("data_repo", winslash = "/"))
  on.exit(options(original_options), add = TRUE)
  pinned <- read_committed_database(revision, repo = "data_repo")
  latest <- read_database()
  latest_preliminary <- read_database(tier = "preliminary")
  working <- read_database("data_repo/database")
  stopifnot(pinned$stocks$scientific_name == "Source species",
    latest$stocks$scientific_name == "Source species",
    identical(latest$commit, latest_preliminary$commit),
    latest$database_selection == "latest_main",
    latest$database_ref == "refs/heads/main",
    pinned$database_selection == "explicit_revision",
    working$stocks$scientific_name == "Uncommitted replacement",
    grepl("^working-[[:xdigit:]]{32}$", working$database_key),
    identical(read_database("data_repo/database")$database_key, working$database_key),
    identical(pinned$commit, revision))
  stocks$scientific_name <- "Another uncommitted value"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  stopifnot(read_database("data_repo/database")$database_key != working$database_key,
    read_database()$stocks$scientific_name == "Source species")
  stocks$scientific_name <- "Uncommitted replacement"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  git("switch", "-c", "experiment")
  git("add", "database")
  git("-c", "user.name=MigrationTest", "-c", "user.email=migration@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", "experiment")
  experimental <- git("rev-parse", "HEAD")[[1]]
  stopifnot(read_database()$stocks$scientific_name == "Source species",
    read_database(revision = experimental)$stocks$scientific_name == "Uncommitted replacement")
  git("switch", "main")
  git("update-ref", "refs/remotes/origin/main", revision)
  stocks$scientific_name <- "Latest committed main"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  git("add", "database")
  git("-c", "user.name=MigrationTest", "-c", "user.email=migration@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", "advance main")
  current <- read_database()
  stopifnot(current$stocks$scientific_name == "Latest committed main",
    current$commit != revision, current$database_key == current$commit,
    read_database(revision = revision)$stocks$scientific_name == "Source species")
  git("update-ref", "refs/remotes/origin/main", experimental)
  divergent <- tryCatch(read_database(), error = identity)
  stopifnot(inherits(divergent, "error"), grepl("behind or diverged", conditionMessage(divergent), fixed = TRUE),
    read_database(revision = revision)$stocks$scientific_name == "Source species")
  git("update-ref", "refs/remotes/origin/main", revision)
  git("switch", "-c", "remote-ahead")
  stocks$scientific_name <- "Remote main ahead"
  write.csv(stocks, "data_repo/database/stocks.csv", row.names = FALSE)
  git("add", "database")
  git("-c", "user.name=MigrationTest", "-c", "user.email=migration@example.invalid",
    "-c", "commit.gpgsign=false", "commit", "-m", "remote main")
  remote_ahead <- git("rev-parse", "HEAD")[[1]]
  git("switch", "main")
  git("update-ref", "refs/remotes/origin/main", remote_ahead)
  behind <- tryCatch(read_database(), error = identity)
  stopifnot(inherits(behind, "error"), grepl("behind or diverged", conditionMessage(behind), fixed = TRUE))
  git("update-ref", "refs/remotes/origin/main", revision)

  entries <- data.frame(assessment_id = c("shared", "shared", "extra"),
    tier = c("curated", "preliminary", "preliminary"), source_commit = revision)
  write.csv(entries, file.path(relative_root, "translations.csv"), row.names = FALSE)
  for (tier in c("curated", "preliminary")) {
    folder <- file.path(relative_root, "scripts", "translation", tier)
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
    .assessment_cache_key(working) != .assessment_cache_key(current))
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
    index_assessment_dashboards(tier, current)
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
cat("Latest committed SAIO reads, explicit snapshots, source-keyed caches, tier dispatch and provenance regression checks passed.\n")
