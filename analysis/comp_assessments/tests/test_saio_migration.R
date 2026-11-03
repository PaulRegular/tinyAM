source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/run_assessment.R")

registry <- assessment_translations()
stopifnot(all(file.exists(registry$script)),
  !anyDuplicated(paste(registry$assessment_id, registry$tier)),
  nrow(assessment_translations("curated")) == 35L,
  !anyNA(.saio_read_table("stocks.csv")$stock_id),
  identical(.database_head(saio_path(), .saio_revision()), .saio_revision()))

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
    result <- system2("git", c("-C", shQuote("data_repo"), ...),
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
  pinned <- read_committed_database(revision, repo = "data_repo")
  working <- read_database("data_repo/database")
  stopifnot(pinned$stocks$scientific_name == "Source species",
    working$stocks$scientific_name == "Uncommitted replacement",
    identical(pinned$commit, revision))

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
    grepl("results/curated/cache/shared$", .assessment_cache("shared", "curated")),
    grepl("results/preliminary/cache/shared$", .assessment_cache("shared", "preliminary")))
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
  rejected <- tryCatch(run_assessments("extra", pinned,
    fit = FALSE, progress = FALSE), error = identity)
  mixed <- tryCatch(.assessment_result_tables(list(preferred, alternative),
    save_results = TRUE), error = identity)
  stopifnot(inherits(rejected, "error"), inherits(mixed, "error"))
  for (tier in c("curated", "preliminary")) {
    folder <- file.path(.assessment_results(tier), "dashboards")
    dir.create(folder, recursive = TRUE)
    writeLines("<html>Fixture</html>", file.path(folder, paste0(tier, "_only.html")))
    index_assessment_dashboards(tier)
  }
  curated_index <- readLines(file.path(.assessment_results("curated"), "index.html"))
  preliminary_index <- readLines(file.path(.assessment_results("preliminary"), "index.html"))
  stopifnot(any(grepl("curated_only.html", curated_index, fixed = TRUE)),
    !any(grepl("preliminary_only.html", curated_index, fixed = TRUE)),
    any(grepl("preliminary_only.html", preliminary_index, fixed = TRUE)),
    !any(grepl("curated_only.html", preliminary_index, fixed = TRUE)))
})
cat("Pinned SAIO reads, tier dispatch, curated defaults and separate output regression checks passed.\n")
