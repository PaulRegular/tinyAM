# Source retrieval is separate from the routine database-to-model workflow.
cache_pacific_source <- function(stock_id, url, name, revision = "") {
  directory <- file.path("analysis/comp_assessments/source_cache", stock_id)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  manifest_file <- file.path(directory, "manifest.csv")
  manifest <- if (file.exists(manifest_file)) {
    read.csv(manifest_file, stringsAsFactors = FALSE)
  } else {
    data.frame(file = character(), url = character(), retrieved = character(),
               md5 = character(), revision = character(), status = character())
  }
  existing <- manifest$file == name & manifest$url == url
  if (any(existing) && file.exists(path) &&
      identical(unname(tools::md5sum(path)), manifest$md5[which(existing)[1L]])) {
    return(path)
  }
  status <- tryCatch({
    utils::download.file(url, path, mode = "wb", quiet = TRUE)
    "downloaded"
  }, error = function(e) conditionMessage(e))
  checksum <- if (status == "downloaded") unname(tools::md5sum(path)) else ""
  manifest <- rbind(manifest[!existing, ], data.frame(
    file = name, url = url, retrieved = format(Sys.Date(), "%Y-%m-%d"),
    md5 = checksum, revision = revision, status = status
  ))
  utils::write.csv(manifest, manifest_file, row.names = FALSE, na = "")
  if (status != "downloaded") cli::cli_abort("Could not cache {url}: {status}")
  path
}
