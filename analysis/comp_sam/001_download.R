# Run from the tinyAM repository root. No stockassessment package is required.
revision <- "c6cfd035c7de59f7b3421dde31901efbca4cb0e8"
here <- file.path("analysis", "comp_sam")
target <- file.path(here, "source")
dir.create(target, recursive = TRUE, showWarnings = FALSE)
inputs <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", "mo.dat", "nm.dat",
            "pf.dat", "pm.dat", "sw.dat", "survey.dat", "script.R", "res.EXP")
paths <- c(paste0("testmore/nscod/", inputs), "stockassessment/tests/nscod/fit.expected.Rdata")
urls <- paste0("https://raw.githubusercontent.com/fishfollower/SAM/", revision, "/", paths)
dest <- file.path(target, basename(paths))
expected <- utils::read.csv(file.path(here, "source_manifest.csv"), stringsAsFactors = FALSE)
for (i in seq_along(paths)) {
  if (!file.exists(dest[i])) utils::download.file(urls[i], dest[i], mode = "wb", quiet = TRUE)
  checksum <- unname(tools::md5sum(dest[i]))
  if (!identical(checksum, expected$md5[match(paths[i], expected$source_path)])) {
    cli::cli_abort("Checksum mismatch for {paths[i]}. Inspect the cache before proceeding.")
  }
}
cat("Verified", length(paths), "public SAM files at", revision, "\n")
