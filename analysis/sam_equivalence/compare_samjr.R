# Run from the repository root with R 4.4 and the pinned SAM library:
# Rscript analysis/sam_equivalence/compare_samjr.R SAM-library [testmore|baseline|all]
args <- commandArgs(TRUE)
if (!length(args)) stop("Supply the pinned stockassessment library directory.")
.libPaths(c(normalizePath(args[1L]), .libPaths()))
library(stockassessment)
pkgload::load_all(".", quiet = TRUE)
source("analysis/sam_equivalence/samjr_helpers.R")
source("analysis/sam_equivalence/samjr_testmore.R")
revision <- "1cc464b80f6f2a4a38b35cce276f13e521761da3"
stopifnot(utils::packageDescription("stockassessment")$RemoteSha == revision)
dir.create("analysis/sam_equivalence/results", recursive = TRUE, showWarnings = FALSE)
root <- normalizePath("analysis/sam_equivalence/results", mustWork = TRUE)
writeLines(c(paste("SAM revision:", revision),
  "SAMjr fixtures: 0384774a3f2f8ec3dc3330375998bd73ef6bbe7c",
  capture.output(sessionInfo())), file.path(root, "sessionInfo.txt"))
selection <- if (length(args) > 1L) match.arg(args[2L], c("testmore", "baseline", "all")) else "all"
failed <- FALSE
if (selection %in% c("testmore", "all")) {
  for (case in c("codIN3Ben", "nscodXtraSd")) {
    folder <- file.path(root, "fixtures", case)
    samjr_fixture(case, folder)
    output <- file.path(root, "testmore", case)
    result <- tryCatch(samjr_testmore(case, folder, output), error = function(e) e)
    if (inherits(result, "error")) {
      writeLines(conditionMessage(result), file.path(output, "failure.txt"))
      cat(case, "FAILED:", conditionMessage(result), "\n")
      failed <- TRUE
    } else {
      unlink(file.path(output, "failure.txt"))
      print(result)
      failed <- failed || any(!result$output_match | !result$numerical_success)
    }
  }
}
if (selection == "testmore") quit(status = as.integer(failed))
cases <- if (selection == "baseline") "nscodbiopro" else c("nscodbiopro", "nscodXtraSd", "codIN3Ben")
for (case in cases) {
  folder <- file.path(root, "fixtures", case)
  samjr_fixture(case, folder)
  warnings <- character()
  native <- tryCatch(withCallingHandlers(
    samjr_source_models(case, folder, stockassessment::sam.fit),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))), error = function(e) e)
  out <- file.path(root, "samjr", case)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (inherits(native, "error")) {
    writeLines(c(warnings, conditionMessage(native)), file.path(out, "source-failure.txt"))
    failed <- TRUE
    next
  }
  if (case == "nscodXtraSd") names(native) <- c("standard", "extra_sd")
  saveRDS(native, file.path(out, "native.rds"))
  unlink(file.path(out, "source-failure.txt"))
  for (name in names(native)) {
    label <- paste(case, name, sep = "-")
    result <- tryCatch(withCallingHandlers(samjr_compare(native[[name]], label),
      warning = function(w) warnings <<- c(warnings, conditionMessage(w))), error = function(e) e)
    saveRDS(result, file.path(out, paste0(name, ".rds")))
    if (inherits(result, "error")) {
      writeLines(conditionMessage(result), file.path(out, paste0(name, "-failure.txt")))
      cat(label, "FAILED:", conditionMessage(result), "\n")
      failed <- TRUE
    } else {
      unlink(file.path(out, paste0(name, "-failure.txt")))
      write.csv(result$checks, file.path(out, paste0(name, "-checks.csv")), row.names = FALSE)
      print(result$checks)
      failed <- failed || any(!result$checks$passed)
    }
  }
  writeLines(unique(warnings), file.path(out, "warnings.txt"))
}
if (failed) quit(status = 1)
