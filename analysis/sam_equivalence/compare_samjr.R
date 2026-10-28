# Run from the repository root with R 4.4 and the pinned SAM library:
# Rscript analysis/sam_equivalence/compare_samjr.R SAM-library [baseline|all]
args <- commandArgs(TRUE)
if (!length(args)) stop("Supply the pinned stockassessment library directory.")
.libPaths(c(normalizePath(args[1L]), .libPaths()))
library(stockassessment)
pkgload::load_all(".", quiet = TRUE)
source("analysis/sam_equivalence/samjr_helpers.R")
revision <- "1cc464b80f6f2a4a38b35cce276f13e521761da3"
stopifnot(utils::packageDescription("stockassessment")$RemoteSha == revision)
root <- normalizePath("analysis/sam_equivalence/results", mustWork = TRUE)
selection <- if (length(args) > 1L) match.arg(args[2L], c("baseline", "all")) else "all"
cases <- if (selection == "baseline") "nscodbiopro" else c("nscodbiopro", "nscodXtraSd", "codIN3Ben")
failed <- FALSE
for (case in cases) {
  folder <- file.path(root, "fixtures", case)
  samjr_fixture(case, folder)
  expressions <- parse(file.path(folder, "script.R"))
  env <- new.env(parent = globalenv())
  old <- setwd(folder)
  warnings <- character()
  native <- tryCatch(withCallingHandlers({
    for (expression in expressions) {
      if (identical(expression, quote(library(samjr)))) next
      eval(expression, env)
      if (exists(if (case == "nscodbiopro") "fit0" else if (case == "nscodXtraSd") "fit2" else "fit", env, inherits = FALSE)) break
    }
    if (case == "nscodbiopro") list(baseline = env$fit0) else if (case == "nscodXtraSd") list(standard = env$fit, extra_sd = env$fit2) else list(baseline = env$fit)
  }, warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) e, finally = setwd(old))
  out <- file.path(root, "samjr", case)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (inherits(native, "error")) {
    writeLines(c(warnings, conditionMessage(native)), file.path(out, "source-failure.txt"))
    failed <- TRUE
    next
  }
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
