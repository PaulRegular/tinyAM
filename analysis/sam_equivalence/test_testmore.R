# Deterministic tests of the SAMjr-style runner, without fitting or downloads.
# Run from the repository root: Rscript analysis/sam_equivalence/test_testmore.R
source("analysis/sam_equivalence/samjr_testmore.R")
folder <- tempfile("samjr-testmore-"); dir.create(folder)
expected <- file.path(folder, "expected")
actual <- file.path(folder, "actual")
writeBin(charToRaw("value 1.000\n"), expected)
for (text in c("value 1.000\n", "value 1.000\r\n")) {
  writeBin(charToRaw(text), actual)
  stopifnot(samjr_output_match(actual, expected))
}
for (text in c("value 1.001\n", "value 1.000 \n", "value 1.000", " value 1.000\n")) {
  writeBin(charToRaw(text), actual)
  stopifnot(!samjr_output_match(actual, expected))
}
fits <- list(Standard = list(opt = list(objective = 10, par = numeric(2))),
  XtraSd = list(opt = list(objective = 9, par = numeric(3))))
table <- samjr_modeltable(fits)
stopifnot(identical(rownames(table), c("XtraSd", "Standard")),
  identical(unname(table[, "#par"]), c(3, 2)),
  identical(unname(table[, "AIC"]), c(24, 24)),
  is.na(table["XtraSd", 4L]), table["Standard", 4L] == 1 - pchisq(2, 1))
old <- options(digits = 3L, width = 120L)
samjr_write_output("codIN3Ben", list(list(opt = list(objective = 751.6052))), actual)
stopifnot(identical(getOption("digits"), 3L), identical(getOption("width"), 120L),
  identical(rawToChar(samjr_output_bytes(actual)), "751.6052 \n"))
options(old)
# The source workflow must use supplied starts, not call a real fitter or run
# downstream forecasts/outputs. It must restore the working directory on error.
writeLines(c("library(samjr)", "par <- list(logN = 123)",
  "fit <- sam.fit(list(obs = 456), list(age = 1), par)",
  "stop('Downstream expression must not run')"), file.path(folder, "script.R"))
cwd <- getwd()
models <- samjr_source_models("codIN3Ben", folder, function(data, conf, parameters) {
  list(data = data, conf = conf, start = parameters)
})
stopifnot(identical(models$baseline$start$logN, 123), identical(getwd(), cwd))
writeLines("stop('source failure')", file.path(folder, "script.R"))
error <- try(samjr_source_models("codIN3Ben", folder, identity), silent = TRUE)
stopifnot(inherits(error, "try-error"), identical(getwd(), cwd))
stopifnot(normalizePath(dirname(normalizePath(folder))) == normalizePath(tempdir()))
unlink(folder, recursive = TRUE)
cat("SAMjr-style runner: deterministic checks passed.\n")
