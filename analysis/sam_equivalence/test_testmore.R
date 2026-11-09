# Deterministic tests of the SAMjr-style runner, without fitting or downloads.
# Run from the repository root: Rscript analysis/sam_equivalence/test_testmore.R
source("analysis/sam_equivalence/samjr_testmore.R")
source("analysis/sam_equivalence/samjr_helpers.R")
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
writeLines(c("log(L) #par AIC", "ESTMATED -276.9702 18 589.9404",
  "ASSIGNEDcor -276.9702 17 587.9404"), expected)
samjr_write_output("nscodcovar", list(ESTMATED = list(opt = list(objective = 276.9702,
  par = numeric(18)))), actual)
stopifnot(samjr_estimated_covariance_match(actual, expected), !samjr_output_match(actual, expected))
for (logit in c(-2, 0, 1.5)) {
  h <- 1e-4
  y <- samjr_rho_distance(logit + c(-h, 0, h))
  derivative <- samjr_rho_distance_derivatives(logit)
  stopifnot(abs(derivative[["first"]] - (y[3L] - y[1L]) / (2 * h)) < 1e-7,
    abs(derivative[["second"]] - (y[1L] - 2 * y[2L] + y[3L]) / h^2) < 1e-5)
}
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
