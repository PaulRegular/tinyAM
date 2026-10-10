args <- commandArgs(TRUE)
n_isolated <- if (length(args)) as.integer(args[1]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2]) else 30L
stopifnot(n_isolated > 0L, n_full > 0L)
pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
source(file.path(root, "q_validation.R"))
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
diagnostics <- parameters <- examples <- list()

checkpoint <- function(complete = FALSE) {
  saveRDS(list(diagnostics = do.call(rbind, diagnostics), parameters = do.call(rbind, parameters),
    examples = examples, n_isolated = n_isolated, n_full = n_full, complete = complete,
    session = sessionInfo()), file.path(root, "results", "q_recovery.rds"))
}

record <- function(result, scope, spec, link, design, replicate, seed) {
  key <- paste(scope, spec, link, design, replicate, sep = "_")
  metadata <- data.frame(scope, spec, link, design, replicate, seed)
  if (!is.null(result$error)) {
    row <- data.frame(optimizer = NA_integer_, message = result$error, gradient = NA_real_,
      pdHess = FALSE, converged = FALSE, q_rmse = NA_real_, q_coverage = NA_real_,
      observation_sd = NA_real_)
  } else {
    row <- result$diagnostics
    parameters[[key]] <<- cbind(metadata, result$parameters)
    if (replicate == 1L) examples[[key]] <<- cbind(metadata, result$curve)
    if (scope == "full" && design == "replicated" && replicate == 1L) {
      saveRDS(result$fit, file.path(root, "results", paste0("q_", key, ".rds")))
    }
  }
  diagnostics[[key]] <<- cbind(metadata, row, warnings = paste(result$warnings, collapse = " | "))
}

for (link in c("log", "logit")) for (spec in q_specs) {
  for (sparse in c(FALSE, TRUE)) for (i in seq_len(n_isolated)) {
    seed <- 210000L + match(link, c("log", "logit")) * 10000L +
      match(spec, q_specs) * 1000L + sparse * 100L + i
    result <- q_capture(q_isolated_recovery, spec, link, sparse, seed)
    record(result, "isolated", spec, link, if (sparse) "sparse" else "replicated", i, seed)
  }
  cli::cli_inform("Isolated q: {spec}, {link} complete.")
  for (i in seq_len(n_full)) {
    seed <- 310000L + match(link, c("log", "logit")) * 10000L + match(spec, q_specs) * 1000L + i
    result <- q_capture(q_full_recovery, spec, link, seed)
    record(result, "full", spec, link, "replicated", i, seed)
    cli::cli_inform("Full q: {spec}, {link}, replicate {i}/{n_full}.")
  }
  checkpoint()
}
for (link in c("log", "logit")) for (spec in c("iid", "ar1")) for (i in seq_len(n_full)) {
  seed <- 410000L + match(link, c("log", "logit")) * 10000L + match(spec, q_specs) * 1000L + i
  result <- q_capture(q_full_recovery, spec, link, seed, sparse = TRUE)
  record(result, "full", spec, link, "sparse", i, seed)
}
checkpoint(complete = TRUE)
