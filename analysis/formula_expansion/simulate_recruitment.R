args <- commandArgs(trailingOnly = TRUE)
n_simple <- if (length(args)) as.integer(args[1L]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2L]) else 30L
workers <- if (length(args) > 2L) as.integer(args[3L]) else 1L
if (any(!is.finite(c(n_simple, n_full, workers))) || any(c(n_simple, n_full, workers) < 1L)) {
  cli::cli_abort("Replicate counts and workers must be positive integers.")
}
pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
source(file.path(root, "recruitment_validation.R"))
dir.create(file.path(root, "results"), showWarnings = FALSE, recursive = TRUE)
path <- file.path(root, "results", "recruitment_recovery.rds")
signature <- tools::md5sum(c(file.path(root, "recruitment_validation.R"),
  file.path(root, "simulate_recruitment.R"), list.files("R", "[.]R$", full.names = TRUE)))
attempts <- list()
study_notes <- character()
if (file.exists(path)) {
  saved <- readRDS(path)
  if (!identical(signature, saved$signature) || saved$n_simple != n_simple || saved$n_full != n_full) {
    cli::cli_abort("The checkpoint uses different source or replicate counts; move it aside before a new study.")
  }
  attempts <- saved$results
  study_notes <- saved$study_notes
}
run <- function(stage, case, design, i) {
  seed <- 580000L + match(case, rec_cases) * 10000L + match(design, c("well informed", "narrow SSB", "correlated covariate")) * 1000L + i
  warnings <- character()
  began <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(
    if (stage == "full") rec_full_recovery(case, seed + 100000L) else rec_isolated_recovery(case, design, seed),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) list(error = conditionMessage(e)))
  diagnostics <- if (is.null(result$error)) result$diagnostics else data.frame(
    optimizer = NA_integer_, message = result$error, objective = NA_real_, gradient = NA_real_, pdHess = FALSE, success = FALSE)
  result$attempt <- cbind(data.frame(stage, case, design, replicate = i, seed = seed + if (stage == "full") 100000L else 0L),
    diagnostics, elapsed = proc.time()[["elapsed"]] - began,
    curve_rmse = if (is.null(result$curve_rmse)) NA_real_ else result$curve_rmse,
    warnings = paste(unique(warnings), collapse = " | "))
  if (!is.null(result$fit) && i == 1L) saveRDS(result$fit, file.path(root, "results", paste0("recruitment_", case, ".rds")))
  result$fit <- NULL
  result
}
key <- function(stage, case, design, i) paste(stage, case, design, i, sep = "_")
checkpoint <- function(complete = FALSE) saveRDS(list(results = attempts, signature = signature,
  n_simple = n_simple, n_full = n_full, workers = workers, complete = complete,
  study_notes = study_notes, session = sessionInfo()), path)
for (case in rec_cases) {
  designs <- if (grepl("^cov", case)) "well informed" else c("well informed", "narrow SSB", "correlated covariate")
  for (design in designs) {
    for (i in seq_len(n_simple)) {
      id <- key("isolated", case, design, i)
      if (!id %in% names(attempts)) attempts[[id]] <- run("isolated", case, design, i)
    }
    checkpoint()
    cli::cli_inform("Isolated recruitment: {case}, {design} complete.")
  }
}
if (workers == 1L) future::plan(future::sequential) else
  future::plan(future::multisession, workers = workers)
for (case in rec_cases) {
  pending <- which(!vapply(seq_len(n_full), function(i) key("full", case, "well informed", i) %in% names(attempts), logical(1)))
  for (batch in split(pending, ceiling(seq_along(pending) / workers))) {
    results <- furrr::future_map(batch, function(i) {
      if (workers > 1L) pkgload::load_all(quiet = TRUE)
      run("full", case, "well informed", i)
    }, .options = furrr::furrr_options(seed = TRUE))
    for (j in seq_along(batch)) attempts[[key("full", case, "well informed", batch[j])]] <- results[[j]]
    checkpoint()
    cli::cli_inform("Full recruitment: {case}, {max(batch)}/{n_full} recorded.")
  }
}
future::plan(future::sequential)
checkpoint(TRUE)
