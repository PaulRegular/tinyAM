args <- commandArgs(trailingOnly = TRUE)
n_simple <- if (length(args)) as.integer(args[1]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2]) else 30L
workers <- if (length(args) > 2L) as.integer(args[3]) else 1L
stopifnot(n_simple > 0, n_full > 0, workers > 0)
pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
source(file.path(root, "sd_validation.R"))
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
path <- file.path(root, "results", "sd_effects_recovery.rds")
source_signature <- tools::md5sum(c(file.path(root, "sd_validation.R"),
  file.path(root, "simulate_sd_effects.R"), list.files("R", pattern = "[.]R$", full.names = TRUE)))
key <- function(stage, type, design, replicate) paste(stage, type, design, replicate, sep = "_")

attempts <- parameters <- curves <- list()
study_notes <- character()
if (file.exists(path)) {
  saved <- readRDS(path)
  stopifnot(saved$n_simple == n_simple, saved$n_full == n_full)
  if (!identical(saved$source_signature, source_signature)) {
    cli::cli_abort("The checkpoint uses different source code. Move it aside before starting a new study.")
  }
  attempts <- split(saved$attempts, with(saved$attempts, key(stage, type, design, replicate)))
  parameters <- split(saved$parameters, with(saved$parameters, key(stage, type, design, replicate)))
  curves <- saved$curves
  study_notes <- saved$study_notes
}
run <- function(stage, process, design, replicate, seed, fun, ...) {
  warnings <- character()
  began <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(fun(..., seed = seed), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  id <- key(stage, process, design, replicate)
  metadata <- data.frame(stage, type = process, design, replicate, seed)
  diagnostics <- if (!is.null(result$error)) data.frame(optimizer = NA_integer_,
    message = result$error, objective = NA_real_, gradient = NA_real_, pdHess = FALSE,
    success = FALSE) else result$diagnostics
  attempt <- cbind(metadata, diagnostics,
    elapsed = proc.time()[["elapsed"]] - began,
    log_sd_rmse = if (is.null(result$error)) result$log_sd_rmse else NA_real_,
    warnings = paste(unique(warnings), collapse = " | "))
  table <- curve <- NULL
  if (is.null(result$error)) {
    if (nrow(result$parameters)) table <- cbind(metadata, result$parameters)
    if (replicate == 1L) {
      curve <- cbind(metadata, result$curve)
      if (stage == "full") saveRDS(result$fit, file.path(root, "results", paste0("sd_", id, ".rds")))
    }
  }
  list(key = id, attempt = attempt, parameters = table, curve = curve)
}
record <- function(result) {
  attempts[[result$key]] <<- result$attempt
  parameters[[result$key]] <<- result$parameters
  curves[[result$key]] <<- result$curve
}
checkpoint <- function(complete = FALSE) {
  saveRDS(list(attempts = do.call(rbind, attempts), parameters = do.call(rbind, parameters),
    curves = curves, n_simple = n_simple, n_full = n_full, complete = complete,
    workers = workers, source_signature = source_signature, study_notes = study_notes,
    session = sessionInfo()), path)
}
for (type in c("iid", "rw", "ar1")) for (repetitions in c(30L, 3L, 1L)) {
  design <- paste(repetitions, "observations per age")
  for (i in seq_len(n_simple)) if (!key("isolated", type, design, i) %in% names(attempts)) {
    record(run("isolated", type, design, i,
      240000 + match(type, c("iid", "rw", "ar1")) * 10000 + repetitions * 100 + i,
      sd_simple_recovery, type = type, repetitions = repetitions))
  }
  checkpoint()
  message("Isolated: ", type, "; observations per age = ", repetitions)
}
for (type in c("common", "quadratic", "iid", "rw", "ar1")) {
  for (i in seq_len(n_simple)) if (!key("tails", type, "30 observations per age", i) %in% names(attempts)) {
    record(run("tails", type, "30 observations per age", i, 380000 + i, sd_tail_recovery, type = type))
  }
  checkpoint()
  message("Noisier tails: ", type)
}
if (workers == 1L) future::plan(future::sequential) else
  future::plan(future::multisession, workers = workers)
for (component in c("catch", "index")) for (type in c("iid", "rw", "ar1")) {
  pending <- which(!key("full", type, component, seq_len(n_full)) %in% names(attempts))
  for (batch in split(pending, ceiling(seq_along(pending) / workers))) {
    results <- furrr::future_map(batch, function(i) {
      if (workers > 1L) pkgload::load_all(quiet = TRUE)
      run("full", type, component, i,
        440000 + match(type, c("iid", "rw", "ar1")) * 10000 + (component == "index") * 1000 + i,
        sd_full_recovery, type = type, component = component)
    }, .options = furrr::furrr_options(seed = TRUE))
    lapply(results, record)
    checkpoint()
    message("Assessment: ", component, "; ", type, "; recorded ", max(batch), "/", n_full)
  }
}
future::plan(future::sequential)
checkpoint(TRUE)
