args <- commandArgs(trailingOnly = TRUE)
n_simple <- if (length(args)) as.integer(args[1L]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2L]) else 30L
workers <- if (length(args) > 2L) as.integer(args[3L]) else 1L
if (any(!is.finite(c(n_simple, n_full, workers))) || any(c(n_simple, n_full, workers) < 1L)) {
  cli::cli_abort("Replicate counts and workers must be positive integers.")
}
pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
source(file.path(root, "recruitment_cross_validation.R"))
signature <- tools::md5sum(c(file.path(root, c("recruitment_validation.R", "recruitment_cross_validation.R", "simulate_recruitment_noise.R")),
  list.files("R", "[.]R$", full.names = TRUE)))
path <- file.path(root, "results", "recruitment_noise_recovery.rds")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
attempts <- list()
if (file.exists(path)) {
  saved <- readRDS(path)
  if (!identical(signature, saved$signature) || saved$n_simple != n_simple || saved$n_full != n_full) {
    cli::cli_abort("The checkpoint uses different source or replicate counts; move it aside before a new study.")
  }
  attempts <- saved$results
}
checkpoint <- function(complete = FALSE) saveRDS(list(results = attempts, signature = signature,
  n_simple = n_simple, n_full = n_full, complete = complete, session = sessionInfo()), path)
if (workers == 1L) future::plan(future::sequential) else future::plan(future::multisession, workers = workers)
for (stage in c("isolated", "full")) for (contrast in if (stage == "isolated") c("narrow", "wider") else c("natural", "wider")) {
  for (sigma in c(.15, .35, .6)) {
    n <- if (stage == "isolated") n_simple else n_full
    pending <- which(!paste(stage, contrast, sigma, seq_len(n), sep = "_") %in% names(attempts))
    for (batch in split(pending, ceiling(seq_along(pending) / workers))) {
      results <- furrr::future_map(batch, function(i) {
        if (workers > 1L) pkgload::load_all(quiet = TRUE)
        seed <- 970000L + match(sigma, c(.15, .35, .6)) * 1000L + i + if (stage == "full") 100000L else 0L
        result <- rec_cross_dataset(stage, "bh_iid", i, sigma = sigma, contrast = contrast,
          fitted_cases = "bh_iid", seed = seed)
        result$results$bh_iid$attempt$sigma <- sigma
        result$results$bh_iid$attempt$contrast <- contrast
        result
      }, .options = furrr::furrr_options(seed = TRUE))
      for (j in seq_along(batch)) attempts[[paste(stage, contrast, sigma, batch[j], sep = "_")]] <- results[[j]]
      checkpoint()
      cli::cli_inform("Recruitment noise: {stage}, {contrast}, sigma {sigma}, {max(batch)}/{n} recorded.")
    }
  }
}
future::plan(future::sequential)
checkpoint(TRUE)
