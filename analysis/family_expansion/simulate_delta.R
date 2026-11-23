# Run from the package root. All generated inputs and attempts are checkpoints.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "report"
if (!phase %in% c("audit", "exact", "spatial", "variable", "report")) {
  stop("Use simulate_delta.R {audit|spatial|variable|exact|report} [results directory].")
}
study <- "analysis/family_expansion"
out <- if (length(args) >= 2L) args[2] else file.path(study, "results", "default-z-fishing-share")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
for (d in c("draws", "fits", "gates")) dir.create(file.path(out, d), showWarnings = FALSE)

key <- function(experiment, driver, effort, id, age = NA_integer_, mortality = NA_character_,
                family = NA_character_) {
  data.frame(experiment, driver, effort, id, age, mortality, family)
}
checkpoint <- function(folder, name, identifier, evaluate) {
  path <- file.path(out, folder, paste0(name, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  warnings <- character()
  began <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(evaluate(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  ok <- !inherits(value, "error")
  result <- list(key = identifier, ok = ok, value = if (ok) value else NULL,
    error = if (ok) "" else conditionMessage(value), warnings = unique(warnings),
    elapsed = proc.time()[["elapsed"]] - began)
  temp <- paste0(path, ".partial")
  saveRDS(result, temp)
  if (!file.rename(temp, path)) stop("Could not checkpoint ", path)
  cat(sprintf("%s %.2fs %s\n", name, result$elapsed, if (ok) "saved" else result$error))
  flush.console()
  result
}
combine <- function(tables) as.data.frame(data.table::rbindlist(tables, fill = TRUE, use.names = TRUE))

summarize <- function() {
  files <- function(d) list.files(file.path(out, d), "[.]rds$", full.names = TRUE)
  gather <- function(folder, field) combine(lapply(files(folder), function(path) {
    x <- readRDS(path)
    result <- x$value[[field]]
    if (is.null(result) || !nrow(result)) return(NULL)
    if (field == "scores") names(result)[names(result) == "age"] <- "obs_age"
    identifiers <- x$key[setdiff(names(x$key), names(result))]
    cbind(identifiers[rep(1L, nrow(result)), , drop = FALSE], result)
  }))
  fits <- lapply(files("fits"), function(path) {
    x <- readRDS(path)
    cbind(x$key, converged = isTRUE(x$value$converged),
      error = x$error, warnings = paste(x$warnings, collapse = " | "), elapsed = x$elapsed,
      gradient = if (is.null(x$value$gradient)) NA_real_ else x$value$gradient,
      failed_checks = paste(x$value$checks$check[x$value$checks$status != "pass"], collapse = ", "))
  })
  available <- function(x) is.finite(x)
  score_by <- function(table, groups) {
    if (!nrow(table)) return(data.frame())
    dt <- data.table::as.data.table(table)
    out <- dt[, {
      pos <- holdout > 0
      z <- available(p_zero)
      ok <- pos & available(error)
      cvr <- pos & !is.na(covered)
      list(n = .N, n_brier = sum(z), n_positive = sum(pos), n_positive_scored = sum(ok),
        n_positive_coverage = sum(cvr), mean_zero = mean(p_zero[z]), observed_zero = mean(holdout == 0),
        brier = mean(((holdout == 0) - p_zero)^2, na.rm = TRUE),
        log_rmse = sqrt(mean(error[ok]^2)), positive_coverage = mean(covered[cvr]))
    }, by = groups]
    as.data.frame(out)
  }
  row_scores <- gather("fits", "scores")
  scores <- if (nrow(row_scores)) score_by(row_scores, c(names(key("", "", "", 0)), "component")) else data.frame()
  profile <- if (nrow(row_scores)) score_by(row_scores, c(names(key("", "", "", 0)),
    "component", "obs_age")) else data.frame()
  draws <- lapply(files("draws"), readRDS)
  reference <- sampling <- observations <- audits <- list()
  for (x in draws) {
    if (!x$ok || x$key$experiment == "calibration") next
    samples <- if (x$key$experiment == "spatial") x$value else setNames(list(x$value), paste(x$key$driver, x$key$effort, sep = "-"))
    for (name in names(samples)) {
      bundle <- samples[[name]]
      cell <- strsplit(name, "-", fixed = TRUE)[[1L]]
      k <- x$key[c("experiment", "id")]
      k$driver <- cell[1]
      k$effort <- cell[2]
      reference[[length(reference) + 1L]] <- cbind(k, bundle$reference)
      if (!is.null(bundle$audit)) {
        audits[[length(audits) + 1L]] <- cbind(k, bundle$audit[setdiff(names(bundle$audit), names(k))])
      }
      if (!is.null(bundle$sampling) && nrow(bundle$sampling)) {
        sampling[[length(sampling) + 1L]] <- cbind(k, bundle$sampling[setdiff(names(bundle$sampling), names(k))])
      }
      observations[[length(observations) + 1L]] <- combine(lapply(c("catch", "index"), function(component) {
        cbind(k, component, bundle$data[[component]][c("year", "age", "obs")])
      }))
    }
  }
  result <- list(fits = combine(fits), recovery = gather("fits", "recovery"),
    parameters = gather("fits", "parameters"), scores = scores, profile = profile,
    reference = unique(combine(reference)), sampling = combine(sampling),
    observations = combine(observations), audits = combine(audits), gates = combine(lapply(files("gates"), readRDS)))
  saveRDS(result, file.path(out, "summary.rds"))
  cat("Summary:", nrow(result$fits), "fit attempts\n")
}
if (phase == "report") {
  summarize()
  html <- rmarkdown::render(file.path(study, "delta_validation_report.Rmd"),
    params = list(result = file.path(out, "summary.rds")), output_dir = out, quiet = TRUE)
  cat("Report:", html, "\n")
  quit(save = "no")
}

pkgload::load_all(Sys.getenv("TINYAM_DELTA_PACKAGE_SOURCE", "."), quiet = TRUE, helpers = FALSE)
source(file.path(study, "delta_validation_fit.R"))
source(file.path(study, "matched-fit.R"))
data.table::setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
options(warn = 1)
spatial <- phase %in% c("audit", "spatial", "variable")
source(file.path(study, "delta_validation.R"))
if (spatial) {
  source("analysis/formula_expansion/observation_error_sim.R")
  source(file.path(study, "delta_fishery_engine.R"))
  source(file.path(study, "matched-spatial.R"))
}
if (phase == "exact") source(file.path(study, "delta_simple.R"))
if (spatial || phase == "exact") {
  sim <- normalizePath(Sys.getenv("TINYAM_DELTA_SIMSURVEY"), winslash = "/", mustWork = TRUE)
  revision <- suppressWarnings(system2("git", c("-C", shQuote(sim), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE))
  if (!delta_simsurvey_revision %in% c(revision, sub("^SimSurvey-", "", basename(sim)))) {
    stop("Set TINYAM_DELTA_SIMSURVEY to SimSurvey at ", delta_simsurvey_revision)
  }
  pkgload::load_all(sim, quiet = TRUE, helpers = FALSE)
  if (spatial) delta_check_allocate_lengths()
}
source_files <- c(list.files("R", "[.]R$", full.names = TRUE), "DESCRIPTION", "NAMESPACE",
  file.path(study, c("delta_validation_fit.R", "matched-fit.R", "simulate_delta.R")))
if (phase == "exact") source_files <- c(source_files,
  file.path(study, c("delta_simple.R", "delta_validation.R")),
  sort(list.files(file.path(sim, "R"), "[.]R$", full.names = TRUE)))
if (spatial) source_files <- c(source_files, "analysis/formula_expansion/observation_error_sim.R",
  file.path(study, c("matched-spatial.R", "delta_validation.R", "delta_fishery_engine.R")),
  sort(list.files(file.path(sim, "R"), "[.]R$", full.names = TRUE)))
manifest <- list(R = getRversion(), RTMB = utils::packageVersion("RTMB"),
  tweedie = utils::packageVersion("tweedie"), sources = tools::md5sum(source_files))
manifest_path <- file.path(out, paste0(phase, "-manifest.rds"))
if (file.exists(manifest_path) && !identical(readRDS(manifest_path), manifest)) {
  stop("Study sources or dependencies differ from the checkpoint manifest. Use a new results directory.")
}
if (!file.exists(manifest_path)) saveRDS(manifest, manifest_path)

drivers <- c("high_f", "low_f")
efforts <- c("low", "standard", "high")
cells <- expand.grid(driver = drivers, effort = efforts, stringsAsFactors = FALSE)
families <- names(matched_families)

phase_status <- function(experiment, status, reason, attempted = NA_integer_, common = NA_integer_, planned = 24L) {
  saveRDS(data.frame(experiment, status, reason, attempted, common_success = common, planned),
    file.path(out, "gates", paste0(experiment, ".rds")))
  cat(experiment, status, reason, "\n")
}

audit_sample <- function(sample) {
  d <- sample$reference
  stopifnot(nrow(d) == 30L * max(d$age), !anyDuplicated(d[c("year", "age")]),
    all(is.finite(d$M_realized)), mean(d$M_realized) > 0,
    max(abs(d$Z - d$F_realized - d$M_realized)) < 1e-10,
    isTRUE(all.equal(d$index_expected, d$N, tolerance = 1e-12)))
  catch <- with(d, N * F_realized * -expm1(-Z) / Z)
  stopifnot(max(abs(catch - d$catch_expected)) <= 1e-8 * max(1, catch))
  for (what in c("data", "holdout")) for (part in c("catch", "index")) {
    obs <- sample[[what]][[part]]
    stopifnot(all(is.na(obs$obs) | is.finite(obs$obs) & obs$obs >= 0),
      any(obs$obs > 0, na.rm = TRUE))
    if (part == "index") stopifnot(all(obs$samp_time == 0))
  }
}

audit_pair <- function(a, b) {
  truth <- c("year", "age", "N", "Z", "ssb")
  ordered <- function(x) x[order(x$year, x$age), truth]
  stopifnot(isTRUE(all.equal(ordered(a$reference), ordered(b$reference), tolerance = 1e-12,
    check.attributes = FALSE)), identical(a$data$index, b$data$index),
    identical(a$holdout$index, b$holdout$index),
    isTRUE(a$audit$valid), isTRUE(b$audit$valid))
  if (a$audit$realized_share <= .5 || b$audit$realized_share >= .5) {
    stop("Realized catches did not retain the intended high and low fishing shares.")
  }
}

fit_cell <- function(draw, experiment, driver, effort, id, age, mortality, arms) {
  records <- setNames(lapply(arms, function(family) {
    identifier <- key(experiment, driver, effort, id, age, mortality, family)
    name <- paste(c(unlist(identifier[1:4]), age, mortality, family), collapse = "-")
    checkpoint("fits", name, identifier, function() {
      if (!draw$ok) stop(draw$error)
      matched_fit(draw$value, age, family, experiment, mortality)
    })
  }), arms)
  scored <- lapply(Filter(function(x) isTRUE(x$value$converged), records), function(x) {
    x$value$scores[c("component", "year", "age", "holdout")]
  })
  if (length(scored) > 1L) {
    stopifnot(all(vapply(scored[-1L], identical, logical(1L), y = scored[[1L]])))
  }
  vapply(records, function(x) isTRUE(x$value$converged), logical(1L))
}

run_spatial <- function(id, fit = TRUE, mortality = "known") {
  draw <- checkpoint("draws", sprintf("spatial-%02d", id), key("spatial", "", "", id),
    function() delta_matched_spatial(id))
  if (draw$ok) {
    stopifnot(all(paste(cells$driver, cells$effort, sep = "-") %in% names(draw$value)))
    lapply(draw$value, audit_sample)
    for (effort in efforts) audit_pair(draw$value[[paste("high_f", effort, sep = "-")]],
      draw$value[[paste("low_f", effort, sep = "-")]])
  } else if (!fit) stop(draw$error)
  if (!fit) return(invisible(NULL))
  selected <- if (mortality == "estimated") subset(cells, effort == "standard" & driver == "low_f") else cells
  outcome <- lapply(seq_len(nrow(selected)), function(j) {
    with(selected[j, ], {
      bundle <- draw
      if (draw$ok) bundle$value <- draw$value[[paste(driver, effort, sep = "-")]]
      primary <- fit_cell(bundle, "spatial", driver, effort, id, 20L, mortality, families)
      if (effort == "standard" && id <= 8L && mortality == "known") {
        fit_cell(bundle, "spatial", driver, effort, id, 14L, mortality, families)
      }
      if (effort == "standard") matrix(primary, nrow = 1L, dimnames = list(driver, names(primary))) else NULL
    })
  })
  do.call(rbind, Filter(Negate(is.null), outcome))
}

run_paired <- function(experiment, run, audit_only = FALSE, low_f_only = FALSE) {
  audit <- tryCatch({ lapply(1:2, run, fit = FALSE); NULL }, error = conditionMessage)
  if (!is.null(audit)) {
    phase_status(experiment, "audit-failed", audit, 0L)
    return(invisible(NULL))
  }
  phase_status(experiment, "audited", "The first two paired histories passed observation and mortality checks.", 0L)
  if (audit_only) return(invisible(NULL))
  pilot <- lapply(1:8, run)
  pair_counts <- Reduce(`+`, lapply(pilot, function(results) results & results[, "lognormal"]))
  worst <- min(pair_counts)
  cat("Common-success counts at standard effort\n")
  print(pair_counts)
  if (low_f_only) {
    phase_status(experiment, "narrowed-complete", paste(
      "Completed eight low-fishing histories. Retained the high-fishing attempts already made.",
      "The original two-share criterion for expansion could not be met."), 8L, worst, planned = 8L)
  } else if (worst >= 7L) {
    phase_status(experiment, "expanded", "Every standard-effort comparison passed the seven-pair gate.", 8L, worst)
    lapply(9:24, run)
    phase_status(experiment, "complete", "Completed the planned expansion.", 24L, worst)
  } else {
    phase_status(experiment, "pilot-limited", "One or more standard-effort comparisons had fewer than seven successful pairs.", 8L, worst)
  }
}

if (phase %in% c("audit", "spatial")) run_paired("spatial", run_spatial, audit_only = phase == "audit")
if (phase == "variable") run_paired("variable-M", function(id, fit = TRUE) {
  run_spatial(id, fit, mortality = "estimated")
}, low_f_only = TRUE)
if (phase == "exact") {
  calibration <- checkpoint("draws", "calibration", key("calibration", "", "", 0L), delta_share_calibration)
  if (!calibration$ok) stop(calibration$error)
  for (scheme in c("poisson-link", "tweedie")) {
    experiment <- if (scheme == "poisson-link") "exact-poisson" else "exact-tweedie"
    for (id in 1:8) {
      identifier <- key(experiment, "low_f", "standard", id)
      draw <- checkpoint("draws", paste(experiment, id, sep = "-"), identifier,
        function() delta_share_exact(id, scheme, calibration$value))
      if (draw$ok) audit_sample(draw$value)
      fit_cell(draw, experiment, "low_f", "standard", id, 20L, "benchmark", families)
      if (scheme == "poisson-link") fit_cell(draw, experiment, "low_f", "standard", id, 20L,
        "ablation", c("poisson_ln", "positive_only"))
    }
    phase_status(experiment, "complete", "Completed the eight-history exact benchmark.", 8L, planned = 8L)
  }
}
summarize()
