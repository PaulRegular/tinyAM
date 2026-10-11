# Source-supported alternatives; stock recipes remain unchanged ----
pkgload::load_all(quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
root <- "analysis/formula_expansion/results/process_expansion/assessments"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
database <- read_committed_database()
attempts <- list()

read_stock <- function(id) {
  e <- list2env(list(source = read_assessment(id, database), do_fit = FALSE,
    silent = TRUE), parent = environment())
  sys.source(file.path("analysis/comp_assessments/scripts/translation/stocks",
    paste0(id, ".R")), envir = e)
  as.list(e)
}

capture_trial <- function(fun) {
  warnings <- character()
  out <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  out$warnings <- warnings
  out
}

reference_for <- function(stock, fitted, ssb_definition = NULL) {
  definitions <- stock$comparison_definitions
  if (!is.null(ssb_definition)) definitions$ssb <- ssb_definition
  database_to_tam_ref(stock$source$assessment$assessment_id,
    stock$comparison_outputs %||% stock$source$outputs,
    obs = stock$obs, years = stock$years, ages = stock$ages,
    terminal_year = stock$source$assessment$terminal_year,
    age_plus_group = stock$age_plus_group,
    comparison_scales = stock$comparison_scales,
    template = fitted, assumptions = stock$source$assumptions,
    comparison_aggregates = stock$comparison_aggregates,
    comparison_age_groups = stock$comparison_age_groups,
    comparison_definitions = definitions)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

for (id in c("afsc_pollock_ebs_2024", "ices_haddock_north_sea_2026", "ices_plaice_north_sea_2026")) {
  cli::cli_inform("Baseline: {id}")
  stock <- read_stock(id)
  baseline <- capture_trial(function() run_assessment(id, database, cache = FALSE))
  if (!is.null(baseline$ref)) baseline$comparison_summary <-
    .assessment_comparison_summary(.assessment_percent_differences(baseline$ref), id)
  attempts[[paste(id, "baseline", sep = "/")]] <- baseline
  saveRDS(attempts, file.path(root, "trials.rds"))
  if (is.null(baseline$fit)) next
  variants <- if (id == "afsc_pollock_ebs_2024") {
    list(spawning_time = list(ssb_settings = list(spawn_time = .25)))
  } else if (id == "ices_haddock_north_sea_2026") {
    stock$obs$weight$N_sd_group <- factor(ifelse(stock$obs$weight$age == 8, "plus", "older"))
    list(correlated_F = list(F_settings = list(process = "cor_rw")),
      grouped_N_SD = list(N_settings = list(process = "iid", init = "free", sd_form = ~ N_sd_group)))
  } else {
    # Published keyVarF: 0 1 2 2 2 2 2 3 3 3 (ages 1:10).
    keys <- c(0, 1, 2, 2, 2, 2, 2, 3, 3, 3)
    stock$obs$weight$F_sd_group <- factor(keys[match(stock$obs$weight$age, 1:10)])
    list(correlated_F = list(F_settings = list(process = "cor_rw", mean_ages = 2:6)),
      grouped_F_SD = list(F_settings = list(process = "rw", mean_ages = 2:6, sd_form = ~ F_sd_group)))
  }
  for (variant in names(variants)) {
    cli::cli_inform("Trial: {id}, {variant}")
    attempts[[paste(id, variant, sep = "/")]] <- capture_trial(function() {
      fit <- do.call(update, c(list(object = baseline$fit, data = stock$obs, silent = TRUE,
        start_par = as.list(baseline$fit$sdrep, "Estimate")), variants[[variant]]))
      ref <- reference_for(stock, fit, if (variant == "spawning_time") list(status = "approximate",
        definition = "Female SSB at year fraction 0.25 in both models",
        reason = "Timing matches the pinned source spawning month. Biological weight and maturity conventions still differ.") else NULL)
      list(fit = fit, ref = ref, comparison_summary = .assessment_comparison_summary(.assessment_percent_differences(ref), id))
    })
    saveRDS(attempts, file.path(root, "trials.rds"))
  }
  if (id != "afsc_pollock_ebs_2024" && all(vapply(names(variants), function(variant) {
    isTRUE(attempts[[paste(id, variant, sep = "/")]]$fit$is_converged)
  }, logical(1)))) {
    second <- variants[[2L]]
    if (id == "ices_haddock_north_sea_2026") second$F_settings <- variants[[1L]]$F_settings else
      second$F_settings$process <- "cor_rw"
    attempts[[paste(id, "combined", sep = "/")]] <- capture_trial(function() {
      fit <- do.call(update, c(list(object = baseline$fit, data = stock$obs, silent = TRUE,
        start_par = as.list(baseline$fit$sdrep, "Estimate")), second))
      ref <- reference_for(stock, fit)
      list(fit = fit, ref = ref, comparison_summary = .assessment_comparison_summary(.assessment_percent_differences(ref), id))
    })
    saveRDS(attempts, file.path(root, "trials.rds"))
  }
  available <- attempts[startsWith(names(attempts), paste0(id, "/"))]
  fits <- c(list(Accepted = baseline$ref), lapply(available, `[[`, "fit"))
  fits <- Filter(Negate(is.null), fits)
  names(fits) <- sub(paste0(id, "/"), "", names(fits), fixed = TRUE)
  rendered <- capture_trial(function() {
    vis_tam(model_list = fits, background = stock$background,
      output_file = file.path(root, paste0(id, ".html")), open_file = FALSE)
    list(success = TRUE)
  })
  saveRDS(rendered, file.path(root, paste0(id, "_render.rds")))
}
saveRDS(attempts, file.path(root, "trials.rds"))
summary <- do.call(rbind, lapply(names(attempts), function(key) {
  a <- attempts[[key]]
  fit <- a$fit
  data.frame(case = key, numerical = !is.null(fit) && fit$is_converged,
    optimizer = if (is.null(fit)) NA else fit$opt$convergence,
    pdHess = !is.null(fit$sdrep) && isTRUE(fit$sdrep$pdHess),
    gradient = if (is.null(fit)) NA else max(abs(fit$gradient)),
    error = a$error %||% "", warnings = paste(a$warnings, collapse = " | "))
}))
utils::write.csv(summary, file.path(root, "attempts.csv"), row.names = FALSE)
