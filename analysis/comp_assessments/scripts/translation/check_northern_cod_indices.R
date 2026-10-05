root <- "analysis/comp_assessments"
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))
database <- read_committed_database()
accepted <- read_assessment("dfo_cod_2j3kl_2025", database)
stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "dfo_cod_2j3kl_2025.R"), stock)

trials <- c("RV", "RV_Smith", "RV_Smith_juveniles")
diagnostics <- lapply(trials, function(trial) {
  translated <- stock$translate_stock(
    accepted, smith_sound = trial != "RV", juveniles = trial == "RV_Smith_juveniles"
  )
  started <- Sys.time()
  fitted <- tryCatch(do.call(tinyAM::fit_tam, c(
    list(obs = translated$obs, years = translated$years, ages = translated$ages,
         silent = TRUE), translated$settings
  )), error = identity)
  failed <- inherits(fitted, "error")
  result <- .assessment_diagnostics(
    accepted$assessment$assessment_id, database,
    if (failed) "fit_failed" else if (isTRUE(fitted$is_converged)) "converged" else "not_converged",
    fit = if (failed) NULL else fitted,
    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")),
    reason = if (failed) conditionMessage(fitted) else ""
  )
  result$reason <- gsub("[[:space:]]+", " ", trimws(result$reason))
  result$trial <- trial
  result$selected <- trial == "RV_Smith"
  print(result)
  result
})
write.csv(do.call(rbind, diagnostics),
           file.path(root, "results", "northern_cod_indices.csv"), row.names = FALSE, na = "")
