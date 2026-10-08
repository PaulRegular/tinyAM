pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

# Compare the two retained changes separately and together.
id <- "afsc_pollock_goa_2024"
database <- read_database()
source_data <- read_assessment(id, database)
stock <- new.env(parent = globalenv())
sys.source(file.path(.assessment_root, "scripts/translation/stocks", paste0(id, ".R")), stock)
translated <- stock$translate_stock(source_data)
previous <- translated$settings
previous$N_settings <- list(process = "iid", init = "exp")
previous$index_settings$q_link <- "log"
models <- list(Previous = previous, `Logit only` = previous,
               `N off only` = previous, Revised = translated$settings)
models$`Logit only`$index_settings$q_link <- "logit"
models$`N off only`$N_settings$process <- "off"
cache <- file.path(.assessment_root, "results/cache", id)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

fits <- references <- diagnostics <- summaries <- list()
for (name in names(models)) {
  elapsed <- system.time(fitted <- do.call(tinyAM::fit_tam, c(
    list(obs = translated$obs, years = translated$years, ages = translated$ages,
         silent = TRUE), models[[name]]
  )))[["elapsed"]]
  fitted <- .assessment_catch_reporting(fitted, translated$catch_reporting)
  fits[[name]] <- fitted
  diagnostics[[name]] <- .assessment_diagnostics(id, database, name, fitted, elapsed)
  if (!isTRUE(fitted$is_converged)) next
  references[[name]] <- database_to_tam_ref(id, database$outputs, obs = translated$obs,
    years = translated$years, ages = translated$ages,
    terminal_year = source_data$assessment$terminal_year,
    age_plus_group = translated$age_plus_group,
    comparison_scales = translated$comparison_scales, template = fitted)
  differences <- .assessment_percent_differences(
    fitted, references[[name]], assumptions = source_data$assumptions)
  summaries[[name]] <- .assessment_comparison_summary(differences, id)
  summaries[[name]]$model <- name
}
diagnostics <- do.call(rbind, diagnostics)
summaries <- do.call(rbind, summaries)
write.csv(diagnostics, file.path(cache, "model_review_diagnostics.csv"), row.names = FALSE)
write.csv(summaries, file.path(cache, "model_review_summary.csv"), row.names = FALSE)
saveRDS(fits, file.path(cache, "model_review_fits.rds"))
print(diagnostics)
print(summaries[summaries$metric %in% c("N", "recruitment", "ssb"),
  c("model", "metric", "mean_absolute_percent_difference", "trend_correlation")])

# Keep the previous recipe visible alongside the revised fit.
if (isTRUE(fits$Previous$is_converged) && isTRUE(fits$Revised$is_converged)) {
  tinyAM::vis_tam(
    model_list = list(Accepted = references$Revised,
                      `tinyAM previous` = fits$Previous, `tinyAM revised` = fits$Revised),
    background = c(translated$background, "",
      "The previous tinyAM recipe uses IID older-age N deviations and a log q link. The revised recipe uses deterministic older-age survival and a logit q link; all observations, fixed M, F settings and catch SD are unchanged. Native SSB definitions differ as described above; numerical review summaries compare accepted N with the same weights and maturity as tinyAM."),
    output_file = file.path(cache, "dashboard_model_review.html"),
    open_file = interactive()
  )
}
