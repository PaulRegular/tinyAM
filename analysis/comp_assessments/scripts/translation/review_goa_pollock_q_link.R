pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

# Fit identical observations and settings with each catchability link.
id <- "afsc_pollock_goa_2024"
database <- read_database()
source_data <- read_assessment(id, database)
stock <- new.env(parent = globalenv())
sys.source(file.path(.assessment_root, "scripts/translation/stocks", paste0(id, ".R")), stock)
translated <- stock$translate_stock(source_data)
cache <- file.path(.assessment_root, "results/cache", id)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

fits <- lapply(setNames(c("log", "logit"), c("log", "logit")), function(link) {
  settings <- translated$settings
  settings$index_settings$q_link <- link
  fitted <- do.call(tinyAM::fit_tam, c(
    list(obs = translated$obs, years = translated$years, ages = translated$ages,
         silent = TRUE), settings
  ))
  .assessment_catch_reporting(fitted, translated$catch_reporting)
})
for (link in names(fits)) saveRDS(fits[[link]], file.path(cache, paste0("fit_q_", link, ".rds")))
diagnostics <- do.call(rbind, lapply(names(fits), function(link) {
  fitted <- fits[[link]]
  out <- .assessment_diagnostics(id, database, link, fitted)
  terminal <- fitted$pop$ssb[fitted$pop$ssb$year == max(translated$years), ]
  out$terminal_ssb_tonnes <- terminal$est / 1000
  out$terminal_ssb_log_se <- terminal$se
  out
}))
write.csv(diagnostics, file.path(cache, "q_link_review.csv"), row.names = FALSE)
print(diagnostics)

# Compare both links with the accepted outputs; preserve the main stock recipe.
ref <- database_to_tam_ref(id, database$outputs, obs = translated$obs,
  years = translated$years, ages = translated$ages,
  terminal_year = source_data$assessment$terminal_year,
  age_plus_group = translated$age_plus_group,
  comparison_scales = translated$comparison_scales, template = fits$logit)
background <- c(translated$background, "",
  "Catchability-link sensitivity: the logit link restricts the full q prediction to between zero and one. Observations and all other settings are identical to the log-link fit. This restriction supplies no prior abundance scale.")
tinyAM::vis_tam(model_list = list(Accepted = ref, `tinyAM log` = fits$log,
                                 `tinyAM logit` = fits$logit),
  background = background, output_file = file.path(cache, "dashboard_logit_q.html"),
  open_file = interactive())
