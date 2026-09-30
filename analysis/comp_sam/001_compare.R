library(tinyAM)

# Inputs ----
case <- "WKCOD_combined_99"
analysis_dir <- file.path("analysis", "comp_sam")
results <- file.path(analysis_dir, "results")
reference <- file.path(analysis_dir, "cache", paste0(case, ".rds"))
dir.create(dirname(reference), recursive = TRUE, showWarnings = FALSE)
dir.create(results, showWarnings = FALSE)

if (!file.exists(reference)) {
  sam_fit <- stockassessment::fitfromweb(case, character.only = TRUE)
  saveRDS(sam_fit, reference)
  write.csv(data.frame(model = case,
    url = paste0("https://stockassessment.org/datadisk/stockassessment/userdirs/user3/", case, "/run/model.RData"),
    downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stockassessment_version = as.character(packageVersion("stockassessment")),
    stockassessment_sha = packageDescription("stockassessment")$RemoteSha),
    file.path(analysis_dir, "reference_provenance.csv"), row.names = FALSE)
}
sam_fit <- readRDS(reference)
write.csv(data.frame(file = basename(reference), md5 = unname(tools::md5sum(reference))),
  file.path(analysis_dir, "reference_checksum.csv"), row.names = FALSE)

years <- 1983:2022
tam_obs <- lapply(sam_to_tam_obs(sam_fit), function(d) d[d$year %in% years, ])
settings <- sam_to_tam_settings(sam_fit, overrides = list(years = years))
write.csv(sam_to_tam_audit(sam_fit, settings), file.path(results, "audit.csv"), row.names = FALSE)

# Fit and compare ----
tam_fit <- do.call(fit_tam, c(list(obs = tam_obs, silent = TRUE), settings))
saveRDS(tam_fit, file.path(results, "tinyAM_fit.rds"))
sam_comparison <- sam_to_tam_comparison(sam_fit)
sam_comparison$dat$years <- years
sam_comparison$dat$is_proj <- rep(FALSE, length(years))
sam_comparison$dat$obs <- tam_obs
models <- list(SAM = sam_comparison, tinyAM = tam_fit)

diagnostics <- data.frame(model = c("SAM", "tinyAM"),
  convergence_code = c(sam_fit$opt$convergence, tam_fit$opt$convergence),
  message = c(sam_fit$opt$message, tam_fit$opt$message),
  max_abs_gradient = c(max(abs(sam_fit$sdrep$gradient.fixed)), max(abs(tam_fit$obj$gr(tam_fit$opt$par)))),
  positive_definite_hessian = c(sam_fit$sdrep$pdHess, tam_fit$sdrep$pdHess))
write.csv(diagnostics, file.path(results, "diagnostics.csv"), row.names = FALSE)
source(file.path(analysis_dir, "002_summarize.R"))
vis_tam(model_list = models, output_file = file.path(results, "SAM_tinyAM_dashboard.html"),
  open_file = FALSE, render_args = list(quiet = TRUE))
print(diagnostics, row.names = FALSE)
