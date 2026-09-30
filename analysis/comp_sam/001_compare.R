library(tinyAM)

# Inputs ----
case <- "WKCOD_combined_99"
analysis_dir <- file.path("analysis", "comp_sam")
results <- file.path(analysis_dir, "results")
tam_fits_dir <- file.path(analysis_dir, "tam_fits")
sam_fits_dir <- file.path(analysis_dir, "sam_fits")
reference <- file.path(sam_fits_dir, paste0(case, ".rds"))
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
    file.path(sam_fits_dir, "sam_source.csv"), row.names = FALSE)
}
sam_fit <- readRDS(reference)

years <- 1983:2022
tam_obs <- lapply(sam_to_tam_obs(sam_fit), function(d) d[d$year %in% years, ])
settings <- sam_to_tam_settings(sam_fit, overrides = list(years = years))
write.csv(sam_to_tam_audit(sam_fit, settings), file.path(results, "audit.csv"), row.names = FALSE)

# Fit and compare ----
tam_fit <- do.call(fit_tam, c(list(obs = tam_obs, silent = TRUE), settings))
saveRDS(tam_fit, file.path(tam_fits_dir, "tinyAM_fit.rds"))
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

# Annual differences ----
pop <- tidy_tam(model_list = models)$pop
metrics <- pop[c("ssb", "recruitment", "N", "F", "M")]
metrics$Fbar <- aggregate(est ~ model + year,
                          data = pop$F[pop$F$age %in% settings$F_settings$mean_ages, ], FUN = mean)

percent_differences <- do.call(rbind, lapply(names(metrics), function(metric) {
  d <- metrics[[metric]]
  if (!"age" %in% names(d)) d$age <- NA_integer_
  sam <- d[d$model == "SAM", c("year", "age", "est")]
  tam <- d[d$model == "tinyAM", c("year", "age", "est")]
  names(sam)[3] <- "SAM"
  names(tam)[3] <- "tinyAM"
  d <- merge(sam, tam, by = c("year", "age"))
  d$metric <- metric
  d$percent_difference <- ifelse(d$SAM == 0, NA_real_, 100 * (d$tinyAM / d$SAM - 1))
  d[c("metric", "year", "age", "SAM", "tinyAM", "percent_difference")]
}))
write.csv(percent_differences, file.path(results, "percent_differences.csv"), row.names = FALSE)

# Summary by metric and age ----
mean_difference <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
groups <- split(percent_differences, paste(percent_differences$metric, percent_differences$age))
summary <- do.call(rbind, lapply(groups, function(d) {
  terminal <- d[which.max(d$year), ]
  data.frame(metric = d$metric[1], age = d$age[1],
             n_years = sum(!is.na(d$percent_difference)),
             mean_percent_difference = mean_difference(d$percent_difference),
             mean_absolute_percent_difference = mean_difference(abs(d$percent_difference)),
             terminal_year = terminal$year, terminal_percent_difference = terminal$percent_difference)
}))
write.csv(summary, file.path(results, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)


vis_tam(model_list = models, output_file = file.path(results, "SAM_tinyAM_dashboard.html"),
  open_file = FALSE, render_args = list(quiet = TRUE))
print(diagnostics, row.names = FALSE)
