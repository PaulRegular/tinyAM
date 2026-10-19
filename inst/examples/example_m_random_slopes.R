# Supply the annual capelin/cod covariate before running this example.
data <- cod_obs
if (!"log_cap_per_cod" %in% names(data$pop_covariates)) {
  cli::cli_abort("Add your supplied log_cap_per_cod column to data$pop_covariates before running this example. No capelin series is substituted or inferred.")
}

breaks <- c(2, 4, 6, 8, 14)
selected <- data$pop_covariates$age %in% 2:14
data$pop_covariates$M_age_group <- NA_character_
data$pop_covariates$M_age_group[selected] <- as.character(
  cut_ages(data$pop_covariates$age[selected], breaks))
data$pop_covariates$M_age_group <- factor(data$pop_covariates$M_age_group)
check_tam_data(data)

# A common food-response slope, followed by partially pooled group slopes.
# mean_ages controls the reported average M, not the ages receiving the effect.
fixed_slope <- fit_tam(data, years = 1983:2025, ages = 2:14,
  N_settings = list(process = "off"),
  F_settings = list(process = "ar1", mu_form = ~ 1, mean_ages = 5:14),
  M_settings = list(process = "ar1", mu_supplied = ~ I(.3),
    mu_form = ~ 0 + log_cap_per_cod, age_breaks = breaks, mean_ages = 5:14),
  index_settings = list(sd_form = ~ 1, q_form = ~ mono(q_block)))
random_slopes <- update(fixed_slope, M_settings = list(
  mu_form = ~ 0 + log_cap_per_cod + iid(M_age_group, by = log_cap_per_cod)))
check_tam(random_slopes)
dashboard <- vis_tam(model_list = list(Fixed = fixed_slope, Random = random_slopes), open_file = FALSE)
