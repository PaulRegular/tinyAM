# Mean fishing patterns ----
monotone <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off"),
  F_settings = list(process = "iid", mu_form = ~ mono(age)),
  M_settings = list(mu_supplied = ~ I(.3)))

logistic_iid <- update(monotone, F_settings = list(mu_form = ~ logistic(age)))
logistic_ar1 <- update(logistic_iid, F_settings = list(process = "ar1"))

# Compare diagnostics and mean curves ----
check_tam(monotone)
check_tam(logistic_iid)
check_tam(logistic_ar1)
head(tidy_pop(logistic_ar1)$mu_F)
if (interactive()) vis_tam(monotone, logistic_iid, logistic_ar1)
