
## Fit two models ----
N_dev <- fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "iid", init = "exp"),
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  silent = TRUE
)
N_dev <- update(N_dev,
  proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1),
  start_par = as.list(N_dev$sdrep, "Estimate")
)
M_dev <- update(
  N_dev,
  N_settings = list(process = "off", init = "exp"),
  M_settings = list(process = "ar1", mu_supplied = ~ I(0.3),
                    age_breaks = c(3, 14))
)

## Combine model objects into named list ----
fits <- list("N_dev" = N_dev, "M_dev" = M_dev)

