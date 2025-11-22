
## Fit two models ----
N_dev <- fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1),
  silent = TRUE
)
M_dev <- update(
  N_dev,
  N_settings = list(process = "off", init_N0 = TRUE),
  M_settings = list(process = "ar1", mu_supplied = ~ I(0.3),
                    age_breaks = c(3, 14))
)

## Combine model objects into named list ----
fits <- list("N_dev" = N_dev, "M_dev" = M_dev)

