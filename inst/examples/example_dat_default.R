# Load package example observations
data(cod_obs)

# Build a TAM data list with the defaults used throughout the documentation
# (explicitly listing the main settings for clarity).
dat <- make_dat(
  cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", assumption = ~ I(0.3)),
  obs_settings = list(sd_catch_form = ~ 1, sd_index_form = ~ 1, q_form = ~ q_block)
)
