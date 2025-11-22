dat <- make_dat(
  cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(sd_form = ~ 1, q_form = ~ q_block)
)
