
## Largely variable names captured by RTMB::get_all()
utils::globalVariables(c(
  # data-frame column names
  "year", "age", "is_proj",
  # captured-from-dat/pars in nll_fun / make_par
  "years", "ages", "proj_settings", "N_settings", "F_settings", "M_settings",
  "catch_settings", "index_settings", "log_sd_r", "log_sd_f", "log_sd_n", "log_sd_m", "log_sd_catch", "log_sd_index",
  "logit_phi_n", "logit_phi_f", "logit_phi_m",
  "log_r0", "log_mu_f", "mu_m", "log_mu_supplied_m", "log_q",
  "F_modmat", "M_modmat", "q_modmat", "sd_catch_modmat", "sd_index_modmat", "obs_map",
  "log_sd_catch_supplied", "log_sd_index_supplied", "fill_missing_map", "any_fill_missing",
  "is_missing", "is_observed", "W", "P"
))
