
## Largely variable names captured by RTMB::get_all()
utils::globalVariables(c(
  # data-frame column names
  "year", "age", "is_proj",
  # captured-from-dat/pars in nll_fun / make_par
  "years", "ages", "proj_settings", "N_settings", "F_settings", "M_settings",
  "log_sd_r", "log_sd_f", "log_sd_n", "log_sd_m", "log_sd_obs",
  "logit_phi_n", "logit_phi_f", "logit_phi_m",
  "log_r0", "log_mu_f", "log_mu_m", "log_mu_assumed_m", "log_q",
  "F_modmat", "M_modmat", "q_modmat", "sd_obs_modmat", "obs_map",
  "is_na_obs"
))
