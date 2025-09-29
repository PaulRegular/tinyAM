
test_that("make_par builds shapes and zeros consistent with dat", {
  dat <- make_dat(
    obs = cod_obs,
    years = 1983:2024,
    ages  = 2:14,
    N_settings = list(process = "iid", init_N0 = TRUE),
    F_settings = list(process = "ar1", mu_form = ~ F_a_block + F_y_block),
    M_settings = list(process = "iid", mu_form = ~ 1, assumption = NULL, age_breaks = seq(2, 14, 2)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )

  par <- make_par(dat)

  # basic presence
  expect_true(all(c("log_sd_r","log_sd_f","log_q","log_sd_obs","log_r","log_f") %in% names(par)))

  # shapes match the design matrices / grids
  expect_length(par$log_q,     ncol(dat$q_modmat))
  expect_length(par$log_sd_obs, ncol(dat$sd_obs_modmat))
  expect_equal(length(par$log_r), length(dat$years))

  expect_equal(dim(par$log_f), c(length(dat$years), length(dat$ages)))

  # because N process = iid (not "off"), log_n and log_sd_n exist
  expect_true("log_sd_n" %in% names(par))
  expect_equal(dim(par$log_n), c(length(dat$years), length(dat$ages) - 1))

  # AR1 for F => logit_phi_f present, length 2
  expect_true("logit_phi_f" %in% names(par))
  expect_length(par$logit_phi_f, 2)

  # M process iid => log_sd_m present; mu_form => log_mu_m present with right length
  expect_true("log_sd_m" %in% names(par))
  expect_true("log_mu_m" %in% names(par))
  expect_length(par$log_mu_m, ncol(dat$M_modmat))

  # init_N0 TRUE => log_r0 exists
  expect_true("log_r0" %in% names(par))

  # missing placeholder length equals NA count in log_obs
  expect_length(par$missing, sum(is.na(dat$log_obs)))
})


