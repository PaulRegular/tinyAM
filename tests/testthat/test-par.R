
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
  expect_length(par$log_q,       ncol(dat$q_modmat))
  expect_length(par$log_sd_obs,  ncol(dat$sd_obs_modmat))
  expect_equal(length(par$log_r), length(dat$years))

  # log_f rows are for historical years only (no projections in this case)
  expect_equal(dim(par$log_f), c(sum(!dat$is_proj), length(dat$ages)))
  expect_equal(sum(!dat$is_proj), length(dat$years))  # sanity for this setup

  # because N process = iid (not "off"), log_n and log_sd_n exist with right dims
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

  # everything numeric is initialized at 0
  numeric_slots <- Filter(is.numeric, par)
  expect_true(all(unlist(numeric_slots) == 0))
})

test_that("make_par shapes adapt correctly when projections are enabled", {
  dat <- make_dat(
    obs = cod_obs,
    years = NULL, ages = NULL,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~ I(0.3), age_breaks = seq(2, 14, 2)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block),
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 1)
  )
  par <- make_par(dat)

  # log_f is only for non-projection years
  expect_equal(nrow(par$log_f), sum(!dat$is_proj))
  expect_equal(ncol(par$log_f), length(dat$ages))

  # log_r has one entry per modeled year (including projections)
  expect_equal(length(par$log_r), length(dat$years))

  # log_n spans all years x (ages-1)
  expect_equal(dim(par$log_n), c(length(dat$years), length(dat$ages) - 1))

  # If M deviations are on, log_m has (years-1) x n_age_blocks
  expect_equal(dim(par$log_m), c(length(dat$years) - 1L, nlevels(dat$M_settings$age_blocks)))
})

test_that("make_par includes/excludes mean-structure parameters appropriately", {
  # With F mu_form present -> log_mu_f created with right length
  datF <- make_dat(
    obs = cod_obs,
    F_settings = list(process = "iid", mu_form = ~ F_a_block + F_y_block),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  parF <- make_par(datF)
  expect_true("log_mu_f" %in% names(parF))
  expect_length(parF$log_mu_f, ncol(datF$F_modmat))

  # Without F mu_form -> log_mu_f absent
  datF0 <- make_dat(
    obs = cod_obs,
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  parF0 <- make_par(datF0)
  expect_false("log_mu_f" %in% names(parF0))

  # With M mu_form present -> log_mu_m created
  datM <- make_dat(
    obs = cod_obs,
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = ~ 1, assumption = NULL, age_breaks = seq(2, 14, 2)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  parM <- make_par(datM)
  expect_true("log_mu_m" %in% names(parM))
  expect_length(parM$log_mu_m, ncol(datM$M_modmat))

  # With only M assumption (no mu_form) -> log_mu_m absent
  datM0 <- make_dat(
    obs = cod_obs,
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~ I(0.3), age_breaks = seq(2, 14, 2)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  parM0 <- make_par(datM0)
  expect_false("log_mu_m" %in% names(parM0))
})

test_that("make_par adds AR1 parameters only for processes set to ar1", {
  # N=ar1, F=iid, M=off
  dat <- make_dat(
    obs = cod_obs,
    N_settings = list(process = "ar1", init_N0 = FALSE),
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  par <- make_par(dat)

  expect_true("logit_phi_n" %in% names(par))
  expect_length(par$logit_phi_n, 2)
  expect_false("logit_phi_f" %in% names(par))
  expect_false("logit_phi_m" %in% names(par))
})

test_that("make_par enforces M mu_form constancy within age_blocks (when provided)", {
  # Construct M age blocks that group multiple ages, and a mu_form that varies by age
  # => should trigger the consistency check abort in make_par
  dat_bad <- make_dat(
    obs = cod_obs,
    M_settings = list(
      process = "iid",
      mu_form = ~ factor(age),  # different fixed level by age
      assumption = NULL,
      age_breaks = c(2, 3, 5, 7, 9, 11, 12, 14)  # at least one block with >1 age
    ),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )

  expect_error(
    make_par(dat_bad),
    "M mean structure varies within M age_blocks"
  )
})
