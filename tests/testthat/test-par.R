
YEARS <- 1983:2024
AGES  <- 2:14

test_dat <- function(obs = cod_obs,
                     years = YEARS,
                     ages  = AGES,
                     N_settings = list(process = "iid", init_N0 = FALSE),
                     F_settings = list(process = "approx_rw",  mu_form = NULL),
                     M_settings = list(process = "off", mu_form = NULL, assumption = ~I(0.3)),
                     obs_settings = list(q_form = ~ q_block, sd_catch_form = ~1,
                                         sd_index_form = ~1, fill_missing = TRUE),
                     proj_settings = NULL) {
  args <- mget(ls())
  do.call(make_dat, args)
}


test_that("make_par builds shapes and zeros consistent with dat", {
  dat <- test_dat(
    N_settings = list(process = "iid", init_N0 = TRUE),
    F_settings = list(process = "ar1", mu_form = ~ F_a_block + F_y_block),
    M_settings = list(process = "iid", mu_form = ~ 1, assumption = NULL, age_breaks = seq(2, 14, 2))
  )
  par <- make_par(dat)

  # basic presence
  expect_true(all(c("log_sd_r","log_sd_f","log_q","log_sd_catch", "log_sd_index","log_r","log_f") %in% names(par)))

  # shapes match the design matrices / grids
  expect_length(par$log_q,      ncol(dat$q_modmat))
  expect_length(par$log_sd_catch, ncol(dat$sd_catch_modmat))
  expect_length(par$log_sd_index, ncol(dat$sd_index_modmat))
  expect_equal(length(par$log_r), length(dat$years))

  # log_f rows are for historical years only (no projections here)
  expect_equal(dim(par$log_f), c(sum(!dat$is_proj), length(dat$ages)))
  expect_equal(sum(!dat$is_proj), length(dat$years))  # sanity

  # N process = iid -> log_n + log_sd_n
  expect_true("log_sd_n" %in% names(par))
  expect_equal(dim(par$log_n), c(length(dat$years), length(dat$ages) - 1))

  # AR1 for F => logit_phi_f present (length 2)
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

  # everything numeric initialized at 0
  numeric_slots <- Filter(is.numeric, par)
  expect_true(all(unlist(numeric_slots) == 0))
})

test_that("make_par shapes adapt when projections are enabled", {
  dat <- test_dat(
    years = NULL, ages = 2:14,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~ I(0.3), age_breaks = seq(2, 14, 2)),
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 1)
  )
  par <- make_par(dat)

  expect_equal(nrow(par$log_f), sum(!dat$is_proj))
  expect_equal(ncol(par$log_f), length(dat$ages))

  expect_equal(length(par$log_r), length(dat$years))
  expect_equal(dim(par$log_n), c(length(dat$years), length(dat$ages) - 1))

})

test_that("make_par includes/excludes mean-structure parameters appropriately", {
  # F mu_form present -> log_mu_f created
  datF  <- test_dat(
    F_settings = list(process = "iid", mu_form = ~ F_a_block + F_y_block),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3))
  )
  parF <- make_par(datF)
  expect_true("log_mu_f" %in% names(parF))
  expect_length(parF$log_mu_f, ncol(datF$F_modmat))

  # F mu_form absent -> log_mu_f absent
  datF0 <- test_dat(
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3))
  )
  parF0 <- make_par(datF0)
  expect_false("log_mu_f" %in% names(parF0))

  # M mu_form present -> log_mu_m created
  datM  <- test_dat(
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = ~ 1, assumption = NULL, age_breaks = seq(2, 14, 2))
  )
  parM <- make_par(datM)
  expect_true("log_mu_m" %in% names(parM))
  expect_length(parM$log_mu_m, ncol(datM$M_modmat))

  # Only M assumption -> log_mu_m absent
  datM0 <- test_dat(
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~ I(0.3), age_breaks = seq(2, 14, 2))
  )
  parM0 <- make_par(datM0)
  expect_false("log_mu_m" %in% names(parM0))
})

test_that("make_par adds AR1 parameters only for processes set to ar1", {
  dat <- test_dat(
    N_settings = list(process = "ar1", init_N0 = FALSE),
    F_settings = list(process = "iid", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3))
  )
  par <- make_par(dat)

  expect_true("logit_phi_n" %in% names(par))
  expect_length(par$logit_phi_n, 2)
  expect_false("logit_phi_f" %in% names(par))
  expect_false("logit_phi_m" %in% names(par))
})

test_that("make_par enforces M mu_form constancy within age_blocks (when provided)", {
  dat_bad <- test_dat(
    M_settings = list(
      process = "iid",
      mu_form = ~ factor(age),  # varies by age
      assumption = NULL,
      age_breaks = c(2, 3, 5, 7, 9, 11, 12, 14)  # some blocks with >1 age
    )
  )
  expect_error(
    make_par(dat_bad),
    "M mean structure varies within M age_blocks"
  )
})

test_that("make_par starts M deviations at year 2 by default", {
  dat <- test_dat(M_settings = list(process = "iid", mu_form = ~1))
  par <- make_par(dat)
  expect_equal(rownames(par$log_m)[1], as.character(dat$years[2]))

  dat <- test_dat(M_settings = list(process = "iid", mu_form = ~1, first_dev_year = 1983))
  par <- make_par(dat)
  expect_equal(rownames(par$log_m)[1], as.character(dat$years[1]))
})
