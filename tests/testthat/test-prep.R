
test_that("cut_ages makes singleton labels when breaks are consecutive", {
  x <- 2:14
  f <- cut_ages(x, 2:14)
  expect_s3_class(f, "factor")
  expect_false(is.ordered(f))
  expect_equal(levels(f), as.character(2:14))
  expect_equal(as.character(unique(f)), as.character(2:14))
})

test_that("cut_ages builds hyphen ranges and closes at max(x)", {
  x <- 2:14
  f <- cut_ages(x, seq(2, 14, by = 2))
  expect_equal(levels(f), c("2-3","4-5","6-7","8-9","10-11","12-14"))
  # spot-check a few assignments
  expect_equal(as.character(f[1:2]), c("2-3","2-3"))
  expect_equal(as.character(f[12:13]), c("12-14","12-14"))
})

test_that("cut_years creates management-era blocks", {
  yrs <- 1983:2025
  f <- cut_years(yrs, c(1983, 1992, 1997, 2003, 2025))
  expect_equal(levels(f), c("1983-1991","1992-1996","1997-2002","2003-2025"))
  # first few should be first block; last few should be last block
  expect_true(all(f[1:3] == "1983-1991"))
  expect_true(all(tail(f, 3) == "2003-2025"))
})

test_that("cut_int respects ordered = TRUE", {
  f <- cut_int(2:10, c(2,5,8,10), ordered = TRUE)
  expect_true(is.ordered(f))
  expect_equal(levels(f), c("2-4","5-7","8-10"))
})

test_that("cut_int input validation errors are informative", {
  expect_error(cut_int(c(1, NA, 3), 1:3), "must be non-NA")
  expect_error(cut_int(c(1, 1.5, 3), 1:3), "integer-valued")
  expect_error(cut_int(2:5, c(2, 2, 5)), "strictly increasing")
  expect_error(cut_int(2:5, 1:4), "first break must equal min")
  expect_error(cut_int(2:5, 2:6), "last break must equal max")
})


testthat::skip_if_not(exists("northern_cod_data"), "northern_cod_data not available")

test_that("make_dat infers years/ages when NULL and builds expected pieces", {
  dat <- make_dat(
    obs = northern_cod_data,
    years = NULL,
    ages  = NULL,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw",  mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )

  expect_type(dat, "list")
  expect_true(all(c("years","ages","obs","SW","MO","obs_map","log_obs",
                    "sd_obs_modmat","q_modmat","F_settings","M_settings","N_settings") %in% names(dat)))

  # dims line up
  ny <- length(dat$years); na <- length(dat$ages)
  expect_equal(dim(dat$SW), c(ny, na))
  expect_equal(dim(dat$MO), c(ny, na))

  # design matrices exist
  expect_true(is.matrix(dat$sd_obs_modmat))
  expect_true(is.matrix(dat$q_modmat))
})

test_that("make_dat builds M age_blocks and handles M assumptions / mu_form", {
  # with age_breaks => grouped age blocks for M deviations
  dat1 <- make_dat(
    obs = northern_cod_data,
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~I(0.3), age_breaks = seq(2, 14, 2))
  )
  expect_true("age_blocks" %in% names(dat1$M_settings))
  expect_s3_class(dat1$M_settings$age_blocks, "factor")

  # If mu_form + assumption => intercept dropped (warning) and M_modmat has no intercept
  expect_warning(
    dat2 <- make_dat(
      obs = northern_cod_data,
      M_settings = list(process = "off", mu_form = ~ age, assumption = ~I(0.2)) # example mu_form with intercept (by default)
    ),
    "Dropping intercept term.*assumed levels"
  )
  expect_true(is.matrix(dat2$M_modmat))
  expect_false("(Intercept)" %in% colnames(dat2$M_modmat))
})

test_that("make_dat stops if neither M assumption nor mu_form is provided", {
  expect_error(
    make_dat(
      obs = northern_cod_data,
      M_settings = list(process = "off", mu_form = NULL, assumption = NULL)
    ),
    "Please supply an assumption or mu_form for M"
  )
})

test_that("If N process off and init_N0 FALSE, init_N0 is forced with warning", {
  expect_warning(
    make_dat(
      obs = northern_cod_data,
      N_settings = list(process = "off", init_N0 = FALSE)
    ),
    "forcing init_N0 to TRUE"
  )
})


test_that("make_par builds shapes and zeros consistent with dat", {
  dat <- make_dat(
    obs = northern_cod_data,
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


