
## cut_* ----

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



## make_dat ----

test_that("make_dat infers years/ages when NULL and builds expected pieces", {
  dat <- make_dat(
    obs = cod_obs,
    years = NULL,
    ages  = NULL,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw",  mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block,
                        fill_missing = TRUE)
  )

  expect_type(dat, "list")
  expect_true(all(c("years","ages","obs","W","P","obs_map","log_obs",
                    "sd_obs_modmat","q_modmat","F_settings","M_settings","N_settings") %in% names(dat)))

  # dims line up
  ny <- length(dat$years); na <- length(dat$ages)
  expect_equal(dim(dat$W), c(ny, na))
  expect_equal(dim(dat$P), c(ny, na))

  # design matrices exist and have columns
  expect_true(is.matrix(dat$sd_obs_modmat))
  expect_true(is.matrix(dat$q_modmat))
  expect_gt(ncol(dat$sd_obs_modmat), 0)
  expect_gt(ncol(dat$q_modmat), 0)
})

test_that("make_dat builds M age_blocks and handles M assumptions / mu_form", {
  # with age_breaks => grouped age blocks for M deviations
  dat1 <- make_dat(
    obs = cod_obs,
    ages = 2:14,
    M_settings = list(process = "iid", mu_form = NULL, assumption = ~ I(0.3), age_breaks = seq(2, 14, 2))
  )
  expect_true("age_blocks" %in% names(dat1$M_settings))
  expect_s3_class(dat1$M_settings$age_blocks, "factor")

  # If mu_form + assumption => intercept dropped (warning) and M_modmat has no intercept
  expect_warning(
    dat2 <- make_dat(
      obs = cod_obs,
      M_settings = list(process = "off", mu_form = ~ age, assumption = ~ I(0.2)) # default mu_form carries intercept
    ),
    "Dropping intercept term.*assumed levels"
  )
  expect_true(is.matrix(dat2$M_modmat))
  expect_false("(Intercept)" %in% colnames(dat2$M_modmat))
})

test_that("make_dat stops if neither M assumption nor mu_form is provided", {
  expect_error(
    make_dat(
      obs = cod_obs,
      M_settings = list(process = "off", mu_form = NULL, assumption = NULL)
    ),
    "Please supply an assumption or mu_form for M"
  )
})

test_that("make_dat forces init_N0 to TRUE, with warning, when N process off and init_N0 FALSE", {
  expect_warning(
    dat <- make_dat(
      obs = cod_obs,
      N_settings = list(process = "off", init_N0 = FALSE)
    ),
    "forcing init_N0 to TRUE"
  )
  expect_true(dat$N_settings$init_N0)
})

test_that("make_dat appends projection years and shapes obs correctly with proj_settings (F_mult API)", {
  # n_proj = 2, n_mean = 3, status-quo F multiplier = 1
  dat <- make_dat(
    obs = cod_obs,
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 1)
  )

  # projected years should follow the terminal observed year
  maxy <- max(unique(cod_obs$catch$year))
  expect_true(all((maxy + 1):(maxy + 2) %in% dat$years))
  expect_equal(dat$proj_years, (maxy + 1):(maxy + 2))

  # is_proj flags are computed from years (not just concatenated)
  expect_true(all(dat$is_proj == (dat$years %in% dat$proj_years)))

  # each core table has is_proj and it matches year %in% proj_years
  for (nm in c("catch", "index", "weight", "maturity")) {
    d <- dat$obs[[nm]]
    expect_true("is_proj" %in% names(d))
    expect_true(all(d$is_proj == (d$year %in% dat$proj_years)))
  }

  # catch/index obs are NA in projection years
  expect_true(all(is.na(subset(dat$obs$catch,  is_proj)$obs)))
  expect_true(all(is.na(subset(dat$obs$index,  is_proj)$obs)))

  # weight/maturity obs in proj years equal the mean of last n_mean years by age
  n_mean <- 3L
  w <- dat$obs$weight
  m <- dat$obs$maturity
  ref_years <- (maxy - n_mean + 1):maxy
  for (a in sort(unique(w$age))) {
    mw <- mean(subset(w, !is_proj & year %in% ref_years & age == a)$obs, na.rm = TRUE)
    mm <- mean(subset(m, !is_proj & year %in% ref_years & age == a)$obs, na.rm = TRUE)
    expect_equal(unique(subset(w, is_proj & age == a)$obs), mw, tolerance = 1e-12)
    expect_equal(unique(subset(m, is_proj & age == a)$obs), mm, tolerance = 1e-12)
  }
})

test_that("make_dat recycles, validates, names proj_settings$F_mult", {
  maxy <- max(unique(cod_obs$catch$year))

  # length-1 F_mult recycled to n_proj and named by proj years
  dat1 <- make_dat(
    obs = cod_obs,
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 0.9)
  )
  expect_equal(length(dat1$proj_settings$F_mult), length(dat1$proj_years))
  expect_identical(names(dat1$proj_settings$F_mult), as.character(dat1$proj_years))
  expect_true(all(abs(dat1$proj_settings$F_mult - 0.9) < 1e-12))

  # exact-length F_mult accepted (and named)
  dat2 <- make_dat(
    obs = cod_obs,
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = c(1.1, 0.8))
  )
  expect_identical(as.numeric(dat2$proj_settings$F_mult), c(1.1, 0.8))
  expect_identical(names(dat2$proj_settings$F_mult), as.character(dat2$proj_years))

  # length mismatch (>1 and != n_proj) -> error
  expect_error(
    make_dat(
      obs = cod_obs,
      proj_settings = list(n_proj = 2, n_mean = 3, F_mult = c(0.8, 0.9, 1.0))
    ),
    "length\\(proj_settings\\$F_mult\\) must equal"
  )

  # zero F_mult triggers warning and replacement with 1e-12
  expect_warning(
    dat3 <- make_dat(
      obs = cod_obs,
      proj_settings = list(n_proj = 2, n_mean = 3, F_mult = c(0, 0.5))
    ),
    "Zero F not supported; replacing 0 with 1e-12"
  )
  expect_true(any(abs(dat3$proj_settings$F_mult - 1e-12) < 1e-18))
  expect_true(any(abs(dat3$proj_settings$F_mult - 0.5)   < 1e-18))
})

test_that("make_dat handles mean_ages correctly", {
  # Defaults to all ages when NULL
  dat1 <- make_dat(
    obs = cod_obs,
    ages = 2:5,
    F_settings = list(process = "iid", mean_ages = NULL),
    M_settings = list(process = "iid", mean_ages = NULL, assumption = ~I(0.2))
  )
  expect_equal(dat1$F_settings$mean_ages, 2:5)
  expect_equal(dat1$M_settings$mean_ages, 2:5)

  # Subset of ages is allowed
  dat2 <- make_dat(
    obs = cod_obs,
    ages = 2:5,
    F_settings = list(process = "iid", mean_ages = c(2, 4)),
    M_settings = list(process = "iid", mean_ages = c(3, 5), assumption = ~I(0.2))
  )
  expect_equal(dat2$F_settings$mean_ages, c(2, 4))
  expect_equal(dat2$M_settings$mean_ages, c(3, 5))

  # Invalid ages -> error
  expect_error(
    make_dat(
      obs = cod_obs,
      ages = 2:5,
      F_settings = list(process = "iid", mean_ages = c(6, 7))
    ),
    "F_settings\\$mean_ages"
  )
})

