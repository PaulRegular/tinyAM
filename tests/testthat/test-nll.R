set.seed(1)

ny <- 40; na <- 20; sd <- 0.5

test_that("rprocess_2d returns matrix with correct dimensions", {
  X <- rprocess_2d(ny, na, sd = 0.3, phi = c(0.5, 0.7))
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(ny, na))
  expect_true(all(is.finite(X)))
})

test_that("rprocess_2d gives approx iid entries with Var ~ sd^2 when phi = (0,0)", {
  X  <- rprocess_2d(ny, na, sd = sd, phi = c(0, 0))

  # Empirical variance of all entries
  v_all <- var(as.numeric(X))
  expect_equal(v_all, sd^2, tolerance = 0.1)

  # Adjacent lag-1 correlations (rows and cols) should be near zero
  r_row <- cor(as.numeric(X[, -na]), as.numeric(X[, -1]))
  r_col <- cor(as.numeric(X[-ny, ]), as.numeric(X[-1, ]))
  expect_equal(r_row, 0, tolerance = 0.1)
  expect_equal(r_col, 0, tolerance = 0.1)
})

test_that("dprocess_2d AR1 density reduces to IID when phi = (0,0)", {
  X  <- matrix(rnorm(ny * na, 0, sd), ny, na)

  lp_iid <- dprocess_2d(X, sd = sd, phi = c(0, 0))
  # Compare to the literal IID log-density
  lp_ref <- sum(dnorm(X, mean = 0, sd = sd, log = TRUE))

  expect_equal(lp_iid, lp_ref, tolerance = 1e-10)
})

test_that("rprocess_2d supplied larger |phi| produces stronger local correlation (empirical check)", {
  X_weak <- rprocess_2d(ny, na, sd = sd, phi = c(0.3, 0.3))
  X_strg <- rprocess_2d(ny, na, sd = sd, phi = c(0.9, 0.9))

  # Lag-1 correlations along rows/cols
  r_row_weak <- cor(as.numeric(X_weak[, -na]), as.numeric(X_weak[, -1]))
  r_col_weak <- cor(as.numeric(X_weak[-ny, ]), as.numeric(X_weak[-1, ]))

  r_row_strg <- cor(as.numeric(X_strg[, -na]), as.numeric(X_strg[, -1]))
  r_col_strg <- cor(as.numeric(X_strg[-ny, ]), as.numeric(X_strg[-1, ]))

  expect_gt(r_row_strg, r_row_weak)
  expect_gt(r_col_strg, r_col_weak)
})

test_that("rprocess_2d self-consistency: sample from AR1 has higher density under matching phi than IID", {
  phi <- c(0.8, 0.7)

  X <- rprocess_2d(ny, na, sd = sd, phi = phi)
  lp_match <- dprocess_2d(X, sd = sd, phi = phi)
  lp_iid   <- dprocess_2d(X, sd = sd, phi = c(0, 0))

  expect_gt(lp_match, lp_iid)
})

test_that("rprocess_2d approximate RW: phi=(0.99,0.99) yields very strong adjacent correlation", {
  X <- rprocess_2d(ny, na, sd = sd, phi = c(0.99, 0.99))

  r_row <- cor(as.numeric(X[, -na]), as.numeric(X[, -1]))
  r_col <- cor(as.numeric(X[-ny, ]), as.numeric(X[-1, ]))

  expect_gt(r_row, 0.9)
  expect_gt(r_col, 0.9)

  # Density is finite and computable
  lp <- dprocess_2d(X, sd = sd, phi = c(0.99, 0.99))
  expect_true(is.finite(lp))
})

test_that("rprocess_2d edge cases: ny=1 or na=1 behave like 1D AR(1)", {
  # Single row (varying age)
  ny <- 1; na <- 80; sd <- 1.0; phi <- c(0.0, 0.8) # phi_age=0, phi_year=0.8 (ignored since ny=1)
  X1 <- rprocess_2d(ny, na, sd = sd, phi = phi)
  expect_equal(dim(X1), c(1, na))
  # Lag-1 corr across ages should be ~phi_age = 0
  r_age <- cor(X1[, -na], X1[, -1])
  expect_equal(as.numeric(r_age), 0, tolerance = 0.15)

  # Single column (varying year)
  ny <- 80; na <- 1; sd <- 1.0; phi <- c(0.8, 0.4) # now only year AR(1) matters
  X2 <- rprocess_2d(ny, na, sd = sd, phi = phi)
  expect_equal(dim(X2), c(ny, 1))
  # Lag-1 corr down years should be ~phi_year = 0.4
  r_year <- cor(X2[-ny, 1], X2[-1, 1])
  expect_equal(as.numeric(r_year), 0.4, tolerance = 0.15)
})



test_that("solve_F_mult returns ~1 when target equals catch at k=1", {
  F <- rep(0.3, 6); M <- rep(0.2, 6); N <- 1e3 * exp(-0.4 * (0:5))
  C_of_k <- function(k) sum(N * (1 - exp(-(k*F + M))) * (k*F / (k*F + M)))
  C1 <- C_of_k(1)
  khat <- solve_F_mult(C1, F, M, N, k_init = 1)
  expect_true(is.finite(khat))
  expect_equal(khat, 1, tolerance = 1e-6)
})

test_that("solve_F_mult scales sensibly with target", {
  F <- seq(0.2, 0.35, length.out = 5); M <- rep(0.25, 5); N <- rep(500, 5)
  C_of_k <- function(k) sum(N * (1 - exp(-(k*F + M))) * (k*F / (k*F + M)))
  C1 <- C_of_k(1)
  C2 <- 1.5 * C1
  k1 <- solve_F_mult(C1, F, M, N, k_init = 1)
  k2 <- solve_F_mult(C2, F, M, N, k_init = 1)
  expect_gt(k2, k1)
})



test_that("nll_fun returns finite JNLL and simulates when requested", {
  dat <- make_dat(
    cod_obs,
    years = 1983:2024,
    ages  = 2:14,
    N_settings = list(process = "iid", init_N0 = TRUE),
    F_settings = list(process = "approx_rw",  mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
    obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
  )
  par <- make_par(dat)

  # Likelihood branch
  jnll <- nll_fun(par, dat, simulate = FALSE)
  expect_true(is.numeric(jnll))
  expect_length(jnll, 1L)
  expect_true(is.finite(jnll))

  # Simulation branch (no RTMB object needed)
  sims <- nll_fun(par, dat, simulate = TRUE)
  expect_type(sims, "list")
  expect_true(all(c("log_f", "log_r", "log_obs") %in% names(sims)))
  if (dat$N_settings$process != "off") {
    expect_true("log_n" %in% names(sims))
  }
  # Simulated obs restored to NA where missing
  is_miss <- is.na(dat$log_obs)
  expect_true(all(is.na(sims$log_obs[is_miss])))
  expect_true(all(is.finite(sims$log_obs[!is_miss])))

  # Optional: build a quick RTMB obj to ensure REPORT/ADREPORT paths are valid
  make_nll_fun <- function(f, d) function(p) f(p, d)
  obj <- RTMB::MakeADFun(make_nll_fun(nll_fun, dat), par, random = c("log_f","log_r","missing"), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 50, eval.max = 50))
  expect_true(is.finite(opt$objective))
  rep <- obj$report()
  expect_true(all(c("N", "F", "M", "Z", "ssb", "log_pred", "log_q_obs") %in% names(rep)))
})

