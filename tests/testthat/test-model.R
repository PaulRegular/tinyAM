
set.seed(1)

test_that("rprocess_2d returns matrices of requested size", {
  X1 <- rprocess_2d(7, 5, sd = 0.3, type = "iid")
  X2 <- rprocess_2d(7, 5, sd = 0.3, type = "rw")
  X3 <- rprocess_2d(7, 5, sd = 0.3, type = "ar1", phi = c(0.4, 0.7))
  expect_equal(dim(X1), c(7, 5))
  expect_equal(dim(X2), c(7, 5))
  expect_equal(dim(X3), c(7, 5))
  expect_true(is.numeric(X1) && is.numeric(X2) && is.numeric(X3))
})

test_that("dprocess_2d IID matches sum of univariate normals", {
  ny <- 6; na <- 4; sd <- 0.5
  X  <- matrix(rnorm(ny * na, 0, sd), ny, na)
  lp_ref <- sum(dnorm(X, mean = 0, sd = sd, log = TRUE))
  lp_fun <- dprocess_2d(X, sd = sd, type = "iid")
  expect_equal(lp_fun, lp_ref, tolerance = 1e-10)
})

test_that("dprocess_2d RW is translation-invariant (intrinsic field)", {
  ny <- 8; na <- 6; sd <- 0.3
  X  <- rprocess_2d(ny, na, sd = sd, type = "rw")
  c0 <- dprocess_2d(X, sd = sd, type = "rw")
  c1 <- dprocess_2d(X + 5, sd = sd, type = "rw")      # add constant
  c2 <- dprocess_2d(X + outer(rep(1, ny), 1:na), sd = sd, type = "rw") # add age slope: diff-in-age unchanged
  expect_equal(c0, c1, tolerance = 1e-10)
  expect_equal(c0, c2, tolerance = 1e-10)
})

test_that("dprocess_2d AR1 with phi=0 reduces to IID with same sd", {
  ny <- 6; na <- 5; sd <- 0.7
  X  <- matrix(rnorm(ny * na, 0, sd), ny, na)
  lp_iid <- dprocess_2d(X, sd = sd, type = "iid")
  lp_ar1 <- dprocess_2d(X, sd = sd, type = "ar1", phi = c(0, 0))
  expect_equal(lp_ar1, lp_iid, tolerance = 1e-10)
})

test_that("dprocess_2d AR1 with phi=1 reduces to a random walk with same sd", {
  ny <- 6; na <- 5; sd <- 0.7
  X  <- matrix(rnorm(ny * na, 0, sd), ny, na)
  lp_rw <- dprocess_2d(X, sd = sd, type = "rw")
  lp_ar1 <- dprocess_2d(X, sd = sd, type = "ar1", phi = c(0.9999, 0.9999))
  expect_equal(lp_ar1, lp_rw, tolerance = 1e-10)
})

test_that("rprocess_2d RW produces correct diff variances (approximate)", {
  # Differences along each axis should have variance ~ sd^2
  ny <- 40; na <- 30; sd <- 0.5
  X  <- rprocess_2d(ny, na, sd = sd, type = "rw")
  var_dy <- var(as.numeric(apply(X, 2, diff)))
  var_da <- var(as.numeric(t(apply(X, 1, diff))))
  expect_equal(var_dy, sd^2, tolerance = 0.1)  # loose tolerance (finite sample)
  expect_equal(var_da, sd^2, tolerance = 0.1)
})

test_that("rprocess_2d AR1 yields sensible empirical correlations", {
  ny <- 30; na <- 30; sd <- 0.4; phi <- c(0.6, 0.8)
  X  <- rprocess_2d(ny, na, sd = sd, type = "ar1", phi = phi)
  # Estimate lag-1 corr down rows (year) using a middle column:
  cy <- cor(X[-1, round(na/2)], X[-ny, round(na/2)])
  # Estimate lag-1 corr across columns (age) using a middle row:
  ca <- cor(X[round(ny/2), -1], X[round(ny/2), -na])
  expect_equal(cy, phi[2], tolerance = 0.1)
  expect_equal(ca, phi[1], tolerance = 0.1)
})
