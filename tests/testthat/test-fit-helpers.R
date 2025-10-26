test_that(".merge_start_par replaces scalars, named vectors, and matrices", {
  par0 <- list(
    scalar = 1,
    vec = c(a = 1, b = 2, c = 3),
    mat = matrix(1:4, nrow = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  )

  start <- list(
    scalar = 10,
    vec = c(b = 20, c = 30),
    mat = matrix(c(9, 8, 7, 6), nrow = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  )

  merged <- tinyAM:::.merge_start_par(par0, start)

  expect_equal(merged$scalar, 10)
  expect_equal(merged$vec, c(a = 1, b = 20, c = 30))
  expect_equal(merged$mat["r1", "c1"], 9)
  expect_equal(merged$mat["r2", "c2"], 6)
})

test_that(".merge_start_par respects vector length and dimension mismatches", {
  par0 <- list(
    vec = c(x = 1, y = 2, z = 3),
    mat = matrix(1:4, nrow = 2)
  )

  start <- list(
    vec = c(u = 9, v = 8),
    mat = matrix(5:6, nrow = 1)
  )

  merged <- tinyAM:::.merge_start_par(par0, start)

  expect_equal(merged$vec, par0$vec)
  expect_equal(merged$mat, par0$mat)
})
