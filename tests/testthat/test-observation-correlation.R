test_that("cor_struct selects fleets and surveys without adding unused parameters", {
  dat <- prepare_tam(cod_obs, years = 2015:2024, ages = 2:7)
  independent <- prepare_tam(cod_obs, years = 2015:2024, ages = 2:7,
    catch_settings = list(cor_struct = "independent"),
    index_settings = list(cor_struct = "independent"))
  expect_null(dat$obs_correlation)
  expect_null(make_par(dat)$logit_rho_obs)
  expect_identical(nll_fun(make_par(dat), dat), nll_fun(make_par(independent), independent))
  expect_error(prepare_tam(cod_obs, catch_settings = list(cor_struct = "exchangeable")),
               "independent.*ar1")
  expect_error(prepare_tam(cod_obs, index_settings = list(cor_struct = c(RV = "ar1"))),
               "unknown series.*RV")
  expect_error(.obs_cor_structure(c("ar1", "ar1"), c("A", "B"), "cor_struct"),
               "scalar or name")
  expect_error(.obs_cor_structure(c(A = "ar1", A = "ar1"), "A", "cor_struct"),
               "unique")

  ar <- prepare_tam(cod_obs, years = 2015:2024, ages = 2:7,
    catch_settings = list(cor_struct = "ar1"),
    index_settings = list(cor_struct = c("Fall 2J3KL" = "ar1")))
  groups <- c("catch:Fishery", "index:Fall 2J3KL")
  expect_identical(ar$obs_correlation$groups, groups)
  par <- make_par(ar)
  expect_identical(par$logit_rho_obs, stats::setNames(c(0, 0), groups))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, ar), par, silent = TRUE)
  expect_true(all(is.finite(obj$gr())))
  expect_true(all(is.finite(obj$he())))
})

test_that("age AR1 matches the observed Gaussian submatrix with unequal SDs and gaps", {
  d <- data.frame(
    type = c(rep("catch", 5), rep("index", 6)),
    fleet = c("A", "A", "A", "B", "B", rep(NA, 6)),
    survey = c(rep(NA, 5), rep("A", 4), "B", "B"),
    year = c(1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2),
    age = c(6, 4, 2, 1, 2, 2, 1, 3, 6, 1, 2)
  )
  dat <- list(obs_map = d, has_fleet = TRUE, fleets = c("A", "B"),
    is_observed = seq_len(nrow(d)) != 2L,
    catch_settings = list(cor_struct = c(A = "ar1")),
    index_settings = list(cor_struct = c(A = "ar1")))
  dat$obs_correlation <- .prepare_obs_correlation(dat, c("A", "B"))
  expect_identical(dat$obs_correlation$groups, c("catch:A", "index:A"))
  expect_equal(dat$obs_correlation$gap, c(4, 1, 3))
  k <- which(dat$is_observed)
  mu <- seq_along(k) / 4
  sd <- seq_along(k) / 10 + .1
  y <- mu + c(-.7, .3, .1, -.2, .5, -.1, .3, -.7, .1, .2)
  reference <- function(rho) {
    covariance <- diag(sd^2)
    for (b in seq_along(rho)) {
      rows <- if (b == 1L) c(1L, 3L) else c(6L, 8L, 9L)
      i <- match(rows, k)
      covariance[i, i] <- outer(sd[i], sd[i]) * rho[b]^abs(outer(d$age[rows], d$age[rows], "-"))
    }
    L <- t(chol(covariance))
    sum(log(diag(L))) + .5 * sum(forwardsolve(L, y - mu)^2) + length(y) / 2 * log(2 * pi)
  }
  par <- list(logit_rho_obs = stats::qlogis(c(.6, .8)))
  actual <- .obs_cor_nll(dat, par, y, mu, sd)
  expect_equal(actual, reference(c(.6, .8)), tolerance = 1e-12)
  expect_equal(.obs_cor_nll(dat, list(logit_rho_obs = c(-Inf, -Inf)), y, mu, sd),
               reference(c(0, 0)), tolerance = 1e-12)
  tape <- RTMB::MakeADFun(function(p) .obs_cor_nll(dat, p, y, mu, sd), par, silent = TRUE)
  expect_equal(tape$fn(), actual, tolerance = 1e-12)
  step <- 1e-5
  shift <- diag(step, 2)
  gradient <- apply(shift, 2, function(h) {
    (reference(stats::plogis(tape$par + h)) - reference(stats::plogis(tape$par - h))) / (2 * step)
  })
  expect_equal(drop(tape$gr()), gradient, tolerance = 1e-7)

  dat$catch_settings$cor_struct <- "independent"
  dat$index_settings$cor_struct <- c(B = "ar1")
  expect_warning(correlation <- .prepare_obs_correlation(dat, c("A", "B")),
                 "no within-year age pairs")
  expect_null(correlation)
})
