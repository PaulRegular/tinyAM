plus_biology_obs <- function() {
  grid <- expand.grid(year = 2000:2003, age = 2:7)
  catch <- index <- weight <- maturity <- grid
  catch$obs <- 20 + grid$age
  index$obs <- 100 + grid$age
  index$survey <- "survey"
  index$samp_time <- .5
  weight$obs <- grid$age^2 / 10 + .05 * (grid$year - 2000)
  maturity$obs <- (grid$age - 2) / 6
  list(catch = catch, index = index, weight = weight, maturity = maturity)
}

plus_biology_dat <- function(process = "iid", ages = 2:4, proj = NULL, init = "exp") {
  suppressWarnings(make_dat(plus_biology_obs(), ages = ages,
    N_settings = list(process = process, init = init),
    F_settings = list(process = "iid"),
    M_settings = list(process = "off", mu_supplied = ~ I(.2)),
    index_settings = list(q_form = ~ 1, sd_form = ~ 1, fill_missing = TRUE),
    proj_settings = proj))
}

plus_biology_par <- function(dat) {
  par <- make_par(dat)
  par$log_r0 <- log(500)
  par$log_r[] <- log(seq(450, 550, length.out = length(par$log_r)))
  par$log_f[] <- log(.2) + seq(0, .3, length.out = length(par$log_f))
  if (!is.null(par$log_n)) par$log_n[] <- log(seq(100, 300, length.out = length(par$log_n)))
  if (!is.null(par$log_n0)) par$log_n0[] <- log(seq(150, 250, length.out = length(par$log_n0)))
  for (nm in grep("^log_sd", names(par), value = TRUE)) par[[nm]][] <- log(.1)
  par
}

plus_biology_report <- function(dat, par = plus_biology_par(dat)) {
  RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)$report()
}

test_that("biological inputs retain all hidden ages and project each age separately", {
  dat <- plus_biology_dat(proj = list(n_proj = 2, n_mean = 2, F_mult = c(.8, 1.2)))
  expect_equal(dat$plus_ages, 4:7)
  expect_equal(dimnames(dat$W_plus_input), list(year = as.character(2000:2005), age = as.character(4:7)))
  obs <- plus_biology_obs()
  for (nm in c("weight", "maturity")) {
    input <- dat[[if (nm == "weight") "W_plus_input" else "P_plus_input"]]
    expected <- matrix(subset(obs[[nm]], age >= 4)$obs, 4, 4)
    expect_equal(unname(input[1:4, ]), expected)
    expect_equal(unname(input[5, ]), colMeans(expected[3:4, ]))
    expect_equal(input[5, ], input[6, ])
    expect_equal(dat$obs[[nm]]$obs[dat$obs[[nm]]$age == 4], unname(input[, 1]))
    expect_true(all(dat$obs[[nm]]$age %in% dat$ages))
  }
  expect_identical(tinyAM:::.plus_fun(obs, 4)$weight, obs$weight)
  expect_identical(tinyAM:::.plus_fun(obs, 4)$maturity, obs$maturity)
  # Retained biological ages must not introduce unused M formula coefficients.
  obs$weight$bio_age <- factor(obs$weight$age)
  modeled <- make_dat(obs, ages = 2:4,
    M_settings = list(process = "off", mu_form = ~ bio_age, mu_supplied = NULL),
    index_settings = list(q_form = ~ 1, sd_form = ~ 1, fill_missing = TRUE))
  expect_equal(levels(modeled$obs$weight$bio_age), as.character(2:4))
  expect_equal(ncol(modeled$M_modmat), 3)
  expect_equal(length(make_par(modeled)$mu_m), 3)
})

test_that("hidden age composition conserves N, biomass and mature biomass", {
  for (process in c("off", "iid", "rw", "ar1")) {
    dat <- plus_biology_dat(process)
    rep <- plus_biology_report(dat)
    total <- rep$N[, "4"]
    expect_equal(rowSums(rep$N_plus), total, tolerance = 1e-12)
    survival <- exp(-rep$Z[1, "4"])
    expected_shares <- c((1 - survival) * survival^(0:2), survival^3)
    expect_equal(unname(rep$N_plus[1, ] / total[1]), expected_shares, tolerance = 1e-12)
    for (y in 2:nrow(rep$N)) {
      pred <- c(rep$N[y - 1, "3"] * exp(-rep$Z[y - 1, "3"]),
                rep$N_plus[y - 1, 1:3] * exp(-rep$Z[y - 1, "4"]))
      pred[4] <- pred[4] + rep$N_plus[y - 1, 4] * exp(-rep$Z[y - 1, "4"])
      ordinary_prediction <- sum(rep$N[y - 1, c("3", "4")] * exp(-rep$Z[y - 1, c("3", "4")]))
      expect_equal(sum(pred), ordinary_prediction, tolerance = 1e-12)
      expect_equal(unname(rep$N_plus[y, ]), unname(pred * total[y] / sum(pred)), tolerance = 1e-12)
    }
    biomass <- rowSums(rep$N_plus * dat$W_plus_input)
    mature <- rowSums(rep$N_plus * dat$W_plus_input * dat$P_plus_input)
    expect_equal(rep$W[, "4"], biomass / total, tolerance = 1e-12)
    expect_equal(rep$P[, "4"], mature / biomass, tolerance = 1e-12)
    expect_true(all(rep$P >= 0 & rep$P <= 1))
    expect_equal(rep$biomass_mat[, "4"], biomass, tolerance = 1e-12)
    expect_equal(rep$ssb_mat[, "4"], mature, tolerance = 1e-12)
    expect_equal(rep$W[, 1:2], dat$W[, 1:2])
    expect_equal(rep$P[, 1:2], dat$P[, 1:2])
    expect_false(isTRUE(all.equal(rep$P[, "4"], rowSums(rep$N_plus * dat$P_plus_input) / total)))
  }
})

test_that("hidden biology uses realized initial and simulated states through projections", {
  for (init in c("exp", "free", "random")) {
    dat <- plus_biology_dat(init = init, proj = list(n_proj = 2, n_mean = 2, F_mult = 1))
    par <- plus_biology_par(dat)
    set.seed(204)
    sim <- nll_fun(par, dat, simulate = TRUE)
    shared <- intersect(names(par), names(sim))
    par[shared] <- sim[shared]
    rep <- plus_biology_report(dat, par)
    expect_equal(rowSums(rep$N_plus), rep$N[, "4"], tolerance = 1e-12)
    expect_equal(rep$ssb_mat[, "4"], rowSums(rep$N_plus * dat$W_plus_input * dat$P_plus_input), tolerance = 1e-12)
    expect_true(all(is.finite(rep$biomass)))
    expect_equal(rownames(rep$N_plus), as.character(dat$years))
    pop <- tidy_rep(list(dat = dat, rep = rep))
    expect_equal(sum(pop$N_plus$is_proj), 2 * length(dat$plus_ages))
  }
})

test_that("no hidden ages leaves existing population and biology equations unchanged", {
  dat <- plus_biology_dat(ages = 2:7)
  expect_null(dat$plus_ages)
  expect_null(dat$W_plus_input)
  expect_null(dat$P_plus_input)
  rep <- plus_biology_report(dat)
  expect_null(rep$N_plus)
  expect_equal(rep$biomass_mat, rep$N * dat$W)
  expect_equal(rep$ssb_mat, rep$N * dat$W * dat$P)
  expect_equal(rep$biomass, rowSums(rep$N * dat$W))
  expect_equal(rep$ssb, rowSums(rep$N * dat$W * dat$P))
})

test_that("hidden terminal group works with one extra age and zero biological weight", {
  dat <- plus_biology_dat(ages = 2:6)
  rep <- plus_biology_report(dat)
  expect_equal(ncol(rep$N_plus), 2)
  expect_equal(rowSums(rep$N_plus), rep$N[, "6"], tolerance = 1e-12)
  dat$W_plus_input[1, ] <- 0
  rep <- plus_biology_report(dat)
  expect_equal(rep$W[1, "6"], 0)
  expect_equal(rep$P[1, "6"], 0)
  expect_equal(rep$ssb_mat[1, "6"], 0)
})

test_that("hidden biology adds no likelihood terms and has valid summary derivatives", {
  dat <- plus_biology_dat("off")
  par <- plus_biology_par(dat)
  plain <- dat
  plain$plus_ages <- plain$W_plus_input <- plain$P_plus_input <- NULL
  expect_identical(make_par(dat), make_par(plain))
  expect_equal(nll_fun(par, dat), nll_fun(par, plain))
  full <- plus_biology_report(dat, par)
  unchanged <- plus_biology_report(plain, par)
  for (nm in c("N", "F", "M", "Z", "log_pred", "sd_obs", "log_q_obs")) {
    expect_equal(full[[nm]], unchanged[[nm]])
  }
  # The uncertainty calculation must differentiate the biological weighting too.
  objective <- function(p) {
    parameters <- par
    parameters$log_f <- par$log_f * 0 + p$theta
    nll_fun(parameters, dat)
  }
  obj <- RTMB::MakeADFun(objective, list(theta = log(.2)), ADreport = TRUE, silent = TRUE)
  h <- 1e-5
  finite_difference <- (obj$fn(obj$par + h) - obj$fn(obj$par - h)) / (2 * h)
  expect_true(all(is.finite(obj$gr(obj$par))))
  expect_equal(c(obj$gr(obj$par)), unname(finite_difference), tolerance = 1e-7)
})
