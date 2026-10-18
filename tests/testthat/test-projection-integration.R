test_that("catch projections align stock-recruit parents and spawning-time biology", {
  data <- cod_obs
  data$maturity$spawn_time <- data$maturity$age / 100
  d <- prepare_tam(data, years = 2000:2006, ages = 2:6,
    N_settings = list(process = "off"),
    rec_settings = list(process = "iid", mu_form = ~ bh(ssb, lag = 2)),
    proj_settings = list(n_proj = 3, catch = c(100, 80, 60)))
  p <- make_par(d)
  p$log_r0 <- log(1e5)
  p$log_r_init[] <- log(1e5)
  p$log_f[] <- log(.15)
  p$log_sd_r <- p$log_sd_f <- log(.05)
  p <- .initialize_rec_curve(p, d, initialize_states = TRUE)
  set.seed(918)
  simulated <- nll_fun(p, d, simulate = TRUE)
  p$log_r[] <- simulated$log_r
  p$log_f[] <- simulated$log_f
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  expect_equal(unname(r$total_yield_pred[d$is_proj]), c(100, 80, 60), tolerance = 1e-7)
  parents <- d$rec$eligible - 2
  expect_true(any(d$is_proj[parents]))
  expect_equal(unname(exp(r$rec_log_parent)), unname(r$ssb[parents]), tolerance = 1e-8)
  expected_mean <- exp(p$log_sr_alpha) * r$ssb[parents] / (1 + exp(p$log_sr_beta) * r$ssb[parents])
  expect_equal(unname(exp(r$rec_log_mean)), unname(expected_mean), tolerance = 1e-8)
  for (y in which(d$is_proj)) {
    last <- length(d$ages)
    plus_ssb <- sum(r$N_plus[y, ] * d$W_plus_input[y, ] * d$P_plus_input[y, ] *
      exp(-d$spawn_time_plus[y, ] * r$Z[y, last]))
    expect_equal(unname(r$ssb_mat[y, last]), unname(plus_ssb), tolerance = 1e-8)
  }
  i <- d$obs_map$type == "index" & d$obs_map$is_proj
  cells <- cbind(as.character(d$obs_map$year[i]), as.character(d$obs_map$age[i]))
  q <- r$log_q_obs[d$obs_map$is_proj[d$obs_map$type == "index"]]
  expect_equal(as.numeric(r$log_pred[i]), as.numeric(q + log(r$N[cells]) - d$obs_map$samp_time[i] * r$Z[cells]))
})
