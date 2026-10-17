test_that("count-scale changes preserve the likelihood and population physics", {
  small <- cod_obs
  large <- small
  for (name in c("catch", "index")) {
    large[[name]]$obs <- large[[name]]$obs * 1000
    large[[name]]$unit <- "fish"
  }
  for (init in c("exp", "latent", "mono", "capped")) {
    a <- prepare_tam(small, years = 2000:2004, ages = 2:6, N_settings = list(init = init))
    b <- prepare_tam(large, years = 2000:2004, ages = 2:6, N_settings = list(init = init))
    p <- make_par(a)
    q <- p
    for (name in intersect(c("log_r0", "log_r", "log_r_init", "log_n0", "log_n"), names(q))) q[[name]][] <- q[[name]] + log(1000)
    oa <- RTMB::MakeADFun(function(p) nll_fun(p, a), p, silent = TRUE)
    ob <- RTMB::MakeADFun(function(p) nll_fun(p, b), q, silent = TRUE)
    expect_equal(oa$fn(), ob$fn(), tolerance = 1e-10)
    expect_equal(oa$gr(), ob$gr(), tolerance = 1e-8)
    expect_equal(oa$report()$N * 1000, ob$report()$N, tolerance = 1e-10)
    expect_equal(oa$report()$F, ob$report()$F)
    expect_equal(oa$report()$M, ob$report()$M)
  }
})
