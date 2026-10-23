test_that("annual fixed recruitment carries its terminal effect into projections", {
  for (form in list(~factor(year), ~0 + as.factor(year), ~0 + base::factor(x = year),
                   ~factor(year, levels = 1992:1983))) {
    args <- list(data = cod_obs, years = 1983:1990, ages = 2:6,
      N_settings = list(process = "off"), rec_settings = list(process = "off", mu_form = form),
      F_settings = list(process = "off", mu_form = ~factor(year) + factor(age)))
    historical <- do.call(prepare_tam, args)
    expect_message(projected <- do.call(prepare_tam, c(args, list(proj_settings = list(n_proj = 2)))),
      "terminal year-specific fixed effect")
    expect_equal(projected$rec$matrix[!projected$is_proj, ], historical$rec$matrix, ignore_attr = TRUE)
    expect_identical(colnames(projected$rec$matrix), colnames(historical$rec$matrix))
    expect_identical(projected$rec$data$year, projected$years)
    expect_equal(projected$rec$matrix[projected$is_proj, ],
      historical$rec$matrix[rep(nrow(historical$rec$matrix), 2), ], ignore_attr = TRUE)
    p <- make_par(projected)
    p$log_r0 <- log(100)
    p$rec_beta[] <- seq_along(p$rec_beta) / 10
    objective <- RTMB::MakeADFun(function(p) nll_fun(p, projected), p, silent = TRUE)
    past_parameters <- .merge_start_par(make_par(historical), p)
    past_objective <- RTMB::MakeADFun(function(p) nll_fun(p, historical), past_parameters, silent = TRUE)
    expect_identical(names(objective$par), names(past_objective$par))
    expect_equal(objective$fn(objective$par), past_objective$fn(past_objective$par), tolerance = 1e-10)
    expect_equal(objective$gr(objective$par), past_objective$gr(past_objective$par), tolerance = 1e-10)
    expect_equal(objective$he(objective$par), past_objective$he(past_objective$par), tolerance = 1e-10)
    recruitment <- objective$report()$recruitment
    expect_equal(unname(tail(recruitment, 2)), rep(unname(recruitment[8]), 2))
    expect_true(all(is.finite(objective$gr(objective$par))))
    simulated <- nll_fun(p, projected, simulate = TRUE)
    expect_equal(unname(simulated$log_r), unname(log(recruitment[-1])))
  }
})

test_that("deterministic numeric time trends keep their projection basis", {
  for (form in list(~I(year - 1983), ~poly(year, 2),
                   ~I(year - 1983) + factor(year, levels = 1983:1992,
                     labels = rep(c("early", "late"), each = 5)))) {
    args <- list(data = cod_obs, years = 1983:1990, ages = 2:6,
      rec_settings = list(process = "off", mu_form = form))
    historical <- do.call(prepare_tam, args)
    projected <- do.call(prepare_tam, c(args, list(proj_settings = list(n_proj = 2))))
    expect_equal(projected$rec$matrix[!projected$is_proj, ], historical$rec$matrix, ignore_attr = TRUE)
    expected <- model.matrix(historical$rec$terms,
      model.frame(historical$rec$terms, data = data.frame(year = 1991:1992)))
    expect_equal(projected$rec$matrix[projected$is_proj, ], expected, ignore_attr = TRUE)
  }
  expect_error(prepare_tam(cod_obs, years = 1983:1990, ages = 2:6,
    rec_settings = list(process = "off", mu_form = ~factor(year) + year),
    proj_settings = list(n_proj = 1)), "redundant")
})
