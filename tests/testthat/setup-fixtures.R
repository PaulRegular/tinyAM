## Common fixtures shared across testthat files

YEARS <- 1983:2024
AGES  <- 2:14

if (!exists("cod_obs", inherits = TRUE)) {
  cod_obs <- tinyAM::cod_obs
}

make_test_dat <- function(...) {
  args <- utils::modifyList(
    list(
      obs = cod_obs,
      years = YEARS,
      ages = AGES
    ),
    list(...)
  )
  do.call(make_dat, args)
}

if (exists("cod_obs", inherits = TRUE)) {
  default_dat <- make_test_dat(
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ I(0.3))
  )
  default_par <- make_par(default_dat)
} else {
  default_dat <- NULL
  default_par <- NULL
}

if (requireNamespace("RTMB", quietly = TRUE) && exists("cod_obs", inherits = TRUE)) {
  set.seed(1)
  default_fit <- fit_tam(
    cod_obs,
    years = YEARS,
    ages = AGES,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
    silent = TRUE,
    grad_tol = 0.1
  )

  set.seed(1)
  N_dev <- update(
    default_fit,
    proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1),
    silent = TRUE
  )

  set.seed(1)
  M_dev <- update(
    N_dev,
    N_settings = list(process = "off", init_N0 = TRUE),
    M_settings = list(
      process = "ar1",
      mu_supplied = ~ I(0.3),
      age_breaks = c(3, 14)
    )
  )
} else {
  default_fit <- NULL
  N_dev <- NULL
  M_dev <- NULL
}
