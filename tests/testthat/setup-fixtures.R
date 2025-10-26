## Common fixtures shared across testthat files

YEARS <- 1983:2024
AGES  <- 2:14

make_test_dat <- function(...) {
  dots <- list(...)
  if (!"obs" %in% names(dots)) {
    if (!exists("cod_obs", inherits = TRUE)) {
      stop("cod_obs not available; supply `obs` explicitly", call. = FALSE)
    }
    dots$obs <- cod_obs
  }
  args <- utils::modifyList(list(years = YEARS, ages = AGES), dots)
  do.call(make_dat, args)
}

if (exists("cod_obs", inherits = TRUE)) {
  default_dat <- make_test_dat(
    N_settings = list(process = "iid", init_N0 = TRUE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3))
  )
  default_par <- make_par(default_dat)
} else {
  default_dat <- NULL
  default_par <- NULL
}

.base_fit_args <- list(
  obs = NULL,
  years = YEARS,
  ages = AGES,
  N_settings = list(process = "iid", init_N0 = FALSE),
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", mu_form = NULL, assumption = ~ I(0.3)),
  obs_settings = list(
    q_form = ~ q_block,
    sd_catch_form = ~ 1,
    sd_index_form = ~ 1,
    fill_missing = TRUE
  ),
  proj_settings = NULL,
  silent = TRUE
)

if (exists("cod_obs", inherits = TRUE)) {
  .base_fit_args$obs <- cod_obs
}

make_test_fit <- function(...) {
  if (!requireNamespace("RTMB", quietly = TRUE)) {
    testthat::skip("RTMB not installed")
  }
  args <- utils::modifyList(.base_fit_args, list(...))
  if (is.null(args$obs)) {
    stop("`obs` must be supplied when cod_obs is unavailable", call. = FALSE)
  }
  do.call(fit_tam, args)
}

if (requireNamespace("RTMB", quietly = TRUE) && !is.null(.base_fit_args$obs)) {
  set.seed(1)
  default_fit <- do.call(fit_tam, .base_fit_args)
} else {
  default_fit <- NULL
}
