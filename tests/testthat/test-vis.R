## vis_tam ----

# Fit two quick models (as in examples)
N_dev <- fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  M_settings = list(process = "off", assumption = ~ I(0.3)),
  proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1),
  silent = TRUE
)
M_dev <- update(
  N_dev,
  N_settings = list(process = "off", init_N0 = TRUE),
  M_settings = list(process = "approx_rw", assumption = ~ I(0.3),
                    age_breaks = c(2, 14))
)
fits <- list("N_dev" = N_dev, "M_dev" = M_dev)


test_that("vis_tam renders cleanly and produces an HTML output", {
  # --- Test: render to temporary file quietly ---
  tmpfile <- tempfile(fileext = ".html")

  expect_no_error({
    vis_tam(fits, output_file = tmpfile, open_file = FALSE, quiet = TRUE)
  })

  expect_true(file.exists(tmpfile))
  expect_match(readLines(tmpfile, n = 1L), "<!DOCTYPE html>", fixed = TRUE)

  # --- Skip interactive browser behavior on CI / non-interactive sessions ---
  skip_if_not(interactive())
  skip_on_cran()
  skip_on_ci()

  expect_no_error({
    vis_tam(fits, output_file = NULL, open_file = TRUE, quiet = TRUE)
  })
})
