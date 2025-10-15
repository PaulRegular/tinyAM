
test_that("vis_tam runs cleanly with valid fits list", {

  # Fit two quick models (as in examples)
  N_dev <- fit_tam(
    cod_obs, years = 1983:2024, ages = 2:14,
    N_settings = list(process = "iid", init_N0 = FALSE),
    M_settings = list(process = "off", assumption = ~ I(0.3)),
    silent = TRUE
  )
  M_dev <- update(
    N_dev,
    N_settings = list(process = "off", init_N0 = TRUE),
    M_settings = list(process = "ar1", assumption = ~ I(0.3),
                      age_breaks = seq(2, 14, by = 6))
  )
  fits <- list("N_dev" = N_dev, "M_dev" = M_dev)

  # --- Test 1: rmarkdown::render (static mode)
  tmpfile <- tempfile(fileext = ".html")
  expect_no_error({
    vis_tam(fits, output_file = tmpfile, quiet = TRUE)
  })
  expect_true(file.exists(tmpfile))

  # --- Test 2: rmarkdown::run (interactive mode)
  skip_on_ci()
  skip_on_cran()
  expect_no_error({
    vis_tam(fits, output_file = NULL, shiny_args = list(quiet = TRUE))
  })

})

