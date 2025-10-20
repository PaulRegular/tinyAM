
## vis_tam ----

source(system.file("examples/example_fits.R", package = "tinyAM"))

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
