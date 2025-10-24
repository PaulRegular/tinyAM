
## vis_tam ----

source(system.file("examples/example_fits.R", package = "tinyAM"))

test_that("vis_tam renders cleanly and produces an HTML output", {
  # --- Test: render to temporary file quietly ---
  tmpfile <- tempfile(fileext = ".html")

  expect_no_error({
    vis_tam(fits, output_file = tmpfile, open_file = FALSE,
            render_args = list(quiet = TRUE))
  })

  expect_true(file.exists(tmpfile))
  expect_match(readLines(tmpfile, n = 1L), "<!DOCTYPE html>", fixed = TRUE)

  # --- Skip interactive browser behavior on CI / non-interactive sessions ---
  skip_if_not(interactive())
  skip_on_cran()
  skip_on_ci()

  expect_no_error({
    vis_tam(fits, output_file = NULL, open_file = TRUE,
            render_args = list(quiet = TRUE))
  })
})

test_that("render_args must be a named list", {
  expect_error(
    vis_tam(fits, output_file = tempfile(fileext = ".html"), open_file = FALSE,
            render_args = TRUE),
    "must be a list"
  )

  expect_error(
    vis_tam(fits, output_file = tempfile(fileext = ".html"), open_file = FALSE,
            render_args = list(TRUE)),
    "named list"
  )
})

test_that("vis_tam accepts models provided via dots", {
  tmpfile <- tempfile(fileext = ".html")

  expect_no_error({
    vis_tam(N_dev, M_dev, output_file = tmpfile, open_file = FALSE,
            render_args = list(quiet = TRUE))
  })

  expect_true(file.exists(tmpfile))
})

test_that("vis_tam accepts list input supplied through dots", {
  tmpfile <- tempfile(fileext = ".html")

  expect_no_error({
    vis_tam(list(N_dev = N_dev, M_dev = M_dev),
            output_file = tmpfile, open_file = FALSE,
            render_args = list(quiet = TRUE))
  })

  expect_true(file.exists(tmpfile))
})

test_that("vis_tam rejects unnamed or duplicated model lists", {
  unnamed <- list(N_dev, M_dev)
  duped <- structure(list(N_dev, M_dev), names = c("N_dev", "N_dev"))

  expect_error(
    vis_tam(model_list = unnamed, output_file = tempfile(fileext = ".html"),
            open_file = FALSE, render_args = list(quiet = TRUE)),
    "named list"
  )

  expect_error(
    vis_tam(model_list = duped, output_file = tempfile(fileext = ".html"),
            open_file = FALSE, render_args = list(quiet = TRUE)),
    "unique"
  )
})

test_that("vis_tam rejects mixed model inputs", {
  expect_error(
    vis_tam(N_dev, model_list = fits,
            output_file = tempfile(fileext = ".html"), open_file = FALSE,
            render_args = list(quiet = TRUE)),
    "Supply models either"
  )
})
