test_that("retrospective and hindcast objects render without changing their fits", {
  r <- fit_retro(default_fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  original <- r$fold_diagnostics
  parameters <- lapply(r$fits[[1]], function(f) f$parameter_values)
  # RTMB environments are shared by R copies; capture their values explicitly.
  objective_state <- function(f) lapply(c("last.par", "last.par.best"),
    function(nm) get(nm, envir = f$obj$env, inherits = FALSE))
  states <- lapply(r$fits[[1]], objective_state)
  file <- tempfile(fileext = ".html")
  vis_tam(r, output_file = file, open_file = FALSE, render_args = list(quiet = TRUE),
    background = "User background retained.")
  expect_true(file.exists(file))
  con <- file(file, "rb"); html <- paste(readLines(con, warn = FALSE), collapse = "\n"); close(con)
  expect_match(html, "Retrospective comparison", fixed = TRUE)
  expect_match(html, "User background retained.", fixed = TRUE)
  expect_match(html, "Full-data[[:space:]]+reference year:")
  expect_false(grepl("(full data)", html, fixed = TRUE))
  expect_false(grepl("(peel)", html, fixed = TRUE))
  expect_identical(r$fold_diagnostics, original)
  expect_identical(lapply(r$fits[[1]], function(f) f$parameter_values), parameters)
  expect_identical(lapply(r$fits[[1]], objective_state), states)
  h <- fit_hindcast(default_fit, folds = 1, start_from_fit = TRUE, progress = FALSE,
    proj_settings = list(F_mult = 1))
  h$fits$alternative <- h$fits[[1]]
  file <- tempfile(fileext = ".html")
  vis_tam(h, model = names(h$fits)[1], output_file = file, open_file = FALSE,
    render_args = list(quiet = TRUE))
  con <- file(file, "rb"); html <- paste(readLines(con, warn = FALSE), collapse = "\n"); close(con)
  expect_match(html, "Origin:", fixed = TRUE)
  expect_match(html, "Forecast; horizon:", fixed = TRUE)
})
