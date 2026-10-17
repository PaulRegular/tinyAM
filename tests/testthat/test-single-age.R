test_that("single-age views are annual Plotly widgets without view controls", {
  d <- data.frame(year = 2000:2003, age = 2, est = 1:4)
  original <- d
  for (frame in list(~age, ~year)) {
    p <- plot_age_year(d, x = if (identical(all.vars(frame), "year")) ~age else ~year, frame = frame)
    expect_s3_class(p, "plotly")
    built <- plotly::plotly_build(p)
    expect_length(built$x$frames, 0)
    expect_equal(built$x$layout$title, "Age 2")
    labels <- unlist(lapply(p$x$layout$updatemenus, function(m) vapply(m$buttons, `[[`, "", "label")))
    expect_false(any(labels %in% c("Annual trends", "Age profiles")))
  }
  expect_identical(d, original)
})

test_that("one catchability estimate per model has points and switchable intervals", {
  d <- data.frame(age = 2, model = c("A", "B", "Accepted"),
    est = c(.2, .3, NA), lwr = c(.1, .2, NA), upr = c(.4, .5, NA))
  original <- d
  p <- suppressWarnings(plot_trend(d, x = ~age, color = ~model))
  traces <- suppressWarnings(plotly::plotly_build(p))$x$data
  expect_true(any(vapply(traces, function(x) identical(x$mode, "markers"), logical(1))))
  expect_true(any(vapply(traces, function(x) isTRUE(x$meta$tam_ci), logical(1))))
  expect_true(any(vapply(p$x$layout$updatemenus, function(m) any(vapply(m$buttons, `[[`, "", "label") == "Hide CI"), logical(1))))
  expect_identical(d, original)
})
