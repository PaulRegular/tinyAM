test_that("ribbons close within valid interval runs", {
  d <- expand.grid(year = 2000:2006, model = c("Official", "tinyAM"))
  d$est <- seq_len(nrow(d))
  d$lwr <- d$est - .5
  d$upr <- d$est + .5
  d[d$model == "Official" & d$year %in% c(2003, 2006), c("lwr", "upr")] <- NA
  original <- d
  p <- plotly::plotly_build(plot_trend(d[nrow(d):1, ], color = ~model))
  ribbons <- Filter(function(z) identical(z$fill, "toself"), p$x$data)
  expect_length(ribbons, 3L)
  expect_true(all(vapply(ribbons, function(z) z$legendgroup %in% levels(d$model), logical(1))))
  expect_true(all(vapply(ribbons, function(z) all(is.finite(z$x) & is.finite(z$y)), logical(1))))
  expect_true(all(vapply(ribbons, function(z) !all(c(2002, 2004) %in% z$x) ||
    min(z$y) > 7, logical(1))))
  expect_identical(d, original)
})

test_that("signed trends include estimates and both interval endpoints", {
  d <- data.frame(year = 2000:2002, est = c(-2, 0, 1), lwr = c(-3, -.2, .5), upr = c(-1, .2, 2))
  p <- plotly::plotly_build(plot_trend(d))
  expect_lt(p$x$layout$yaxis$range[1], -3)
  expect_gt(p$x$layout$yaxis$range[2], 2)
  d[c("est", "lwr", "upr")] <- lapply(d[c("est", "lwr", "upr")], abs)
  expect_equal(plotly::plotly_build(plot_trend(d))$x$layout$yaxis$range[1], 0)
})
