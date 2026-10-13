age_year_data <- function() {
  d <- expand.grid(year = 2000:2002, age = 1:3, model = c("A", "B"))
  d$est <- d$year - 1999 + d$age
  d$lwr <- d$est - .5
  d$upr <- d$est + .5
  d$is_proj <- d$year == 2002
  d
}

test_that("age-year views honor the requested frame and preserve inputs", {
  d <- age_year_data()[18:1, ]
  original <- d
  calls <- list()
  plot <- function(data, x, xlab, frame, ...) {
    calls[[length(calls) + 1L]] <<- list(data = data, x = x, xlab = xlab, frame = frame)
    plot_trend(data, x = x, xlab = xlab, frame = frame, ...)
  }
  annual <- tinyAM:::.plot_age_year(d, plot = plot, frame = ~age, color = ~model)
  expect_identical(d, original)
  expect_length(calls, 2L)
  expect_identical(all.vars(calls[[1]]$x), "year")
  expect_identical(all.vars(calls[[1]]$frame), "age")
  expect_identical(all.vars(calls[[2]]$x), "age")
  expect_identical(all.vars(calls[[2]]$frame), ".view_year")
  expect_identical(calls[[1]]$data$year, sort(d$year))
  expect_equal(levels(calls[[2]]$data$.view_year),
               c("2000", "2001", "2002 (projection)"))
  html <- htmltools::renderTags(annual)$html
  expect_match(html, 'data-default-view="annual"', fixed = TRUE)
  expect_match(html, "Annual trends", fixed = TRUE)
  expect_match(html, "Age profiles", fixed = TRUE)
  profile <- tinyAM:::.plot_age_year(d, x = ~age, xlab = "Age", frame = ~year)
  expect_match(htmltools::renderTags(profile)$html, 'data-default-view="profile"', fixed = TRUE)
})

test_that("both views preserve ribbons, model groups, and projection semantics", {
  plots <- list()
  plot <- function(...) {
    p <- plot_trend(...)
    plots[[length(plots) + 1L]] <<- plotly::plotly_build(p)
    p
  }
  tinyAM:::.plot_age_year(age_year_data(), plot = plot, frame = ~age,
                        color = ~model, legendgroup = ~model)
  expect_length(plots, 2L)
  expect_length(plots[[1]]$x$layout$shapes, 1L)
  expect_length(plots[[2]]$x$layout$shapes, 0L)
  for (p in plots) {
    expect_length(p$x$frames, 3L)
    expect_true(any(vapply(p$x$data, function(tr) !is.null(tr$fill), logical(1))))
    expect_match(paste(unlist(p$x$layout$updatemenus), collapse = " "), "Hide CI", fixed = TRUE)
    expect_setequal(unique(vapply(p$x$data, function(tr) as.character(tr$legendgroup),
                                 character(1))), c("A", "B"))
  }
})

test_that("unframed and block-indexed plots are left alone", {
  d <- age_year_data()
  expect_s3_class(tinyAM:::.plot_age_year(d), "plotly")
  d$age_block <- d$age
  expect_s3_class(tinyAM:::.plot_age_year(d, frame = ~age_block), "plotly")
  for (rows in list(d[d$age == 1, ], d[d$year == 2000, ])) {
    expect_no_error(htmltools::renderTags(tinyAM:::.plot_age_year(rows, frame = ~age)))
  }
})

test_that("observed-predicted plots accept an age axis without changing defaults", {
  d <- age_year_data()
  d$obs <- d$est
  d$pred <- d$est + .1
  original <- d
  annual <- plotly::plotly_build(plot_obs_pred(d, frame = ~age))
  profile <- plotly::plotly_build(plot_obs_pred(d, x = ~age, xlab = "Age", frame = ~year))
  expect_identical(d, original)
  expect_equal(annual$x$layout$xaxis$title, "Year")
  expect_equal(profile$x$layout$xaxis$title, "Age")
  expect_length(annual$x$layout$shapes, 1L)
  expect_length(profile$x$layout$shapes, 0L)
  for (p in list(annual, profile)) {
    expect_true(any(vapply(p$x$data, function(tr) tr$mode == "markers", logical(1))))
    expect_true(any(vapply(p$x$data, function(tr) tr$mode == "lines", logical(1))))
  }
})
