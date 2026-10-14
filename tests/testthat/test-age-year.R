dashboard_helpers <- local({
  lines <- readLines(system.file("rmd", "vis_tam.Rmd", package = "tinyAM"))
  first <- match("```{r helpers, include = FALSE}", lines) + 1L
  last <- first + match("```", lines[first:length(lines)]) - 2L
  env <- new.env()
  # Evaluate helper definitions without fitting or rendering a dashboard.
  env$fits <- list()
  eval(parse(text = lines[first:last]), envir = env)
  env
})

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
  annual <- tinyAM::plot_age_year(d, plot = plot, frame = ~age, color = ~model)
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
  profile <- tinyAM::plot_age_year(d, x = ~age, xlab = "Age", frame = ~year)
  expect_match(htmltools::renderTags(profile)$html, 'data-default-view="profile"', fixed = TRUE)
})

test_that("both views preserve ribbons, model groups, and projection semantics", {
  plots <- list()
  plot <- function(...) {
    p <- plot_trend(...)
    plots[[length(plots) + 1L]] <<- plotly::plotly_build(p)
    p
  }
  tinyAM::plot_age_year(age_year_data(), plot = plot, frame = ~age,
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
  expect_s3_class(tinyAM::plot_age_year(d, frame = NULL), "plotly")
  d$age_block <- d$age
  expect_s3_class(tinyAM::plot_age_year(d, frame = ~age_block), "plotly")
  for (rows in list(d[d$age == 1, ], d[d$year == 2000, ])) {
    expect_no_error(htmltools::renderTags(tinyAM::plot_age_year(rows, frame = ~age)))
  }
  expect_s3_class(tinyAM::plot_age_year(d), "shiny.tag")
})

test_that("exported plotting functions handle unavailable reporting tables", {
  for (d in list(NULL, data.frame(), transform(age_year_data(), est = NA_real_))) {
    expect_s3_class(plot_trend(d), "plotly")
    expect_no_warning(plotly::plotly_build(plot_trend(d)))
    expect_no_warning(htmltools::renderTags(tinyAM::plot_age_year(d)))
  }
  d <- transform(age_year_data(), obs = NA_real_, pred = NA_real_)
  for (rows in list(NULL, d[0, ], d)) {
    expect_no_warning(plotly::plotly_build(plot_obs_pred(rows)))
  }
})

test_that("exported controls include ribbons added after the original trend", {
  d <- age_year_data()
  p <- plot_trend(d, frame = ~age, color = ~model, add_intervals = FALSE,
                  add_buttons = FALSE)
  p <- plotly::add_ribbons(p, ymin = ~lwr, ymax = ~upr, showlegend = FALSE)
  p <- plotly::plotly_build(tinyAM::add_plot_buttons(p))
  ribbons <- which(vapply(p$x$data, function(trace) !is.null(trace$fill), logical(1))) - 1L
  menu <- p$x$layout$updatemenus[[2]]
  expect_identical(menu$buttons[[2]]$label, "Hide CI")
  expect_equal(menu$buttons[[2]]$args[[2]], ribbons)
  expect_length(p$x$frames, 3L)
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

test_that("mortality means retain ones without borrowing annual uncertainty", {
  d <- age_year_data()
  mean <- transform(d, est = 1, lwr = .8, upr = 1.2)
  fits <- list(A = structure(list(dat = list(F_settings = list(mu_form = ~1))), class = "tam_fit"),
               B = structure(list(dat = list(F_settings = list(mu_form = NULL))), class = "tam_fit"))
  joined <- dashboard_helpers$mortality_plot_data(d, mean, fits, "F")
  expect_identical(joined$est, d$est)
  expect_identical(joined$lwr, d$lwr)
  expect_identical(joined$upr, d$upr)
  expect_true(all(joined$mean_est[joined$model == "A"] == 1))
  expect_true(all(is.na(joined$mean_est[joined$model == "B"])))
  joined$lwr <- joined$upr <- NA_real_
  p <- plotly::plotly_build(dashboard_helpers$plot_mortality(joined, frame = ~age,
    color = ~model, legendgroup = ~model))
  ribbons <- Filter(function(tr) !is.null(tr$fill), p$x$data)
  expect_gt(length(ribbons), 0L)
  expect_true(all(vapply(ribbons, function(tr) all(is.na(tr$y) | tr$y <= 1.2), logical(1))))
  mean_lines <- Filter(function(tr) !is.null(tr$hovertemplate), p$x$data)
  expect_gt(length(mean_lines), 0L)
  expect_true(all(vapply(mean_lines, function(tr) identical(tr$line$dash, "dot"), logical(1))))
  custom <- plotly::plotly_build(dashboard_helpers$plot_mortality(joined, frame = ~age,
    color = ~model, mean_line = list(dash = "dash", width = 2)))
  custom_lines <- Filter(function(tr) !is.null(tr$hovertemplate), custom$x$data)
  expect_true(all(vapply(custom_lines, function(tr)
    identical(tr$line$dash, "dash") && tr$line$width == 2, logical(1))))
  expect_match(paste(unlist(p$x$layout$updatemenus), collapse = " "), "Hide CI", fixed = TRUE)
})

test_that("sparse reference surfaces retain model traces in both orientations", {
  d <- age_year_data()
  d$est[d$model == "A" & (d$age == 2 | d$year == 2001)] <- NA_real_
  d$lwr <- d$upr <- NA_real_
  original <- d
  plots <- list()
  plot <- function(...) {
    p <- plot_trend(...)
    plots[[length(plots) + 1L]] <<- plotly::plotly_build(p)
    p
  }
  tinyAM::plot_age_year(d, plot = plot, frame = ~age, color = ~model, legendgroup = ~model)
  expect_identical(d, original)
  for (p in plots) {
    expect_length(p$x$data, 2L)
    for (frame in p$x$frames) {
      expect_equal(frame$traces, 0:1)
      expect_length(frame$data, 2L)
    }
  }
})

test_that("different projection horizons do not split a common year into frames", {
  d <- age_year_data()
  d$is_proj[d$model == "A"] <- FALSE
  views <- list()
  plot <- function(data, ...) {
    views[[length(views) + 1L]] <<- data
    plot_trend(data, ...)
  }
  tinyAM::plot_age_year(d, plot = plot, frame = ~age, color = ~model)
  expect_equal(levels(views[[2]]$.view_year), c("2000", "2001", "2002 (includes projections)"))
  expect_equal(length(unique(views[[2]]$.view_year[views[[2]]$year == 2002])), 1L)
  expect_match(views[[2]]$.view_hover[views[[2]]$is_proj][1], "(projection)", fixed = TRUE)
  expect_false(any(grepl("(projection)", views[[2]]$.view_hover[!views[[2]]$is_proj], fixed = TRUE)))
})

test_that("view buttons inherit native Plotly styling and control spacing", {
  d <- age_year_data()
  for (intervals in c(TRUE, FALSE)) {
    p <- plot_trend(d, frame = ~age, add_intervals = intervals)
    native <- plotly::plotly_build(p)$x$layout$updatemenus
    for (view in c("annual", "profile")) {
      menus <- plotly::plotly_build(tinyAM:::.age_year_menu(p, view))$x$layout$updatemenus
      menu <- menus[[length(menus)]]
      controls <- Filter(function(menu) any(unlist(lapply(menu$buttons, `[[`, "label")) %in%
        c("Linear", "Log", "Show CI", "Hide CI")), native)
      previous <- controls[[length(controls)]]
      expect_length(menus, length(native) + 1L)
      expect_identical(menu[c("type", "x", "y", "xanchor", "yanchor")],
                       previous[c("type", "x", "y", "xanchor", "yanchor")])
      expect_identical(menu$pad$r, previous$pad$r)
      expect_equal(menu$pad$t, if (intervals) 150 else 75)
      expect_equal(menu$active, if (view == "annual") 0 else 1)
      expect_identical(vapply(menu$buttons, `[[`, character(1), "label"),
                       c("Annual trends", "Age profiles"))
      expect_true(all(vapply(menu$buttons, function(button) identical(button$execute, FALSE), logical(1))))
    }
  }
  d$obs <- d$est
  d$pred <- d$est
  menus <- plotly::plotly_build(tinyAM:::.age_year_menu(plot_obs_pred(d), "annual"))$x$layout$updatemenus
  expect_length(menus, 2L)
  expect_equal(menus[[2]]$pad$t, 75)
})
