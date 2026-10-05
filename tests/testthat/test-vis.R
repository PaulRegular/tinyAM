
## vis_tam ----

suppressWarnings(
  source(system.file("examples/example_fits.R", package = "tinyAM"))
)

dashboard_plots <- function(file) {
  html <- paste(readLines(file, warn = FALSE), collapse = "\n")
  scripts <- regmatches(html, gregexpr(
    '(?s)<script type="application/json"[^>]*>.*?</script>', html, perl = TRUE
  ))[[1L]]
  widgets <- lapply(scripts, function(script) {
    json <- sub("^<script[^>]*>", "", script)
    json <- sub("</script>$", "", json)
    jsonlite::fromJSON(json, simplifyVector = FALSE)
  })
  lapply(Filter(function(widget) !is.null(widget$x$layout), widgets), `[[`, "x")
}

dashboard_residuals <- function(plots) {
  Filter(function(plot) {
    title <- unlist(c(plot$layout$title, plot$layout$yaxis$title))
    any(grepl("residuals", title, fixed = TRUE))
  }, plots)
}

test_that("vis_tam renders cleanly and produces an HTML output", {
  # --- Test: render to temporary file quietly ---
  tmpfile <- tempfile(fileext = ".html")

  expect_no_error({
    vis_tam(fits, output_file = tmpfile, open_file = FALSE,
            render_args = list(quiet = TRUE))
  })

  expect_true(file.exists(tmpfile))
  expect_match(readLines(tmpfile, n = 1L), "<!DOCTYPE html>", fixed = TRUE)

  residuals <- dashboard_residuals(dashboard_plots(tmpfile))
  expect_gt(length(residuals), 0L)
  for (plot in residuals) {
    expect_gt(length(plot$data), 0L)
    expect_true(any(is.finite(unlist(lapply(plot$data, `[[`, "y")))))
    frame_names <- vapply(plot$frames, `[[`, character(1), "name")
    expect_setequal(frame_names, names(fits))
    expect_true(all(vapply(plot$frames, function(frame) length(frame$traces) > 0L,
                           logical(1))))
  }

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

test_that("vis_tam errors when supplied objects are not TAM fits", {
  fake <- list(bad = list(not = "a fit"))
  expect_error(
    vis_tam(model_list = fake,
            output_file = tempfile(fileext = ".html"), open_file = FALSE,
            render_args = list(quiet = TRUE)),
    "All supplied models must be"
  )
})


test_that("tam_ref objects accept a Background page", {
  blank_values <- function(x, keys) {
    for (name in setdiff(names(x), keys)) x[[name]][] <- NA
    x
  }
  source_pop <- lapply(N_dev$pop, blank_values,
                       keys = c("year", "age", "is_proj"))
  source_pop$ssb$est <- N_dev$pop$ssb$est + 1
  source_pop$N$est <- N_dev$pop$N$est + 1
  source_pop$F$est <- N_dev$pop$F$est + 0.1
  source_pop$ssb$unit <- "t"
  source_pop$N$unit <- "thousand fish"
  source_pop$recruitment$unit <- "thousand fish"
  source_pop$recruitment$notes <- rep("Recruitment ages 3\u20138",
                                      nrow(source_pop$recruitment))
  source_obs_pred <- N_dev$obs_pred
  source_obs_pred <- lapply(source_obs_pred, blank_values,
                            keys = c("year", "age", "survey", "fleet", "samp_time",
                                     "q_block", "q_key", "is_proj", "obs"))
  source_obs_pred <- lapply(source_obs_pred, function(x) { x$osa_res <- NA_real_; x })
  source <- N_dev
  source$call <- quote(database_to_tam_ref("fixture"))
  source$pop <- source_pop
  source$obs_pred <- source_obs_pred
  source$fixed_par <- blank_values(N_dev$fixed_par, c("par", "coef", "age"))
  source$random_par <- lapply(N_dev$random_par, blank_values,
                             keys = c("par", "coef", "age", "year", "is_proj"))
  source$rep <- lapply(N_dev$rep, function(x) { x[] <- NA; x })
  source$is_converged <- NA
  source$grad_tol <- NA_real_
  source$comparison_scales <- c(ssb = 1e-3, N = 1e-3, recruitment = 1e-3)
  class(source) <- c("tam_ref", "list")

  expect_equal(names(source$pop), names(N_dev$pop))
  expect_equal(names(source$rep), names(N_dev$rep))
  expect_equal(names(source$obs_pred), names(N_dev$obs_pred))
  expect_equal(source$random_par$log_f$year, N_dev$random_par$log_f$year)
  expect_true(all(is.na(source$random_par$log_f$est)))
  expect_equal(source$pop$N[c("year", "age", "is_proj")],
               N_dev$pop$N[c("year", "age", "is_proj")])
  expect_true(all(is.finite(source$pop$N$est)))

  tabs <- tidy_tam(model_list = list(Assessment = source, tinyAM = N_dev))
  expect_true(all(c("Assessment", "tinyAM") %in% tabs$pop$ssb$model))
  expect_equal(tabs$pop$ssb$est[tabs$pop$ssb$model == "Assessment"], source_pop$ssb$est)
  expect_equal(tabs$pop$ssb$est[tabs$pop$ssb$model == "tinyAM"],
               N_dev$pop$ssb$est * 1e-3)
  expect_equal(unique(tabs$pop$ssb$unit), "t")
  ssb_plot <- plot_trend(tabs$pop$ssb, color = ~model, add_buttons = FALSE)
  ssb_traces <- plotly::plotly_build(ssb_plot)$x$data
  accepted_ssb <- Filter(function(trace) identical(trace$mode, "lines") &&
                           identical(trace$name, "Assessment") && is.null(trace$fill),
                         ssb_traces)
  expect_length(accepted_ssb, 1L)
  expect_equal(as.numeric(accepted_ssb[[1L]]$y), source_pop$ssb$est)

  file <- tempfile(fileext = ".html")
  vis_tam(model_list = list(Assessment = source, tinyAM = N_dev),
          background = c("## Assessment assumptions", "",
                         "| Component | Accepted assessment | tinyAM representation | Reason |",
                         "|---|---|---|---|",
                         "| Years | 2000\u20132024 and ages are retained. | Same years | Common period |"),
          output_file = file, open_file = FALSE,
          render_args = list(quiet = TRUE))
  con <- file(file, "rb")
  html <- readLines(con, warn = FALSE)
  close(con)
  expect_true(any(grepl("Background", html, fixed = TRUE)))
  expect_true(any(grepl("Assessment assumptions", html, fixed = TRUE)))
  expect_true(any(grepl("ages are retained", html, fixed = TRUE)))
  expect_true(any(grepl('id="assessment-background"', html, fixed = TRUE)))
  expect_true(any(grepl("<th>Component</th>", html, fixed = TRUE)))
  expect_true(any(grepl("<td>Years</td>", html, fixed = TRUE)))
  expect_true(any(grepl("<h1>Fishery</h1>", html, fixed = TRUE)))
  expect_true(any(grepl("<h1>Trends</h1>", html, fixed = TRUE)))
  expect_true(any(grepl("<h1>Parameters</h1>", html, fixed = TRUE)))
  expect_true(any(grepl("<h1>Inputs</h1>", html, fixed = TRUE)))
  expect_true(any(grepl("Spawning stock biomass", html, fixed = TRUE)))
  expect_true(any(grepl("Average F", html, fixed = TRUE)))
  expect_true(any(grepl("Average M", html, fixed = TRUE)))
  expect_true(any(grepl("Metrics-at-age", html, fixed = TRUE)))
  for (field in c("N", "F", "M")) {
    expect_true(any(grepl(paste0(">", field, "</th>"), html, fixed = TRUE)))
  }
  expect_false(any(grepl(".x</th>", html, fixed = TRUE)))
  expect_false(any(grepl(".y</th>", html, fixed = TRUE)))
  expect_true(any(grepl("Recruitment ages 3", html, fixed = TRUE)))
  expect_false(any(grepl("Comparing SAM and tinyAM", html, fixed = TRUE)))
  visible_html <- paste(html, collapse = "\n")
  visible_html <- gsub("(?is)<(script|style)\\b.*?</\\1>", "", visible_html,
                        perl = TRUE)
  expect_false(grepl("unavailable", tolower(visible_html), fixed = TRUE))

  plots <- dashboard_plots(file)
  residuals <- dashboard_residuals(plots)
  expect_gt(length(residuals), 0L)
  for (plot in residuals) {
    expect_gt(length(plot$data), 0L)
    expect_true(any(is.finite(unlist(lapply(plot$data, `[[`, "y")))))
    expect_length(plot$frames, 0L)
  }
  for (metric in c("N", "F")) {
    comparisons <- Filter(function(plot) {
      identical(unname(unlist(plot$layout$yaxis$title)), metric)
    }, plots)
    expect_gt(length(comparisons), 0L)
    for (plot in comparisons) {
      traces <- c(plot$data, unlist(lapply(plot$frames, `[[`, "data"), recursive = FALSE))
      reference <- Filter(function(trace) identical(trace$name, "Assessment") &&
                            identical(trace$mode, "lines"), traces)
      expect_gt(length(reference), 0L)
      expect_true(any(is.finite(unlist(lapply(reference, `[[`, "y")))))
    }
  }
})
