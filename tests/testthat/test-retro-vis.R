test_that("S3 dispatch preserves legacy fit, reference and list labels", {
  testthat::local_mocked_bindings(.render_vis_tam = function(fits, ...) names(fits),
    .package = "tinyAM")
  first <- second <- default_fit
  reference <- default_fit; class(reference) <- c("tam_ref", "list")
  expect_identical(vis_tam(first), "first")
  expect_identical(vis_tam(first, other = second), c("first", "other"))
  expect_identical(vis_tam(reference, first), c("reference", "first"))
  expect_identical(vis_tam(list(a = first, b = second)), c("a", "b"))
  expect_identical(vis_tam(model_list = list(a = first)), "a")
})

test_that("dashboard dispatch retains reproducible process controls", {
  testthat::local_mocked_bindings(.render_vis_tam = function(fits, ...,
    process_seed = 123456L, retro_origins = NULL) {
      list(seed = process_seed, origins = retro_origins)
    }, .package = "tinyAM")
  r <- structure(list(fits = list(base = list(`2024` = default_fit)),
    fold_diagnostics = data.frame(model = "base", fold = 2024L,
      is_reference = TRUE, is_converged = TRUE, error = NA_character_)),
    class = "tam_retro")
  expect_identical(vis_tam(default_fit, process_seed = 42L)$seed, 42L)
  expect_identical(vis_tam(model_list = list(base = default_fit), process_seed = 43L)$seed, 43L)
  view <- vis_tam(r, process_seed = 44L)
  expect_identical(view$seed, 44L)
  expect_identical(view$origins, c(`2024` = 2024L))
})

test_that("retrospective model selection, omissions and ordering are explicit", {
  r <- structure(list(fits = list(base = list(`2022` = default_fit,
    `2024` = default_fit, `2023` = NULL)), fold_diagnostics = data.frame(
      model = "base", fold = 2022:2024, is_reference = c(FALSE, FALSE, TRUE),
      is_converged = c(TRUE, FALSE, TRUE), error = c(NA, "failed fold", NA))), class = "tam_retro")
  original <- r
  view <- tinyAM:::.retro_dashboard(r)
  expect_identical(names(view$fits), c("2024", "2022"))
  expect_match(paste(view$background, collapse = " "), "Full-data reference year: \\*\\*2024")
  expect_match(paste(view$background, collapse = " "), "2023: failed fold")
  expect_identical(r, original)
  missing_reference <- r
  missing_reference$fits$base$`2024` <- NULL
  without_reference <- tinyAM:::.retro_dashboard(missing_reference)
  expect_identical(names(without_reference$fits), "2022")
  expect_match(paste(without_reference$background, collapse = " "), "reference is unavailable")
  r$fits$alternative <- r$fits$base
  expect_error(tinyAM:::.retro_dashboard(r, is_interactive = FALSE), "base, alternative")
  expect_null(suppressMessages(tinyAM:::.retro_dashboard(r, is_interactive = TRUE,
    choose = function(...) "")))
  expect_identical(suppressMessages(tinyAM:::.retro_dashboard(r, is_interactive = TRUE,
    choose = function(...) "base"))$origins, view$origins)
  expect_error(tinyAM:::.retro_dashboard(r, "absent"), "exact name")
  expect_error(vis_tam(r, default_fit), "alone")
  expect_error(vis_tam(default_fit, r), "alone")
  expect_error(vis_tam(r, model_list = list(base = default_fit)), "alone")
  r$fold_diagnostics$is_converged <- FALSE
  expect_error(tinyAM:::.retro_dashboard(r, "base"), "fold_diagnostics")
})

test_that("forecast metadata are public plot defaults and mixed origins omit boundaries", {
  d <- data.frame(year = rep(2020:2023, 2), fold = rep(c(2021, 2022), each = 4),
    model = rep(c("a", "b"), each = 4), age = 2, est = 1:8, obs = 1:8, pred = 1:8)
  d$horizon <- pmax(0, d$year - d$fold)
  d$is_proj <- d$horizon > 0
  original <- d
  for (plot in list(plot_trend, plot_obs_pred)) {
    p <- plotly::plotly_build(plot(d, color = ~model, add_buttons = FALSE))
    expect_length(p$x$layout$shapes, 0)
    expect_match(paste(unlist(lapply(p$x$data, `[[`, "text")), collapse = " "), "Origin: 2021")
    expect_match(paste(unlist(lapply(p$x$data, `[[`, "text")), collapse = " "), "Forecast; horizon: 1")
    p <- plotly::plotly_build(plot(d[d$model == "a", ], add_buttons = FALSE))
    expect_equal(p$x$layout$shapes[[1]]$x0, 2021)
    p <- plotly::plotly_build(plot(d, text = "custom", hovertemplate = "%{text}", add_buttons = FALSE))
    expect_false(any(grepl("Origin", unlist(lapply(p$x$data, `[[`, "text")))))
  }
  expect_identical(d, original)
  tabs <- tinyAM:::.retro_plot_metadata(list(N = d), c(a = 2021, b = 2022))
  expect_equal(tabs$N$fold, d$fold)
})
