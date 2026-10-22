test_that("development hindcasts refresh stale worker namespaces", {
  skip_if_not_installed("pkgload")
  skip_if_not(pkgload::is_dev_package("tinyAM"), "Requires the development-source workflow")
  previous <- future::plan()
  withr::defer(future::plan(previous))
  future::plan(future::multisession, workers = I(1))

  worker_pid <- future::value(future::future({
    ns <- asNamespace("tinyAM")
    if (exists(".retro_fold", ns, inherits = FALSE)) {
      assignInNamespace(".retro_fold", function(...) stop("Stale worker code"), ns = "tinyAM")
    }
    Sys.getpid()
  }))
  expect_false(identical(worker_pid, Sys.getpid()))
  before <- default_fit$parameter_values
  result <- fit_hindcast(default_fit, folds = 0, progress = FALSE)
  expect_true(all(result$fold_diagnostics$is_converged))
  expect_true(all(is.na(result$fold_diagnostics$error)))
  expect_identical(default_fit$parameter_values, before)
})
