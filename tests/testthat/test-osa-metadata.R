test_that("tidying preserves retained OSA metadata across intervals and model stacks", {
  f <- default_fit
  local_mocked_bindings(oneStepPredict = function(obj, subset, ...) {
    data.frame(residual = subset / 100)
  }, .package = "RTMB")
  f <- add_osa_residuals(f, method = "oneStepGeneric", subset = seq_len(3L), trace = FALSE)
  f$obs_pred$catch <- f$obs_pred$catch[rev(seq_len(nrow(f$obs_pred$catch))), ]
  f$obs_pred$catch$custom_note <- as.character(seq_len(nrow(f$obs_pred$catch)))
  original <- f$obs_pred
  one <- tidy_tam(f, interval = .8, process_residuals = FALSE)
  expect_identical(one$osa$f, f$osa)
  expect_identical(attr(one$obs_pred, "osa")$f, f$osa)
  expect_identical(unique(one$obs_pred$catch$osa_method), "oneStepGeneric")
  expected <- .observation_match(one$obs_pred$catch, original$catch)
  expect_identical(one$obs_pred$catch$osa_res, original$catch$osa_res[expected])
  expect_identical(one$obs_pred$catch$custom_note, original$catch$custom_note[expected])
  expect_identical(f$obs_pred, original)
  both <- tidy_tam(Retained = f, Uncomputed = default_fit, process_residuals = FALSE)
  expect_identical(names(both$osa), "Retained")
  expect_identical(both$osa$Retained, f$osa)
  expect_identical(unique(both$obs_pred$catch$osa_method[both$obs_pred$catch$model == "Retained"]), "oneStepGeneric")
  expect_true(all(is.na(both$obs_pred$catch$osa_method[both$obs_pred$catch$model == "Uncomputed"])))
})
