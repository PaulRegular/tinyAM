test_that("survey and catchability units are concise without inventing dimensions", {
  expect_equal(tinyAM:::.q_unit("thousand fish", "thousand fish"), "unitless")
  expect_equal(tinyAM:::.q_unit("thousand fish per tow", "thousand fish"), "per tow")
  expect_equal(tinyAM:::.q_unit("thousand fish / net", "thousand fish"), "per net")
  expect_true(is.na(tinyAM:::.q_unit("native survey-index units", "thousand fish")))
  expect_equal(tinyAM:::.unit_label(data.frame(est = 1, unit = "thousand fish"), "Abundance"),
    "Abundance (thousand fish)")
  expect_equal(tinyAM:::.unit_label(data.frame(est = 1, unit = "native survey-index units"), "Abundance index"),
    "Abundance index")
  expect_equal(tinyAM:::.count_unit("mean numbers per tow"), list(scale = 1, effort = "tow"))
})
