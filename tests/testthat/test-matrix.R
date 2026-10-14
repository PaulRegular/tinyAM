matrix_data <- function() {
  d <- expand.grid(year = 2000:2002, age = 1:3, model = c("tinyAM", "Accepted"),
                   metric = c("Abundance", "Biomass"))
  d$est <- seq_len(nrow(d))
  d$unit <- ifelse(d$metric == "Abundance", "fish", "kg")
  d$is_proj <- d$year == 2002
  d
}
matrix_spec <- function(p) p$jsHooks$render[[1L]]$data

test_that("matrix coordinates, ranges and metadata preserve sparse inputs", {
  d <- matrix_data()[c(30:1, 33:36), ]
  d$est[1] <- 0
  d$est[2] <- NA_real_
  d$definition <- "Source definition"
  original <- d
  p <- plot_matrix(d, metric = ~metric, frame = ~model)
  expect_s3_class(p, "plotly")
  expect_identical(d, original)
  spec <- matrix_spec(p)
  expect_identical(spec$models, c("tinyAM", "Accepted"))
  for (m in spec$metrics) {
    expect_identical(m$years, sort(unique(d$year)))
    expect_identical(m$ages, sort(unique(d$age)))
    expect_length(m$cells, 2L)
    for (i in seq_along(spec$models)) {
      z <- do.call(rbind, m$cells[[i]]$linear)
      logged <- do.call(rbind, m$cells[[i]]$log)
      expect_identical(dim(z), c(3L, 3L))
      for (a in seq_along(m$ages)) for (y in seq_along(m$years)) {
        row <- d[d$year == m$years[y] & d$age == m$ages[a] &
                   d$metric == m$name & d$model == spec$models[i], ]
        expected <- if (nrow(row)) row$est else NA_real_
        expect_equal(z[a, y], expected)
        expect_equal(logged[a, y], if (is.na(expected) || expected <= 0) NA_real_ else log10(expected))
      }
    }
    values <- d$est[d$metric == m$name]
    expect_lte(m$linear[1], min(values, na.rm = TRUE))
    expect_gte(m$linear[2], max(values, na.rm = TRUE))
    expect_true(any(grepl("projection", unlist(lapply(m$cells, `[[`, "text")))))
    expect_true(any(grepl("Source definition", unlist(lapply(m$cells, `[[`, "text")))))
  }
})

test_that("matrix plots handle constants, gaps, single cells and validation", {
  d <- data.frame(year = 2000, age = 1, est = 0)
  for (type in c("heatmap", "surface", "contour")) {
    p <- plot_matrix(d, type = type)
    s <- matrix_spec(p)
    expect_identical(s$type, type)
    expect_true(all(is.finite(s$metrics[[1]]$linear)))
    expect_equal(s$metrics[[1]]$log, c(0, 1))
    expect_equal(s$metrics[[1]]$cells[[1]]$linear, list(0))
    expect_true(is.na(s$metrics[[1]]$cells[[1]]$log[[1]]))
  }
  d$est <- NA_real_
  expect_equal(matrix_spec(plot_matrix(d))$metrics[[1]]$linear, c(0, 1))
  expect_error(plot_matrix(rbind(d, d)), "Duplicate")
  expect_error(plot_matrix(d, metric = ~abs(age)), "existing column")
  expect_error(plot_matrix(d, frame = ~missing), "existing column")
  expect_error(plot_matrix(d, frame = ~year + age), "existing column")
  d$year <- NA_real_
  expect_error(plot_matrix(d), "finite")
  d <- matrix_data()
  d$unit[2] <- "tonnes"
  expect_error(plot_matrix(d, metric = ~metric, frame = ~model), "different units")
  d <- matrix_data()
  d$model <- factor(d$model, levels = c("Accepted", "tinyAM"))
  expect_identical(matrix_spec(plot_matrix(d, metric = ~metric, frame = ~model))$models,
                   c("Accepted", "tinyAM"))
})

test_that("pattern tables align source aliases and independent matrix units", {
  tab <- data.frame(year = 2000:2001, age = 1, est = c(2, 3), unit = "thousand fish", se = .2)
  pop <- list(N = tab, biomass_at_age = transform(tab, unit = "tonnes"),
              ssb_mat = transform(tab, unit = "kt"), F = transform(tab, unit = "per year"))
  source <- pop
  d <- tinyAM:::.tidy_patterns(pop, reference = TRUE)
  expect_identical(pop, source)
  expect_equal(d$population$est, c(2000, 3000, 2000, 3000, 2e6, 3e6))
  expect_equal(d$population$unit, c("fish", "fish", rep("kg", 4)))
  expect_equal(d$mortality$est, c(2, 3))
  expect_false("se" %in% names(d$population))
  expect_null(tinyAM:::.tidy_patterns(list(N = tab), reference = TRUE)$mortality)
  pop$ssb_mat$unit <- "unknown"
  d <- tinyAM:::.tidy_patterns(pop, reference = TRUE)
  expect_true(all(is.na(d$population$est[d$population$metric == "SSB"])))
  pop <- list(N = tab, biomass_mat = tab, ssb_mat = tab)
  d <- tinyAM:::.tidy_patterns(pop, ssb_definition = "Age/year-specific spawning times")
  expect_equal(d$population$est, rep(c(2, 3), 3))
  expect_equal(d$population$definition[d$population$metric == "SSB"],
               rep("Age/year-specific spawning times", 2))
})

test_that("tidy patterns retain model labels separately from aggregate scales", {
  fit <- default_fit
  reference <- fit
  class(reference) <- c("tam_ref", "list")
  reference$comparison_scales <- c(ssb = 1e-3, N = 1e-3)
  reference$pop$N$unit <- "thousand fish"
  reference$pop$N$est <- reference$pop$N$est / 1000
  reference$pop$biomass_at_age <- reference$pop$biomass_mat
  reference$pop$biomass_at_age$unit <- "tonnes"
  reference$pop$biomass_at_age$est <- reference$pop$biomass_at_age$est / 1000
  reference$pop$ssb_mat$unit <- "kg"
  original <- list(Accepted = reference, tinyAM = fit)
  tabs <- tidy_tam(model_list = original, label_type = "factor")
  expect_setequal(as.character(tabs$patterns$population$model), names(original))
  d <- tabs$patterns$population
  for (metric in c("Abundance", "Biomass", "SSB")) {
    expect_equal(d$est[d$model == "tinyAM" & d$metric == metric],
                 d$est[d$model == "Accepted" & d$metric == metric])
  }
  expect_identical(fit, original$tinyAM)
  expect_identical(reference, original$Accepted)
})

test_that("matrix widget controls coordinate every selection independently", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is unavailable for JavaScript coordinator checks")
  result <- system2(node, c(shQuote(test_path("matrix-controls.js")),
    shQuote(system.file("js", "matrix.js", package = "tinyAM"))), stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "All matrix control combinations passed")
})
