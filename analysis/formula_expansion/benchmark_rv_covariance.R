# Compare covariance models using only saved standard-effort RV monitoring draws.
# From the package root:
# Rscript analysis/formula_expansion/benchmark_rv_covariance.R [checkpoint-directory]

study <- file.path("analysis", "formula_expansion")
default_result <- file.path(study, "results", "observation_error_z_covariance", "study-grid40")
args <- commandArgs(trailingOnly = TRUE)
result <- if (length(args)) args[[1L]] else default_result

rv_center_log_errors <- function(train, test) {
  cells <- split(train$log_error, train$cell)
  moments <- function(y) {
    y <- y[is.finite(y)]
    n <- length(y)
    if (n < 5L) return(c(n = n, mean = NA_real_, sd = NA_real_))
    c(n = n, mean = mean(y), sd = stats::sd(y))
  }
  m <- t(vapply(cells, moments, numeric(3)))
  j <- match(train$cell, rownames(m))
  i <- match(test$cell, rownames(m))
  # Account for loss of variance from estimating a calibration cell mean.
  n <- m[j, "n"]
  train$centered <- (train$log_error - m[j, "mean"]) *
    sqrt(ifelse(n >= 2L, n / (n - 1), NA_real_))
  test$centered <- test$log_error - m[i, "mean"]
  train$cell_sd <- unname(m[j, "sd"])
  test$cell_sd <- unname(m[i, "sd"])
  list(train = train, test = test)
}

# Benchmark independent blocks by survey draw. Include history in the block ID
# because replicate numbers and years are reused across histories.
rv_blocks <- function(data) {
  data$draw <- data$replicate
  data$replicate <- interaction(data$history_id, data$draw, drop = TRUE)
  d <- data[is.finite(data$centered) & is.finite(data$cell_sd) & data$cell_sd > 0, ]
  d <- d[order(d$replicate, d$year, d$age), ]
  d$z <- d$centered / d$cell_sd
  block <- interaction(d$replicate, d$year, drop = TRUE)
  following <- which(block[-1L] == head(block, -1L)) + 1L
  moments <- rowsum(cbind(n = 1, sum_z = d$z, sum_z2 = d$z^2, log_sd = log(d$cell_sd)),
    block, reorder = FALSE)
  keys <- c("scenario", "spatial", "mortality", "phi_year", "history_id", "year", "draw")
  list(d = d, block = block, following = following, moments = moments,
    gap = d$age[following] - d$age[following - 1L],
    key = d[!duplicated(block), keys, drop = FALSE])
}

rv_exchangeable_loss <- function(moments, rho) {
  n <- moments[, "n"]
  z_sum_sq <- moments[, "sum_z"]^2 / n
  centered_ss <- pmax(0, moments[, "sum_z2"] - z_sum_sq)
  determinant <- (n - 1) * log1p(-rho) + log1p((n - 1) * rho)
  quadratic <- centered_ss / (1 - rho) + z_sum_sq / (1 + (n - 1) * rho)
  unname((n * log(2 * pi) + determinant + quadratic) / 2 + moments[, "log_sd"])
}

rv_ar1_loss <- function(x, rho) {
  d <- x$d
  loss <- -stats::dnorm(d$z, log = TRUE) + log(d$cell_sd)
  current <- x$following
  phi <- rho^x$gap
  loss[current] <- -stats::dnorm(d$z[current], phi * d$z[current - 1L],
    sqrt(1 - phi^2), log = TRUE) + log(d$cell_sd[current])
  drop(rowsum(loss, x$block, reorder = FALSE))
}

rv_best_rho <- function(loss) {
  grid <- seq(0, .99, length.out = 41L)
  value <- vapply(grid, loss, numeric(1))
  i <- which.min(value)
  interval <- grid[pmax(1L, pmin(length(grid), i + c(-1L, 1L)))]
  optimum <- stats::optimize(loss, interval, tol = 1e-8)
  if (optimum$objective < value[i]) optimum$minimum else grid[i]
}

rv_benchmark_check <- function() {
  x <- data.frame(scenario = "test", spatial = "test", mortality = "test",
    phi_year = 0, history_id = rep(c("a", "b", "c"), c(3L, 3L, 1L)), year = 5,
    replicate = 1, age = c(1, 3, 4, 1, 3, 4, 9), centered = c(.4, -.2, .1, .1, .6, -.3, .1),
    cell_sd = c(.2, .5, .7, .4, .2, .4, .6))
  blocks <- rv_blocks(x)
  dense_loss <- function(rho, model) {
    vapply(split(x, x$history_id), function(d) {
      R <- if (model == "ar1") rho^abs(outer(d$age, d$age, "-")) else
        diag(1 - rho, nrow(d)) + rho
      L <- chol(R * outer(d$cell_sd, d$cell_sd))
      (nrow(d) * log(2 * pi) + sum(backsolve(L, d$centered, transpose = TRUE)^2)) / 2 +
        sum(log(diag(L)))
    }, numeric(1))
  }
  for (rho in c(0, .2, .8)) {
    stopifnot(isTRUE(all.equal(unname(rv_ar1_loss(blocks, rho)),
      unname(dense_loss(rho, "ar1")), tolerance = 1e-12)),
      isTRUE(all.equal(rv_exchangeable_loss(blocks$moments, rho),
        unname(dense_loss(rho, "exchangeable")), tolerance = 1e-12)))
  }
}

rv_benchmark_fold <- function(data, heldout) {
  split <- ifelse(data$replicate <= 10L, 1L, 2L)
  centered <- rv_center_log_errors(data[split != heldout, , drop = FALSE],
    data[split == heldout, , drop = FALSE])
  train <- rv_blocks(centered$train)
  test <- rv_blocks(centered$test)
  stopifnot(length(train$following) > 0L, nrow(test$key) > 0L)
  age_rho <- rv_best_rho(function(rho) sum(rv_ar1_loss(train, rho)))
  shared_rho <- rv_best_rho(function(rho) sum(rv_exchangeable_loss(train$moments, rho)))
  n <- test$moments[, "n"]
  iid <- rv_ar1_loss(test, 0) / n
  stopifnot(isTRUE(all.equal(unname(iid), rv_exchangeable_loss(test$moments, 0) / unname(n),
    tolerance = 1e-12)))
  ar1 <- if (age_rho == 0) iid else rv_ar1_loss(test, age_rho) / n
  exchangeable <- if (shared_rho == 0) iid else rv_exchangeable_loss(test$moments, shared_rho) / n
  data.frame(test$key, fold = heldout, n_age = unname(n),
    iid = unname(iid), ar1 = unname(ar1), exchangeable = unname(exchangeable),
    ar1_rho = age_rho, exchangeable_rho = shared_rho,
    training_age_pairs = length(train$following))
}

rv_benchmark_summary <- function(history, keys) {
  metrics <- c("ar1_minus_iid", "exchangeable_minus_iid", "exchangeable_minus_ar1")
  groups <- split(history, history[keys], drop = TRUE)
  do.call(rbind, lapply(groups, function(d) {
    stopifnot(!anyDuplicated(d$history_id))
    n <- nrow(d)
    do.call(rbind, lapply(metrics, function(metric) {
      mean <- mean(d[[metric]])
      half <- if (n > 1L) stats::qt(.975, n - 1) * stats::sd(d[[metric]]) / sqrt(n) else NA_real_
      data.frame(d[1L, keys, drop = FALSE], comparison = metric,
        estimate = mean, lower = mean - half, upper = mean + half, n_histories = n)
    }))
  }))
}

rv_benchmark_check()
manifest_file <- file.path(result, "manifest.rds")
manifest <- readRDS(manifest_file)
stopifnot(identical(manifest$settings$phase, "study"),
  identical(manifest$settings$repeats, 20L), identical(manifest$settings$effort, "standard"))
case_files <- as.vector(outer(manifest$scenarios$scenario, manifest$histories$history_id,
  function(scenario, history) file.path(result, "checkpoints", scenario, history, "case.rds")))
cases <- lapply(case_files, function(path) {
  d <- as.data.frame(readRDS(path)$rows)
  d[d$component == "RV" & d$bin == "modeled", ]
})
rv <- do.call(rbind, cases)
stopifnot(nrow(rv) > 0L, all(rv$effort == "standard"),
  !anyDuplicated(rv[c("scenario", "history_id", "year", "age", "replicate")]))
rv$cell <- interaction(rv$history_id, rv$year, rv$age, drop = TRUE)
positive <- is.finite(rv$obs) & rv$obs > 0 & is.finite(rv$truth) & rv$truth > 0
rv$log_error <- NA_real_
rv$log_error[positive] <- log(rv$obs[positive] / rv$truth[positive])

draws <- do.call(rbind, lapply(split(rv, rv$scenario), function(data) {
  rbind(rv_benchmark_fold(data, 1L), rv_benchmark_fold(data, 2L))
}))
draws$ar1_minus_iid <- draws$ar1 - draws$iid
draws$exchangeable_minus_iid <- draws$exchangeable - draws$iid
draws$exchangeable_minus_ar1 <- draws$exchangeable - draws$ar1
metrics <- c("iid", "ar1", "exchangeable", "ar1_minus_iid",
  "exchangeable_minus_iid", "exchangeable_minus_ar1")
keys <- c("scenario", "spatial", "mortality", "phi_year")
snapshots <- stats::aggregate(draws[metrics], draws[c(keys, "history_id", "year")], mean)
history <- stats::aggregate(snapshots[metrics], snapshots[c(keys, "history_id")], mean)
by_regime <- rv_benchmark_summary(history, keys)
core <- history[history$mortality != "temporal_only", ]
spatial_history <- stats::aggregate(core[metrics], core[c("spatial", "history_id")], mean)
by_spatial <- rv_benchmark_summary(spatial_history, "spatial")
correlations <- unique(draws[c(keys, "fold", "ar1_rho", "exchangeable_rho", "training_age_pairs")])
support <- data.frame(n_planned_draws = nrow(manifest$scenarios) * nrow(manifest$histories) *
    length(manifest$settings$years) * manifest$settings$repeats,
  n_scored_draws = nrow(draws), n_input_age_rows = nrow(rv),
  n_positive_age_rows = sum(positive), n_scored_age_rows = sum(draws$n_age))
result_folder <- file.path(default_result, "rv_covariance_benchmark")
dir.create(result_folder, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(by_regime = by_regime, by_spatial = by_spatial, history = history,
  draws = draws, correlations = correlations, support = support,
  provenance = list(manifest = tools::md5sum(manifest_file), checkpoints = tools::md5sum(case_files),
    script = tools::md5sum(file.path(study, "benchmark_rv_covariance.R")),
    R = as.character(getRversion()))),
  file.path(result_folder, "benchmark.rds"))
for (name in c("by_regime", "by_spatial", "history", "correlations", "support")) {
  utils::write.csv(get(name), file.path(result_folder, paste0(name, ".csv")), row.names = FALSE)
}
print(by_spatial, row.names = FALSE, digits = 3)
print(support, row.names = FALSE)
