# Use separate monitoring draws to diagnose observation means on three stocks.
.oe_pilot_center <- function(calibration, train) {
  d <- subset(calibration, is.finite(truth) & truth > 0 & is.finite(obs) & obs >= 0)
  d$ratio <- d$obs / d$truth
  positive <- subset(d, ratio > 0)
  positive$log_ratio <- log(positive$ratio)
  cells <- stats::aggregate(log_ratio ~ year + age, positive,
    function(x) if (length(x) >= 3L) mean(x) else NA_real_)
  ages <- stats::aggregate(log_ratio ~ age, positive, mean)
  shift <- cells$log_ratio[match(.oe_fit_key(train), .oe_fit_key(cells))]
  fallback <- !is.finite(shift)
  shift[fallback] <- ages$log_ratio[match(train$age[fallback], ages$age)]
  shift[!is.finite(shift)] <- mean(positive$log_ratio)
  stopifnot(all(is.finite(shift)))
  adjusted <- train
  adjusted$obs <- train$obs * exp(-shift)
  adjusted$se <- train$se * exp(-shift)
  stopifnot(identical(is.na(train$obs), is.na(adjusted$obs)),
    identical(train$obs == 0, adjusted$obs == 0))
  a <- stats::aggregate(ratio ~ age, d, mean)
  a$geometric_ratio <- exp(ages$log_ratio[match(a$age, ages$age)])
  a$n_fallback <- as.integer(table(factor(train$age[fallback], levels = a$age)))
  list(train = adjusted, ratios = a)
}

.oe_pilot_diagnostic_stock <- function(stock, cache, index) {
  stock_cache <- file.path(out, stock$population_id)
  input <- required(readRDS(file.path(stock_cache, "inputs.rds")))
  original <- readRDS(file.path(stock_cache, "result.rds"))
  pop <- oe_pilot_population(stock$population_seed, stock$field_seed, grid)
  n <- 20L
  calibration <- vector("list", length(pop$years) * n)
  for (i in seq_along(pop$years)) {
    admit()
    year <- pop$years[i]
    snapshot <- oe_snapshot(pop, year)
    removals <- required(readRDS(file.path(stock_cache, sprintf("removals-%02d.rds", year))))
    for (replicate in seq_len(n)) {
      base_seed <- 910000000L + 10000L * index + 50L * (year - 1L) + replicate
      draw <- .oe_pilot_draw(snapshot, removals,
        c(Gillnet = base_seed, RV = base_seed + 2000L), cache, sprintf("cal-%02d", replicate))
      calibration[(i - 1L) * n + replicate] <- list(draw$value)
    }
    rm(snapshot, removals)
    invisible(gc())
  }
  train <- input$histories$primary$fit
  centered <- lapply(names(train), function(nm) {
    x <- .oe_pilot_center(.oe_pilot_stack(calibration, nm), train[[nm]])
    x$ratios$component <- if (nm == "catch") "Gillnet" else "RV"
    x
  })
  names(centered) <- names(train)
  fit <- checkpoint(file.path(cache, "baseline-centered.rds"), function() {
    .oe_pilot_fit(input$base, lapply(centered, "[[", "train"), NULL, input$truth, "constant")
  })
  before <- if (NROW(original$annual)) {
    subset(original$annual, model == "constant" & draw == "primary" & quantity == "ssb")
  }
  after <- if (isTRUE(fit$value$usable)) subset(fit$value$annual, quantity == "ssb")
  metrics <- function(d, data) {
    valid <- NROW(d) == length(pop$years) && all(is.finite(d$relative_error))
    data.frame(data = data, usable = valid,
      abs_error = if (valid) 100 * mean(abs(d$relative_error)) else NA_real_,
      bias = if (valid) 100 * mean(d$relative_error) else NA_real_,
      coverage = if (valid) 100 * mean(d$covered) else NA_real_)
  }
  list(ratios = do.call(rbind, lapply(centered, "[[", "ratios")),
    recovery = rbind(metrics(before, "Primary"), metrics(after, "Known log centers")),
    successful_draws = sum(!vapply(calibration, is.null, logical(1))), attempted_draws = length(calibration))
}

oe_pilot_mean_diagnostic <- function() {
  for (file in c("observation_error_pilot_populations.R", "observation_error_fit.R")) {
    source(file.path(study, file), encoding = "UTF-8")
  }
  manifest_path <- file.path(out, "population-manifest.rds")
  manifest <- readRDS(manifest_path)$seed_design$populations[1:3, ]
  root <- file.path(out, "mean-diagnostic")
  dir.create(root, showWarnings = FALSE)
  verify_manifest(list(source = tools::md5sum(file.path(study, "observation_error_pilot_means.R")),
    pilot = tools::md5sum(manifest_path)), file.path(root, "manifest.rds"))
  results <- lapply(seq_len(nrow(manifest)), function(i) {
    stock <- manifest[i, , drop = FALSE]
    cache <- file.path(root, stock$population_id)
    dir.create(cache, showWarnings = FALSE)
    ans <- checkpoint(file.path(cache, "result.rds"), function() {
      .oe_pilot_diagnostic_stock(stock, cache, i)
    })
    if (!is.null(ans$value)) {
      ans$value$ratios$regime <- ans$value$recovery$regime <- stock$regime
      ans$value$ratios$population_id <- ans$value$recovery$population_id <- stock$population_id
    }
    ans
  })
  values <- lapply(results, "[[", "value")
  save_checkpoint(list(ratios = .oe_pilot_stack(values, "ratios"),
    recovery = .oe_pilot_stack(values, "recovery"),
    results = results, pilot = unname(tools::md5sum(manifest_path))),
    file.path(out, "mean-diagnostic.rds"))
}
