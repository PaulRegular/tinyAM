pld_treatments <- data.frame(
  model = c("Lognormal", "Survey PLD", "Catch PLD", "Both PLD", "Fixed detection", "No occurrence"),
  catch = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  index = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE),
  fixed = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  ablated = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))

# Keep the production objective and fitter. Change one summand in private
# copies of the functions so the ablation uses identical parameter maps,
# optimizer, starts, positive density and conditional lognormal SDs.
pld_occurrence_ablated_fit <- function() {
  fit <- tinyAM::fit_tam
  nll <- tinyAM::nll_fun
  occurrence <- list(quote(pld$log_p_zero[dat$pld$zero]), quote(pld$log_p_positive[positive]))
  replaced <- integer(length(occurrence))
  drop_occurrence <- function(code) {
    index <- which(vapply(occurrence, identical, logical(1), y = code))
    if (length(index)) {
      replaced[index] <<- replaced[index] + 1L
      return(0)
    }
    if (!is.call(code)) return(code)
    entries <- as.list(code)
    for (i in seq_along(entries)[-1L]) {
      if (is.call(entries[[i]])) {
        entries[i] <- list(drop_occurrence(entries[[i]]))
      }
    }
    as.call(entries)
  }
  body(nll) <- drop_occurrence(body(nll))
  if (any(replaced != 1L)) stop("The observation likelihood has changed. Audit the ablation.")
  private <- new.env(parent = environment(fit))
  private$nll_fun <- nll
  environment(fit) <- private
  fit
}

pld_get_family <- function(use_pld, fixed, detection, kind) {
  family <- if (use_pld) "pld_lognormal" else "lognormal"
  constructor <- getExportedValue("tinyAM", family)
  if (!use_pld) return(constructor())
  constructor(y50 = if (fixed) unname(detection[[kind]]) else NULL)
}

pld_weighted_mortality <- function(rows, values) {
  use <- rows$age %in% 5:10 & is.finite(values) & rows$N > 0
  stopifnot(any(use))
  d <- stats::aggregate(list(N = rows$N[use], weighted = (rows$N * values)[use]),
    list(year = rows$year[use]), sum)
  data.frame(year = d$year, est = d$weighted / d$N)
}

pld_mortality_events <- function(annual) {
  event <- function(quantity, component, truth) {
    d <- annual[annual$quantity == quantity, ]
    x <- if (truth) d$truth else d$est
    baseline <- mean(x[d$year %in% pld_windows[[paste0(component, "_baseline")]]])
    height <- mean(x[d$year %in% pld_windows[[paste0(component, "_peak")]]]) - baseline
    event <- d$year %in% pld_windows[[paste0(component, "_event")]]
    weight <- pmax(x[event] - baseline, 0)
    center <- if (sum(weight) > 0) sum(d$year[event] * weight) / sum(weight) else NA_real_
    c(height = height, center = center, late = mean(x[d$year %in% pld_windows$late]))
  }
  natural <- event("M_bar", "M", FALSE)
  natural_truth <- event("M_bar", "M", TRUE)
  fishing <- event("F_bar", "F", FALSE)
  fishing_truth <- event("F_bar", "F", TRUE)
  data.frame(height = natural[["height"]], height_truth = natural_truth[["height"]],
    center = natural[["center"]], center_truth = natural_truth[["center"]],
    M_late = natural[["late"]], M_late_truth = natural_truth[["late"]],
    F_height = fishing[["height"]], F_height_truth = fishing_truth[["height"]],
    F_center = fishing[["center"]], F_center_truth = fishing_truth[["center"]],
    F_late = fishing[["late"]], F_late_truth = fishing_truth[["late"]])
}

pld_metrics <- function(fit, annual, reference, plus_age) {
  fit_mortality <- function(x, values, weighted_interval) {
    reference_rows <- merge(reference, x[c("year", "age", "est")],
      by = c("year", "age"))
    score <- pld_weighted_mortality(reference_rows, reference_rows$est)
    truth <- pld_weighted_mortality(reference, values)
    score$truth <- truth$est[match(score$year, truth$year)]
    # Fbar intervals from tinyAM use fitted N weights. They do not cover
    # this true-N-weighted Fbar. The model shares one M across all ages.
    interval <- weighted_interval[match(score$year, weighted_interval$year), ]
    cbind(score, interval[c("lwr", "upr")])
  }
  f <- tinyAM::get_f(fit)
  m <- tinyAM::get_m(fit)
  unavailable <- data.frame(year = annual$year, lwr = NA_real_, upr = NA_real_)
  fbar <- fit_mortality(f, reference$F_realized, unavailable)
  mbar <- fit_mortality(m, reference$M_realized, tinyAM::get_mbar(fit))
  population <- list(ssb = tinyAM::get_ssb(fit),
    recruitment = tinyAM::get_recruitment(fit),
    abundance = fit$pop$abundance,
    M_bar = mbar,
    F_bar = fbar)
  annual <- do.call(rbind, lapply(names(population), function(name) {
    p <- population[[name]][c("year", "est", "lwr", "upr")]
    p$quantity <- name
    p$truth <- if (name %in% c("F_bar", "M_bar")) population[[name]]$truth else
      annual[[name]][match(p$year, annual$year)]
    p
  }))
  f <- f[f$age < plus_age, c("year", "age", "est", "lwr", "upr")]
  f <- merge(f, data.frame(year = reference$year, age = reference$age,
    truth = reference$F_realized), by = c("year", "age"))
  truth_age <- stats::aggregate(list(truth = reference$N / 1000),
    list(year = reference$year, age = pmin(reference$age, plus_age)), sum)
  n <- fit$pop$N[c("year", "age", "est", "unit")]
  stopifnot(all(n$unit == "thousand fish"), nrow(n) == nrow(truth_age),
    !anyDuplicated(n[c("year", "age")]))
  n <- merge(n, truth_age, by = c("year", "age"))
  n$lwr <- n$upr <- NA_real_
  n_annual <- stats::aggregate(est ~ year, n, sum)
  abundance <- annual[annual$quantity == "abundance", ]
  stopifnot(nrow(n) == nrow(truth_age), isTRUE(all.equal(n_annual$est,
    abundance$est[match(n_annual$year, abundance$year)], tolerance = 1e-10)))
  list(annual = annual, F = f,
    population_age = n, pulse = pld_mortality_events(annual))
}

pld_fit_model <- function(data, annual, reference, model, plus_age, detection) {
  fitter <- if (model$ablated) pld_occurrence_ablated_fit() else tinyAM::fit_tam
  fit <- fitter(data, ages = seq_len(plus_age), silent = TRUE,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "rw"),
    F_settings = list(process = "off", mu_form = ~ factor(year) + normal(age), mean_ages = 5:10),
    M_settings = list(process = "ar1", mu_form = ~ 1, mu_supplied = NULL,
      first_dev_year = 1L, age_breaks = c(1L, plus_age + 1L), mean_ages = 5:10),
    catch_settings = list(sd_form = ~ 1,
      family = pld_get_family(model$catch, model$fixed, detection, "catch")),
    index_settings = list(q_form = ~ logistic(age), sd_form = ~ 1,
      family = pld_get_family(model$index, model$fixed, detection, "index")))
  columns <- c("year", "age", "obs", "pred", "sd", "family", "p_zero", "log_pred_positive")
  obs <- lapply(fit$obs_pred[c("catch", "index")], function(d) {
    d[, intersect(columns, names(d)), drop = FALSE]
  })
  half <- fit$fixed_par[fit$fixed_par$par %in% c("y50_catch", "y50_index"), , drop = FALSE]
  list(converged = fit$is_converged, numerical = fit$diagnostics$numerical,
    objective = fit$opt$objective, optimizer_code = fit$opt$convergence,
    gradient = fit$diagnostics$max_gradient,
    metric = if (fit$is_converged) pld_metrics(fit, annual, reference, plus_age),
    y50 = half, calibration = if (model$fixed) detection, obs = obs)
}

pld_calibration_rows <- function(draw, kind) {
  series <- draw$series[["standard"]]
  raw <- series$observations[[paste0("native_", kind)]][c("year", "age", "obs")]
  truth <- if (kind == "catch") {
    ref <- draw$reference
    data.frame(year = ref$year, age = ref$age, obs = ref$catch_realized / 1000)
  } else {
    transform(series$expected_index, obs = mu)[c("year", "age", "obs")]
  }
  pooled <- merge(pld_pool(raw, 14L), pld_pool(truth, 14L), by = c("year", "age"),
    suffixes = c("", "_truth"))
  stopifnot(nrow(pooled) == 14L * length(pld_years), !any(pooled$obs > 0 & pooled$obs_truth == 0, na.rm = TRUE))
  data.frame(history = draw$id, pooled)
}

pld_calibrate <- function(files) {
  completed <- lapply(files, readRDS)
  stopifnot(length(completed) == 6L, all(vapply(completed, `[[`, logical(1), "ok")))
  histories <- lapply(completed, `[[`, "value")
  for (id in seq_along(histories)) {
    stopifnot(histories[[id]]$id == id, identical(histories[[id]]$seeds, pld_seeds(id + 100L)))
  }
  estimate <- function(kind) {
    rows <- do.call(rbind, lapply(histories, pld_calibration_rows, kind = kind))
    fit <- function(d) {
      d <- d[is.finite(d$obs) & is.finite(d$obs_truth) & d$obs_truth > 0, ]
      zero <- d$obs == 0
      stopifnot(any(zero), any(!zero))
      range <- range(log(d$obs_truth)) + c(-20, 20)
      result <- stats::optimize(function(log_y50) {
        lambda <- exp(log(log(2)) + log(d$obs_truth) - log_y50)
        sum(lambda[zero]) - sum(log(-expm1(-lambda[!zero])))
      }, range)
      c(y50 = exp(result$minimum), zeros = sum(zero), observed = nrow(d))
    }
    data.frame(component = kind, as.list(fit(rows)),
      omitted_at_zero_mean = sum(is.finite(rows$obs) & rows$obs_truth == 0),
      missing = sum(is.na(rows$obs)))
  }
  do.call(rbind, lapply(c("catch", "index"), estimate))
}

pld_fit_input <- function(input, out, checkpoint, detection = pld_y50) {
  record <- readRDS(input)
  if (!record$ok) return(invisible(NULL))
  draw <- record$value
  prefix <- sub("[.]rds$", "", basename(input))
  if (startsWith(prefix, "controlled")) {
    settings <- expand.grid(scenario = c("exact", "multinomial"), treatment = 1:6,
      stringsAsFactors = FALSE)
    settings <- settings[!(settings$scenario == "multinomial" & settings$treatment > 4L), ]
    settings$name <- paste0(settings$scenario, "-a14")
    settings$age <- 14L
    reference <- transform(draw$prediction, F_realized = pld_f(age, year), M_realized = pld_m(year))
    outcomes <- draw[settings$scenario]
  } else {
    settings <- expand.grid(scenario = names(draw$series), treatment = 1:6, age = 14L,
      stringsAsFactors = FALSE)
    settings <- settings[settings$scenario == "standard" | settings$treatment %in% c(1L, 4L), ]
    settings$name <- paste0(settings$scenario, "-a", settings$age)
    outcomes <- lapply(seq_len(nrow(settings)), function(i) {
      draw$series[[settings$scenario[i]]]$data[[1L]]
    })
    reference <- draw$reference
  }
  for (i in seq_len(nrow(settings))) {
    row <- settings[i, ]
    fit_id <- paste(prefix, row$name, sprintf("model%02d", row$treatment), sep = "-")
    file <- file.path(out, paste0("fit-", fit_id, ".rds"))
    checkpoint(file, function() pld_fit_model(outcomes[[i]], draw$annual,
      reference, pld_treatments[row$treatment, ], row$age, detection))
  }
  invisible(NULL)
}
