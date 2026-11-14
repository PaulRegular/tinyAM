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

pld_get_family <- function(use_pld, fixed, kind) {
  family <- if (use_pld) "pld_lognormal" else "lognormal"
  constructor <- getExportedValue("tinyAM", family)
  if (!use_pld) return(constructor())
  constructor(y50 = if (fixed) unname(pld_y50[[kind]]) else NULL)
}

pld_metrics <- function(fit, annual, reference, plus_age) {
  population <- list(ssb = tinyAM::get_ssb(fit),
    recruitment = tinyAM::get_recruitment(fit),
    abundance = fit$pop$abundance,
    M_bar = tinyAM::get_mbar(fit),
    F_bar = tinyAM::get_fbar(fit))
  annual <- do.call(rbind, lapply(names(population), function(name) {
    p <- population[[name]][c("year", "est", "lwr", "upr")]
    p$quantity <- name
    p$truth <- if (name == "F_bar") {
      reference$F_bar[match(p$year, reference$year)]
    } else annual[[name]][match(p$year, annual$year)]
    p
  }))
  f <- tinyAM::get_f(fit)
  f <- f[f$age < plus_age, c("year", "age", "est", "lwr", "upr")]
  f$truth <- pld_f(f$age)
  pulse <- annual[annual$quantity == "M_bar", ]
  pulse_height <- function(x) mean(x[pulse$year %in% 9:12]) -
    mean(x[pulse$year %in% c(1:4, 17:20)])
  pulse_center <- function(x) {
    base <- mean(x[pulse$year %in% c(1:4, 17:20)])
    weight <- pmax(x - base, 0)
    if (sum(weight) > 0) sum(pulse$year * weight) / sum(weight) else NA_real_
  }
  list(annual = annual, F = f,
    pulse = data.frame(height = pulse_height(pulse$est), height_truth = pulse_height(pulse$truth),
      center = pulse_center(pulse$est), center_truth = pulse_center(pulse$truth)))
}

pld_fit_model <- function(data, annual, reference, model, plus_age) {
  fitter <- if (model$ablated) pld_occurrence_ablated_fit() else tinyAM::fit_tam
  fit <- fitter(data, ages = seq_len(plus_age), silent = TRUE,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "rw"),
    F_settings = list(process = "off", mu_form = ~ normal(age), mean_ages = 5:10),
    M_settings = list(process = "ar1", mu_form = ~ 1, mu_supplied = NULL,
      first_dev_year = 1L, age_breaks = c(1L, plus_age + 1L), mean_ages = 5:10),
    catch_settings = list(sd_form = ~ 1, family = pld_get_family(model$catch, model$fixed, "catch")),
    index_settings = list(q_form = ~ logistic(age), sd_form = ~ 1,
      family = pld_get_family(model$index, model$fixed, "index")))
  columns <- c("year", "age", "obs", "pred", "sd", "family", "p_zero", "log_pred_positive")
  obs <- lapply(fit$obs_pred[c("catch", "index")], function(d) {
    d[, intersect(columns, names(d)), drop = FALSE]
  })
  half <- fit$fixed_par[fit$fixed_par$par %in% c("y50_catch", "y50_index"), , drop = FALSE]
  list(converged = fit$is_converged, numerical = fit$diagnostics$numerical,
    objective = fit$opt$objective, optimizer_code = fit$opt$convergence,
    gradient = fit$diagnostics$max_gradient,
    metric = if (fit$is_converged) pld_metrics(fit, annual, reference, plus_age),
    y50 = half, obs = obs)
}

pld_reference_fbar <- function(rows, values = pld_f(rows$age)) {
  data <- transform(rows, FN = values * N)
  use <- data$age %in% 5:10
  d <- stats::aggregate(data[use, c("N", "FN")], list(year = data$year[use]), sum)
  data.frame(year = d$year, F_bar = d$FN / d$N)
}

pld_fit_input <- function(input, out, checkpoint) {
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
    reference <- draw$prediction
    reference$N <- reference$index / pld_q(reference$age)
    fbar <- pld_reference_fbar(reference)
    outcomes <- draw[settings$scenario]
  } else {
    settings <- expand.grid(scenario = names(draw$series), treatment = 1:4, age = c(14L, 20L),
      stringsAsFactors = FALSE)
    settings <- settings[settings$age == 14L |
      (settings$scenario == "standard" & draw$id <= 8L), ]
    settings$name <- paste0(settings$scenario, "-a", settings$age)
    outcomes <- lapply(seq_len(nrow(settings)), function(i) {
      draw$series[[settings$scenario[i]]]$data[[match(settings$age[i], c(14L, 20L))]]
    })
    fbar <- pld_reference_fbar(draw$reference)
  }
  for (i in seq_len(nrow(settings))) {
    row <- settings[i, ]
    fit_id <- paste(prefix, row$name, sprintf("model%02d", row$treatment), sep = "-")
    file <- file.path(out, paste0("fit-", fit_id, ".rds"))
    checkpoint(file, function() pld_fit_model(outcomes[[i]], draw$annual,
      fbar, pld_treatments[row$treatment, ], row$age))
  }
  invisible(NULL)
}
