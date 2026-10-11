#' Check numerical convergence and model diagnostics
#'
#' @description
#' Reviews a fitted model's optimization, uncertainty, formula designs and
#' residual patterns. Numerical convergence, known structural problems and
#' advisory findings are reported separately. A converged model can still have
#' weakly estimated parameters or unsuitable biological assumptions.
#'
#' @details
#' Numerical checks require optimizer success, finite estimates and predictions,
#' a sufficiently small gradient, and successful uncertainty calculation with a
#' positive-definite Hessian. At an active parameter bound, the gradient is
#' adjusted for the direction in which the parameter is allowed to move. Both
#' raw and adjusted gradients are retained.
#'
#' Structural checks look for redundant formula coefficients. Passing these
#' checks does not prove that abundance, mortality and catchability are
#' distinguishable from the available data. Advisory findings include strong
#' parameter correlations, variance or correlation boundaries, sparse data and
#' conditional residual patterns. These do not change numerical convergence.
#' Residual summaries exclude projections and missing or zero observations.
#' Conditional residuals depend on fitted states and their shrinkage; they are
#' not independent tests of the assumed observation distribution.
#'
#' By default the function uses information already available in the fit and
#' does not refit the model. `detailed = TRUE` also examines the fixed-parameter
#' covariance eigenvalues and weakly supported parameter combinations. Profiles,
#' one-step residuals, simulations and retrospectives are separate investigations.
#' Mortality-mean advisories highlight jointly estimated M mean/residual
#' variation and F/M mean variation. Formula-effect advisories flag 95% AR1
#' correlation intervals wider than 0.5 and process-SD intervals spanning more
#' than a factor of ten, including catchability effects. Catchability groups
#' with fewer than five observed effect levels are flagged when process SD or
#' correlation is estimated. Observation-specific IID effects with estimated
#' observation SD receive a caution under the logit link; exact log-link
#' variance aliases are rejected before fitting.
#' Logistic curves are flagged when the fitted midpoint is outside the observed
#' age/size range and its 95% interval is wider than that range.
#' Stock-recruit advisories flag limited parent-SSB contrast, wide curve or AR1
#' intervals, and absolute correlation above 0.5 between adjacent recruitment
#' innovations when at least 15 historical innovations are available. The first
#' stationary AR1 residual is excluded. Innovations depend on fitted states;
#' this check can reveal an unsuitable residual process but does not establish
#' or rule out a biological stock-recruit relationship.
#' These descriptive thresholds identify imprecise estimates; they do not prove
#' a structural problem or change the convergence criteria.
#'
#' @param fit A fitted `tam_fit` object.
#' @param grad_tol Positive gradient tolerance. `NULL` uses the tolerance stored
#'   in the fit, or `1e-2` if none is stored.
#' @param detailed Examine covariance eigenvalues and weak parameter directions?
#' @return A `tam_check` list with `is_converged`, `numerical`, `structural`,
#'   `structural_status`, `advisories`, `residuals`, `data`, and `curvature`.
#'   `numerical` records pass, fail or not-assessed status for each check.
#'   Calling this function does not modify the fit or emit warnings; its print
#'   method displays the findings.
#' @examples
#' \dontrun{
#' checks <- check_tam(fit)
#' checks$is_converged
#' checks$residuals
#' check_tam(fit, detailed = TRUE)
#' }
#' @seealso [fit_tam()], [sim_tam()], [fit_retro()], [tidy_obs_pred()]
#' @export
check_tam <- function(fit, grad_tol = NULL, detailed = FALSE) {
  fit <- .require_tam_fit(fit, arg = "fit")
  if (is.null(grad_tol)) grad_tol <- if (is.null(fit$grad_tol)) 1e-2 else fit$grad_tol
  if (!is.numeric(grad_tol) || length(grad_tol) != 1L ||
      !is.finite(grad_tol) || grad_tol <= 0) {
    cli::cli_abort("{.arg grad_tol} must be one positive finite number.")
  }
  if (!is.logical(detailed) || length(detailed) != 1L || is.na(detailed)) {
    cli::cli_abort("{.arg detailed} must be TRUE or FALSE.")
  }

  opt <- if (is.list(fit$opt)) fit$opt else NULL
  sdr <- if (is.list(fit[["sdrep"]]) && !inherits(fit[["sdrep"]], "try-error")) fit[["sdrep"]] else NULL
  gradient <- if (!is.null(sdr$gradient.fixed)) sdr$gradient.fixed else fit$gradient
  parameters <- if (!is.null(opt$par)) opt$par else sdr$par.fixed
  bounds <- fit$bounds
  if (is.null(bounds) && !is.null(parameters)) {
    lower <- rep(-Inf, length(parameters))
    lower[names(parameters) == "dq"] <- 0
    bounds <- list(lower = lower, upper = rep(Inf, length(parameters)))
  }
  projected <- .project_gradient(gradient, parameters, bounds)
  max_gradient <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (any(!is.finite(x))) return(Inf)
    if (length(x)) max(abs(x)) else 0
  }
  raw_max <- max_gradient(gradient)
  adjusted_max <- max_gradient(projected$gradient)
  status <- function(value) {
    if (length(value) != 1L || is.na(value)) "not_assessed" else if (value) "pass" else "fail"
  }
  row <- function(check, value, detail) {
    data.frame(check = check, status = status(value), detail = detail,
               stringsAsFactors = FALSE)
  }
  finite_parameters <- if (is.null(parameters)) NA else all(is.finite(parameters))
  if (!is.null(sdr$par.random)) finite_parameters <- isTRUE(finite_parameters) &&
    all(is.finite(sdr$par.random))
  covariance <- sdr$cov.fixed
  uncertainty_ok <- if (is.null(sdr)) NA else {
    if (is.null(covariance)) NA else all(is.finite(covariance)) &&
      all(diag(covariance) >= 0)
  }
  numerical <- rbind(
    row("optimizer", if (is.null(opt$convergence)) NA else opt$convergence == 0,
        if (is.null(opt)) "Optimizer result unavailable." else
          paste0("Code ", opt$convergence, ": ", opt$message)),
    row("objective", if (is.null(opt$objective)) NA else is.finite(opt$objective),
        paste("Objective:", if (is.null(opt$objective)) "unavailable" else signif(opt$objective, 7))),
    row("parameters", finite_parameters, "Fixed estimates and available random states must be finite."),
    row("gradient", if (is.na(adjusted_max)) NA else
          is.finite(adjusted_max) && adjusted_max <= grad_tol,
        sprintf("Adjusted maximum %s; raw maximum %s; tolerance %s.",
                signif(adjusted_max, 3), signif(raw_max, 3), grad_tol)),
    row("uncertainty", uncertainty_ok, if (!is.null(fit$sdreport_error))
          fit$sdreport_error else "Fixed-parameter covariance must be available and finite."),
    row("hessian", if (is.null(sdr$pdHess)) NA else isTRUE(sdr$pdHess),
        "The Hessian must be positive definite."),
    row("predictions", if (is.null(fit$rep$log_pred) || is.null(fit$rep$sd_obs)) NA else
          all(is.finite(fit$rep$log_pred)) &&
          all(is.finite(fit$rep$sd_obs) & fit$rep$sd_obs > 0),
        "Log predictions and positive observation SDs must be finite.")
  )
  structural <- .check_tam_structure(fit$dat, fit$parameter_map)
  structural_status <- if (!nrow(structural) || all(structural$status == "not_assessed")) {
    "not_assessed"
  } else if (any(structural$status == "fail")) "issues" else "no_known_issues"
  residuals <- .check_tam_residuals(fit$obs_pred)
  data <- .check_tam_data_summary(fit$dat)
  advisories <- .check_tam_advisories(fit, projected$active, residuals, data)
  curvature <- if (detailed) .check_tam_curvature(covariance, names(parameters)) else NULL
  structure(list(
    is_converged = all(numerical$status == "pass"),
    numerical = numerical, max_gradient = adjusted_max, raw_max_gradient = raw_max,
    grad_tol = grad_tol, optimizer_code = opt$convergence,
    optimizer_message = opt$message, pd_hessian = sdr$pdHess,
    active_bounds = projected$active, structural = structural,
    structural_status = structural_status, advisories = advisories,
    residuals = residuals, data = data, curvature = curvature
  ), class = "tam_check")
}

.project_gradient <- function(gradient, parameters, bounds) {
  active <- rep(FALSE, length(gradient))
  if (!is.null(gradient) && !is.null(parameters) &&
      length(gradient) == length(parameters) && !is.null(bounds)) {
    at_lower <- is.finite(bounds$lower) & !is.na(parameters) & parameters == bounds$lower
    at_upper <- is.finite(bounds$upper) & !is.na(parameters) & parameters == bounds$upper
    active <- at_lower | at_upper
    lower <- at_lower & is.finite(gradient)
    upper <- at_upper & is.finite(gradient)
    gradient[lower] <- pmin(gradient[lower], 0)
    gradient[upper] <- pmax(gradient[upper], 0)
  } else if (!is.null(gradient) && !is.null(parameters) && length(gradient) != length(parameters)) {
    gradient[] <- NA_real_
  }
  list(gradient = gradient, active = active)
}

.check_tam_structure <- function(dat, parameter_map = list()) {
  out <- data.frame(check = character(), status = character(), detail = character())
  designs <- c(q_modmat = "index", sd_index_modmat = "index",
               sd_catch_modmat = "catch", F_modmat = "catch", M_modmat = "weight")
  parameters <- c(q_modmat = if (identical(dat$index_settings$q_link, "logit")) "logit_q" else "log_q",
                  sd_index_modmat = "log_sd_index", sd_catch_modmat = "log_sd_catch",
                  F_modmat = "log_mu_f", M_modmat = "mu_m")
  for (nm in names(designs)) {
    x <- dat[[nm]]
    d <- dat$obs[[designs[[nm]]]]
    if (!is.matrix(x) || is.null(d)) next
    mapped <- parameter_map[[parameters[[nm]]]]
    if (!is.null(mapped)) x <- x[, !is.na(mapped), drop = FALSE]
    historical <- if (is.null(d$is_proj)) rep(TRUE, nrow(d)) else !d$is_proj
    if (nm %in% c("q_modmat", "sd_index_modmat", "sd_catch_modmat")) {
      historical <- historical & is.finite(d$obs) & d$obs > 0
    }
    x <- x[historical, , drop = FALSE]
    if (nm == "q_modmat" && !is.null(dat$q_mono_modmat)) {
      x <- cbind(x, dat$q_mono_modmat[historical, , drop = FALSE])
    }
    if (nm == "F_modmat" && dat$F_settings$process %in% c("rw", "cor_rw") && ncol(x)) {
      ny <- sum(!dat$is_proj)
      x <- vapply(seq_len(ncol(x)), function(j) {
        z <- matrix(x[, j], ny, length(dat$ages))
        as.vector(z[-1, , drop = FALSE] - z[-ny, , drop = FALSE])
      }, numeric((ny - 1L) * length(dat$ages)))
    }
    assessable <- nrow(x) > 0L && all(is.finite(x))
    rank <- if (!ncol(x)) 0L else if (assessable) {
      scale <- sqrt(colSums(x^2))
      scale[scale == 0] <- 1
      qr(sweep(x, 2, scale, "/"))$rank
    } else NA_integer_
    out <- rbind(out, data.frame(check = nm,
      status = if (!assessable) "not_assessed" else if (rank == ncol(x)) "pass" else "fail",
      detail = if (!assessable) "No finite informative design available." else
        sprintf("Rank %d of %d fitted columns.", rank, ncol(x))))
  }
  out
}

.check_tam_residuals <- function(obs_pred) {
  out <- data.frame(type = character(), survey = character(), age = numeric(),
                    n = integer(), mean = numeric(), sd = numeric(),
                    max_abs = numeric(), lag1 = numeric())
  for (type in c("catch", "index")) {
    d <- obs_pred[[type]]
    if (is.null(d) || !all(c("std_res", "year", "age") %in% names(d))) next
    historical <- if (is.null(d$is_proj)) rep(TRUE, nrow(d)) else !d$is_proj
    observed <- if (is.null(d$obs)) rep(TRUE, nrow(d)) else is.finite(d$obs) & d$obs > 0
    d <- d[historical & observed & is.finite(d$std_res), , drop = FALSE]
    if (!nrow(d)) next
    if (is.null(d$survey)) d$survey <- ""
    groups <- split(d, interaction(d$survey, d$age, drop = TRUE))
    for (z in groups) {
      z <- z[order(z$year), , drop = FALSE]
      adjacent <- which(diff(z$year) == 1L)
      lag1 <- if (length(adjacent) >= 5L &&
                   stats::sd(z$std_res[adjacent]) > 0 &&
                   stats::sd(z$std_res[adjacent + 1L]) > 0) {
        stats::cor(z$std_res[adjacent], z$std_res[adjacent + 1L])
      } else NA_real_
      out <- rbind(out, data.frame(type = type, survey = as.character(z$survey[1]),
        age = z$age[1], n = nrow(z), mean = mean(z$std_res), sd = stats::sd(z$std_res),
        max_abs = max(abs(z$std_res)), lag1 = lag1))
    }
  }
  out
}

.check_tam_data_summary <- function(dat) {
  out <- data.frame(type = character(), survey = character(), observed = integer(),
                    zero = integer(), missing = integer(), filled = integer(),
                    years = integer())
  for (type in c("catch", "index")) {
    d <- dat$obs[[type]]
    if (is.null(d) || !all(c("obs", "year") %in% names(d))) next
    if (!is.null(d$is_proj)) d <- d[!d$is_proj, , drop = FALSE]
    if (is.null(d$survey)) d$survey <- ""
    for (z in split(d, d$survey)) {
      observed <- is.finite(z$obs) & z$obs > 0
      out <- rbind(out, data.frame(type = type, survey = as.character(z$survey[1]),
        observed = sum(observed), zero = sum(z$obs == 0, na.rm = TRUE),
        missing = sum(is.na(z$obs)),
        filled = if (isTRUE(dat[[paste0(type, "_settings")]]$fill_missing)) sum(!observed) else 0L,
        years = length(unique(z$year[observed]))))
    }
  }
  out
}

.q_process_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  if (!length(dat$q_terms)) return(out)
  d <- dat$obs$index
  if (!is.data.frame(d)) return(out)
  observed <- !d$is_proj & is.finite(d$obs) & d$obs > 0
  for (term in dat$q_terms) {
    if (term$type == "logistic") next
    if (!is.null(term$sd_parameter) || !is.null(term$phi_parameter)) {
      for (group in term$groups) {
        rows <- group$rows[observed[group$rows] & term$multiplier[group$rows] != 0]
        levels <- length(unique(d[[term$variable]][rows]))
        if (levels < 5L) out <- rbind(out, data.frame(issue = "q_effect_support",
          detail = paste0(term$id, ", group ", group$group, ": ", length(rows),
            " observations cover only ", levels, " effect levels. This group's effect has limited observed support; inspect its uncertainty or consider a simpler grouping/formula or an externally supported SD.")))
      }
    }
    independent <- term$type == "iid" || (term$type == "ar1" && identical(term$phi, 0))
    if (!independent || is.null(term$sd_parameter) ||
        !identical(dat$index_settings$q_link, "logit")) next
    z <- .q_term_design(term, nrow(d))[observed, , drop = FALSE]
    z <- z[, colSums(abs(z)) > 0, drop = FALSE]
    s <- dat$sd_index_modmat[observed, , drop = FALSE]
    active <- rowSums(abs(z)) > 0
    if (ncol(z) && all(colSums(z != 0) == 1L) && ncol(s) &&
        any(s[active, , drop = FALSE] != 0)) {
      out <- rbind(out, data.frame(issue = "q_variance_separation",
        detail = paste0(term$id, ": each effect level has only one informative observation, and observation SD is also estimated. The logit link does not guarantee that these sources of variation can be separated. Use replicated levels, simplify the effect or supply one SD.")))
    }
  }
  out
}

.formula_uncertainty_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  sdr <- fit[["sdrep"]]
  if (!length(.formula_terms(fit$dat)) || !is.list(sdr) || !isTRUE(sdr$pdHess)) return(out)
  p <- sdr$par.fixed
  covariance <- sdr$cov.fixed
  if (!is.matrix(covariance) || nrow(covariance) != length(p)) return(out)
  variance <- diag(covariance)
  variance[!is.finite(variance) | variance < 0] <- NA_real_
  half_width <- stats::qnorm(.975) * sqrt(variance)
  wide_sd <- character()
  for (component in c("q", "F", "M", "sd_catch", "sd_index")) {
    terms <- fit$dat[[paste0(component, "_terms")]]
    if (!length(terms)) next
    sd_names <- c(unlist(lapply(terms, `[[`, "sd_parameter")),
      if (component %in% c("F", "M") && identical(fit$dat[[paste0(component, "_settings")]]$process, "iid"))
        paste0("log_sd_", tolower(component)))
    for (term in terms) {
      i <- match(term$phi_parameter, names(p))
      if (!length(i) || is.na(i) || !is.finite(half_width[i]) || !is.finite(p[i])) next
      width <- diff(stats::plogis(p[i] + c(-1, 1) * half_width[i]))
      if (width > .5) out <- rbind(out, data.frame(issue = "formula_AR1_uncertainty",
        detail = paste0(term$id, ": the 95% AR1 correlation interval spans more than 0.5. ",
          "Persistence is weakly estimated; compare with IID/RW effects or a scientifically supported fixed phi before interpreting it.")))
    }
    i <- match(sd_names, names(p))
    i <- i[!is.na(i)]
    wide <- is.finite(half_width[i]) & is.finite(p[i]) & 2 * half_width[i] > log(10)
    wide_sd <- c(wide_sd, names(p)[i[wide]])
  }
  if (length(wide_sd)) out <- rbind(out, data.frame(issue = "formula_SD_uncertainty",
    detail = paste0("95% intervals span more than tenfold for ", paste(unique(wide_sd), collapse = ", "),
      ". These variance components are weakly estimated; compare simpler formulas or use externally supported process SDs.")))
  d <- fit$dat$obs$index
  for (term in Filter(function(x) x$type == "logistic", fit$dat$q_terms)) {
    indices <- which(names(p) == paste0("q_a50_", term$id))
    if (length(indices) != length(term$groups)) next
    for (g in seq_along(term$groups)) {
      group <- term$groups[[g]]
      rows <- group$rows[!d$is_proj[group$rows] & is.finite(d$obs[group$rows]) & d$obs[group$rows] > 0]
      if (!length(rows)) next
      support <- range(d[[term$variable]][rows])
      i <- indices[g]
      if (is.finite(p[i]) && is.finite(half_width[i]) &&
          (p[i] < support[1] || p[i] > support[2]) && 2 * half_width[i] > diff(support)) {
        out <- rbind(out, data.frame(issue = "logistic_support",
          detail = paste0(term$id, ", group ", group$group,
            ": the fitted midpoint lies outside the observed range (", paste(support, collapse = " to "),
            "), and its 95% interval is wider than that range. The data may not distinguish maximum catchability from the age/size curve. Consider a simpler curve or observations covering its rising part and plateau.")))
      }
    }
  }
  out
}

.check_tam_advisories <- function(fit, active, residuals, data) {
  out <- rbind(.mean_process_advisories(fit$dat), .q_process_advisories(fit$dat),
               .sd_process_advisories(fit$dat),
               .rec_process_advisories(fit$dat), .rec_fit_advisories(fit),
               .formula_uncertainty_advisories(fit), .cor_rw_advisories(fit),
               .process_sd_advisories(fit$dat), .process_sd_fit_advisories(fit))
  add <- function(issue, detail) {
    out <<- rbind(out, data.frame(issue = issue, detail = detail))
  }
  if (any(active)) add("active_bounds", paste(sum(active),
    "parameter(s) at a bound; symmetric uncertainty intervals may be inappropriate."))
  sdr <- fit[["sdrep"]]
  p <- if (is.list(sdr)) sdr$par.fixed else NULL
  if (!is.null(p)) {
    small_sd <- startsWith(names(p), "log_sd_") & exp(p) < 1e-6
    high_phi <- startsWith(names(p), "logit_phi_") & stats::plogis(p) > .98
    if (any(small_sd)) add("small_process_sd", "One or more SD estimates are below 1e-6; examine support for those variance components.")
    if (any(high_phi)) add("high_process_correlation", "One or more AR1 correlations exceed 0.98; examine the mean and temporal process together.")
  }
  covariance <- if (is.list(sdr)) sdr$cov.fixed else NULL
  if (is.matrix(covariance) && all(is.finite(covariance)) &&
      all(diag(covariance) > 0) && nrow(covariance) > 1L) {
    correlation <- stats::cov2cor(covariance)
    pairs <- which(upper.tri(correlation) & abs(correlation) > .95, arr.ind = TRUE)
    if (nrow(pairs)) add("parameter_correlation", paste(nrow(pairs),
      "fixed-parameter pair(s) have absolute correlation above 0.95; inspect detailed checks."))
  }
  if (identical(fit$dat$index_settings$q_link, "logit") &&
      any(fit$obs_pred$index$q > 1 - 1e-6, na.rm = TRUE)) {
    add("q_boundary", "Some q estimates are within 1e-6 of one; inspect their uncertainty and information at those ages.")
  }
  if (nrow(data) && any(data$years < 5L)) {
    add("sparse_series", "Some observation series have fewer than five observed years.")
  }
  if (nrow(data) && any(data$filled > data$observed)) {
    add("filled_observations", "Some series contain more filled than observed values; filled values do not provide additional observations.")
  }
  if (identical(fit$dat$N_settings$init, "random") && length(fit$dat$ages) < 10L) {
    add("initial_abundance", "Random initial-abundance variability is estimated from fewer than nine age transitions.")
  }
  supported <- residuals$n >= 10L
  if (any(supported & abs(residuals$mean) > .5)) {
    add("residual_location", "Some age/series groups have conditional mean residuals outside -0.5 to 0.5.")
  }
  if (any(supported & (residuals$sd < .5 | residuals$sd > 1.5), na.rm = TRUE)) {
    add("residual_spread", "Some age/series groups have conditional residual SDs outside 0.5 to 1.5; fitted-state shrinkage also affects spread.")
  }
  if (any(residuals$max_abs > 3)) add("residual_outliers", "Some conditional residuals exceed three SDs in absolute value.")
  if (any(abs(residuals$lag1) > .5, na.rm = TRUE)) {
    add("residual_correlation", "Some age/series groups have absolute adjacent-year residual correlation above 0.5.")
  }
  out
}

.check_tam_curvature <- function(covariance, labels) {
  if (!is.matrix(covariance) || !nrow(covariance) || any(!is.finite(covariance))) {
    return(list(status = "not_assessed"))
  }
  eig <- eigen(covariance, symmetric = TRUE)
  if (any(eig$values <= 0)) return(list(status = "fail", covariance_eigenvalues = eig$values))
  if (is.null(labels)) labels <- paste0("parameter_", seq_len(nrow(covariance)))
  labels <- make.unique(labels)
  correlation <- stats::cov2cor(covariance)
  pairs <- which(upper.tri(correlation) & abs(correlation) > .95, arr.ind = TRUE)
  combinations <- data.frame(parameter = labels, loading = eig$vectors[, 1])
  combinations <- combinations[order(abs(combinations$loading), decreasing = TRUE), , drop = FALSE]
  list(status = "pass", hessian_condition = max(eig$values) / min(eig$values),
       hessian_eigenvalues = 1 / eig$values,
       weak_direction = utils::head(combinations, 5L),
       correlations = data.frame(parameter1 = labels[pairs[, 1]],
         parameter2 = labels[pairs[, 2]], correlation = correlation[pairs]))
}

#' @rdname check_tam
#' @param x A `tam_check` object.
#' @param ... Unused.
#' @export
print.tam_check <- function(x, ...) {
  inform <- function(text, symbol = "i") {
    cli::cli_inform(stats::setNames("{text}", symbol))
  }
  inform(paste("Numerical convergence:", if (x$is_converged) "passed" else "not passed"),
         if (x$is_converged) "v" else "x")
  inform(paste0("Optimizer: ", if (is.null(x$optimizer_code)) "unavailable" else x$optimizer_code,
                if (!is.null(x$optimizer_message)) paste0(" (", x$optimizer_message, ")")))
  inform(sprintf("Max |gradient|: %s adjusted; %s raw (tolerance %s)",
                 signif(x$max_gradient, 3), signif(x$raw_max_gradient, 3), x$grad_tol))
  inform(paste("Positive-definite Hessian:", if (is.null(x$pd_hessian)) "not assessed" else
    if (isTRUE(x$pd_hessian)) "yes" else "no"),
    if (is.null(x$pd_hessian)) "i" else if (isTRUE(x$pd_hessian)) "v" else "x")
  inform(paste("Structural checks:", switch(x$structural_status,
    no_known_issues = "no known redundancies detected; biological identifiability is not established",
    issues = "issues detected", not_assessed = "not assessed")),
    if (x$structural_status == "issues") "x" else "i")
  failures <- x$numerical[x$numerical$status != "pass", , drop = FALSE]
  for (detail in failures$detail) inform(detail, "x")
  for (detail in x$structural$detail[x$structural$status == "fail"]) inform(detail, "x")
  inform(paste("Advisory findings:", nrow(x$advisories)))
  for (detail in x$advisories$detail) inform(detail)
  invisible(x)
}
