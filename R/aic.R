#' Compare fitted models using information criteria
#'
#' Compare models fitted to the same catch and survey observations. Smaller
#' values indicate a better balance between fit and model complexity.
#' Marginal AIC assesses the model after integrating over its random effects.
#' Conditional AIC (cAIC) assesses new observations sharing the fitted population
#' and other random effects. Neither replaces hindcast evaluation of forecasting
#' skill.
#'
#' @param object A fitted [fit_tam()] model.
#' @param ... Additional `tam_fit` models. For `logLik()`, currently unused.
#' @param k Penalty per parameter for marginal AIC; the usual AIC uses `2`.
#' @param details Return a table including conditional log-likelihood, parameter
#'   counts and effective degrees of freedom (EDF)?
#'
#' @details
#' **Likelihood scale.** Densities are evaluated on log observations, as in
#' [nll_fun()]. No natural-scale lognormal Jacobian is added. Absolute values
#' should not be compared with criteria calculated on a different response
#' scale. Multiple-model calls require identical included observations and units.
#'
#' **Marginal AIC.** `logLik()` returns the optimized Laplace approximation to
#' the marginal log-likelihood. Its `df` counts all estimated outer parameters,
#' including process SDs and correlations, after parameter maps are applied.
#' Integrated latent states are not included in that count.
#'
#' **Conditional AIC.** Equation 6 of Zheng, Cadigan and Thorson (2026) gives
#' \deqn{\mathrm{cAIC}=-2\ell_c+2p_c+2q-2\mathrm{tr}(H_j^{-1}H_r).}
#' Here \eqn{\ell_c} is the observation log-likelihood at the fitted conditional
#' mode, \eqn{q} is the number of integrated random coordinates, and \eqn{H_j}
#' and \eqn{H_r} are random-coordinate Hessians of the joint and process/prior
#' negative log-likelihoods. RTMB supplies derivatives; sparse Cholesky solves
#' calculate the trace without forming a dense inverse.
#'
#' `edf_fixed` is \eqn{p_c}, the number of estimated outer parameters entering
#' the conditional observation likelihood, identified using structural
#' dependencies of observation predictions and SDs. It need not equal
#' `n_fixed`: process-only parameters are excluded. The supplementary Gaussian
#' example similarly excludes a prior-only parameter. This differs from the
#' all-outer-parameter convention currently used by sdmTMB and WHAM.
#' Counting follows tinyAM's fitted coordinates. Reparameterizations that move
#' variance parameters into observation predictions can change this approximation.
#'
#' `edf_random` equals \eqn{q-\mathrm{tr}(H_j^{-1}H_r)} and `EDF` equals
#' `edf_fixed + edf_random`. Random EDF is partitioned by fitted coordinates:
#' recruitment (`log_r`), N (`log_n`), F (`log_f`), M (`log_m`), and other
#' coordinates, including latent initial abundance and Gaussian formula effects.
#' These contributions sum to `edf_random`. Their allocation depends on the
#' parameterization; they are accounting contributions rather than independent
#' measures of support. Unusual contributions are warned about, not truncated.
#'
#' **Limitations.** Intrinsic F/M/N random walks and latent initial abundance
#' have unpenalized directions. Criteria are returned with warnings: marginal
#' likelihood normalization is arbitrary, and comparisons between different
#' unpenalized structures are particularly uncertain. Recruitment and formula
#' random walks anchored at a fixed starting value are not intrinsically
#' improper. Equation 6 assumes jointly Gaussian random effects; nonlinear
#' stock-recruit and active cohort-survival constructions can violate this
#' assumption despite Gaussian innovations. Their cAIC is an exploratory local
#' approximation. Parameter boundaries and poor numerical convergence also
#' limit the usual approximations. A singular or indefinite joint random
#' Hessian prevents cAIC calculation; a singular process Hessian alone does not.
#'
#' @return `logLik()` returns a `logLik` object with `df` and `nobs` attributes.
#'   `AIC()` returns a number for one fit, or a standard `df`/`AIC` table for
#'   multiple fits. `cAIC()` returns a number for one fit, or a table sorted by
#'   increasing cAIC with `EDF`, `cAIC` and `delta`. With `details = TRUE`, it
#'   also includes conditional `logLik`, `n_fixed`, `n_random`, `edf_fixed`,
#'   `edf_rec`, `edf_N`, `edf_F`, `edf_M`, `edf_other` and `edf_random`.
#'   Single-model detailed output has `delta = 0`.
#'
#' @references
#' Zheng, N., Cadigan, N. G. and Thorson, J. T. (2026).
#' A unifying theory for the evaluation of conditional Akaike information for
#' mixed-effects models. *Journal of Statistical Computation and Simulation*.
#' \doi{10.1080/00949655.2026.2736799}. Equation 6 and supplementary R/TMB code.
#'
#' @seealso [check_tam()], [fit_hindcast()]
#'
#' @examples
#' \dontrun{
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' alternative <- update(fit, F_settings = list(process = "ar1", mu_form = ~1))
#' AIC(fit, alternative)
#' cAIC(fit, alternative)
#' cAIC(fit, details = TRUE)
#' }
#' @name information_criteria
NULL

.aic_validate_fit <- function(object) {
  if (!inherits(object, "tam_fit") || !is.list(object$opt) ||
      length(object$opt$objective) != 1L || !is.finite(object$opt$objective)) {
    cli::cli_abort("Information criteria require a fitted {.cls tam_fit} with a finite optimized likelihood.")
  }
  invisible(TRUE)
}

.aic_warnings <- function(object, conditional = FALSE) {
  dat <- object$dat
  intrinsic <- vapply(c("N", "F", "M"), function(x) {
    dat[[paste0(x, "_settings")]]$process %in% c("rw", "cor_rw")
  }, logical(1))
  if (any(intrinsic) || identical(dat$N_settings$init, "latent")) {
    cli::cli_warn(c("Information criterion uses an improper process distribution.",
      "i" = "Unpenalized process levels or latent initial abundances make marginal normalization arbitrary; cAIC is exploratory. Comparisons across different unpenalized structures are uncertain."
    ))
  }
  if (conditional && (!is.null(dat$rec$curve) || dat$N_settings$process != "off")) {
    cli::cli_warn(c("cAIC is exploratory for this nonlinear latent-state distribution.",
      "i" = "Gaussian process errors do not guarantee the jointly Gaussian random effects assumed by Equation 6."
    ))
  }
  if (!isTRUE(object$is_converged)) {
    cli::cli_warn("Information criterion is based on a fit that has not passed numerical convergence checks; inspect {.code check_tam(fit)}.")
  }
  lower <- object$bounds$lower
  upper <- object$bounds$upper
  if (length(lower) && any(object$opt$par <= lower + 1e-7 | object$opt$par >= upper - 1e-7)) {
    cli::cli_warn("Estimated parameters are on a bound; ordinary information-criterion penalties are approximate.")
  }
}

.aic_observations <- function(object) {
  dat <- object$dat
  cols <- c("type", "survey", "year", "age", "samp_time", "unit")
  d <- dat$obs_map[dat$is_observed, intersect(cols, names(dat$obs_map)), drop = FALSE]
  for (nm in setdiff(cols, names(d))) d[[nm]] <- NA_character_
  d <- d[, cols, drop = FALSE]
  d[] <- lapply(d, as.character)
  d$log_obs <- as.numeric(dat$observed)
  d <- d[do.call(order, c(unname(d[cols]), list(na.last = TRUE))), , drop = FALSE]
  rownames(d) <- NULL
  d
}

.aic_models <- function(object, dots, labels) {
  models <- c(list(object), dots)
  if (any(!vapply(models, inherits, logical(1), "tam_fit"))) {
    cli::cli_abort("All models must be fitted {.cls tam_fit} objects.")
  }
  lapply(models, .aic_validate_fit)
  reference <- .aic_observations(object)
  for (model in models[-1L]) {
    if (!identical(reference, .aic_observations(model))) {
      cli::cli_abort("Models must use the same included catch and survey observations, identities and units. Compare models fitted to a common dataset.")
    }
  }
  names(models) <- make.unique(labels)
  models
}

#' @rdname information_criteria
#' @export
logLik.tam_fit <- function(object, ...) {
  .aic_validate_fit(object)
  .aic_warnings(object)
  structure(-as.numeric(object$opt$objective), df = length(object$opt$par),
    nobs = sum(object$dat$is_observed), class = "logLik")
}

#' @rdname information_criteria
#' @export
AIC.tam_fit <- function(object, ..., k = 2) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) cli::cli_abort("{.arg k} must be one finite number.")
  expressions <- as.list(substitute(list(object, ...)))[-1L]
  models <- .aic_models(object, list(...), vapply(expressions, function(x) paste(deparse(x), collapse = ""), character(1)))
  ll <- lapply(models, stats::logLik)
  result <- data.frame(df = vapply(ll, attr, numeric(1), "df"),
    AIC = vapply(ll, function(x) -2 * as.numeric(x) + k * attr(x, "df"), numeric(1)),
    row.names = names(models))
  if (length(models) == 1L) result$AIC else result
}

# Solve in column batches so memory does not grow as a dense q-by-q inverse.
# This covers the entire prior pattern, including entries outside the joint
# Hessian's stored sparsity pattern.
.aic_trace <- function(joint, prior) {
  n <- nrow(joint)
  if (!n) return(numeric())
  if (!identical(dim(joint), dim(prior)) ||
      any(!is.finite(joint@x)) || any(!is.finite(prior@x))) {
    cli::cli_abort("cAIC requires finite, aligned joint and process Hessians.")
  }
  tryCatch({
    factor <- Matrix::Cholesky(joint, LDL = FALSE)
    pivots <- Matrix::diag(Matrix::expand(factor)$L)
    if (min(pivots) < 1e-6 * max(pivots)) {
      cli::cli_warn("Joint-Hessian Cholesky pivots differ greatly; cAIC may be numerically sensitive. Inspect poorly supported random effects.")
    }
    contribution <- numeric(n)
    for (start in seq.int(1L, n, by = 64L)) {
      cols <- seq.int(start, min(n, start + 63L))
      solution <- Matrix::solve(factor, prior[, cols, drop = FALSE], system = "A")
      contribution[cols] <- solution[cbind(cols, seq_along(cols))]
    }
    if (any(!is.finite(contribution))) cli::cli_abort("cAIC sparse solves produced non-finite values.")
    contribution
  }, error = function(e) {
    cli::cli_abort(c("Cannot calculate cAIC: the joint random-effect Hessian must be positive definite and numerically usable.",
      "i" = "Inspect check_tam(fit), poorly supported effects and unidentifiable states.",
      "x" = conditionMessage(e)))
  })
}


.aic_conditional_count <- function(object, component, random) {
  dat <- object$dat
  dat$.likelihood_component <- "design"
  parameters <- component$env$parList(par = object$parameter_values)
  design_object <- RTMB::MakeADFun(function(p) {
    aic_predictors <- nll_fun(p, dat)
    RTMB::ADREPORT(aic_predictors)
    0
  }, parameters, map = object$parameter_map, ADreport = TRUE, silent = TRUE)
  tape <- RTMB::GetTape(design_object)
  design <- methods::as(tape$jacfun(sparse = TRUE)(object$parameter_values), "CsparseMatrix")
  # Stored entries include structural dependencies whose derivatives are zero
  # at this particular mode. Never test the numerical values to count them.
  dependencies <- which(diff(design@p) > 0L)
  sum(!random[dependencies])
}

.caic_fit <- function(object) {
  .aic_warnings(object, conditional = TRUE)
  mode <- object$parameter_values
  random <- object$obj$env$lrandom()
  if (length(mode) != length(random) || any(!is.finite(mode))) {
    cli::cli_abort("Stored fitted parameters are unavailable or inconsistent. Refit the model before calculating cAIC.")
  }
  process <- .fitted_objective(object, "process")
  observation <- .fitted_objective(object, "observation")
  if (!isTRUE(all.equal(process$env$last.par, mode, check.attributes = FALSE)) ||
      !identical(names(process$env$last.par), names(mode))) {
    cli::cli_abort("Mapped parameter ordering differs from the fitted objective; cAIC cannot align the random effects.")
  }
  if (any(random)) {
    joint_hessian <- object$obj$env$spHess(mode, random = TRUE)
    prior_hessian <- process$env$spHess(mode, random = TRUE)
    negative_edf <- .aic_trace(joint_hessian, prior_hessian)
  } else {
    negative_edf <- numeric()
  }
  categories <- c(log_r = "rec", log_n = "N", log_f = "F", log_m = "M")
  group <- unname(categories[names(mode)[random]])
  group[is.na(group)] <- "other"
  edf <- vapply(c("rec", "N", "F", "M", "other"), function(x) {
    sum(1 - negative_edf[group == x])
  }, numeric(1))
  if (any(edf < -1e-6) || any(1 - negative_edf > 1 + 1e-6)) {
    cli::cli_warn("cAIC has unusual random-effect EDF contributions; inspect model curvature and applicability. Contributions have not been truncated.")
  }
  fixed <- .aic_conditional_count(object, process, random)
  ll <- -as.numeric(observation$env$f(mode))
  result <- data.frame(logLik = ll, n_fixed = length(object$opt$par), n_random = sum(random),
    edf_fixed = fixed, edf_rec = edf["rec"], edf_N = edf["N"], edf_F = edf["F"],
    edf_M = edf["M"], edf_other = edf["other"], edf_random = sum(edf),
    EDF = fixed + sum(edf), row.names = NULL)
  result$cAIC <- -2 * ll + 2 * result$EDF
  if (any(!is.finite(unlist(result)))) cli::cli_abort("cAIC calculation produced non-finite likelihood or EDF values.")
  result
}

#' @rdname information_criteria
#' @export
cAIC <- function(object, ..., details = FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details)) cli::cli_abort("{.arg details} must be TRUE or FALSE.")
  expressions <- as.list(substitute(list(object, ...)))[-1L]
  models <- .aic_models(object, list(...), vapply(expressions, function(x) paste(deparse(x), collapse = ""), character(1)))
  results <- do.call(rbind, lapply(models, .caic_fit))
  results$delta <- results$cAIC - min(results$cAIC)
  results <- results[order(results$cAIC), , drop = FALSE]
  if (details) return(results)
  if (length(models) == 1L) return(results$cAIC)
  results[, c("EDF", "cAIC", "delta"), drop = FALSE]
}
