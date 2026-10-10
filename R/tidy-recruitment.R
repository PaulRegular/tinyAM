#' Tidy recruitment relationships and residuals
#'
#' @description
#' Collect recruitment predictions, process residuals, and (when specified)
#' parent SSB/recruitment pairs and the fitted stock-recruit curve. Early fixed
#' boundary recruitment states are excluded from process diagnostics.
#'
#' @param fit A fitted `tam_fit` object.
#' @param interval Confidence level, default `0.95`.
#' @return A named list of tables: `mean`, `prediction`, `residual`, `innovation`,
#'   and, for stock-recruit models, `pairs`, `curve`, and `reference`.
#' @details Recruitment and curve estimates are medians with log-scale SEs.
#' Residuals are signed log-scale deviations with SEs on their reported scale.
#' RW residuals are increments; AR1 residuals are deviations from the baseline,
#' while innovations remove the preceding residual's predicted contribution.
#' The curve is evaluated at median historical numeric covariates and reference
#' factor levels. `reference` identifies those values. Curve intervals use the
#' exact curve derivatives and RTMB's fixed-parameter covariance (the delta
#' method); pair/prediction/residual uncertainty comes from ADREPORT. Missing
#' uncertainty remains `NA`. Curve intervals do not include new process noise.
#' @examples
#' \dontrun{
#' fit <- fit_tam(cod_obs, N_settings = list(process = "off",
#'   rec_form = ~ bh(ssb) + iid(year)))
#' tidy_recruitment(fit)$pairs
#' }
#' @seealso [recruitment_formulas], [tidy_tam()], [vis_tam()]
#' @export
tidy_recruitment <- function(fit, interval = .95) {
  fit <- .require_tam_fit(fit, arg = "fit")
  rec <- fit$dat$rec
  if (is.null(rec)) return(list())
  z <- stats::qnorm(.5 + interval / 2)
  sdr <- fit[["sdrep"]]
  errors <- if (is.null(sdr)) list() else as.list(sdr, "Std. Error", report = TRUE)
  table <- function(name, logarithmic = FALSE) {
    value <- unname(fit$rep[[name]])
    se <- errors[[name]]
    if (is.null(se)) se <- rep(NA_real_, length(value))
    d <- data.frame(year = fit$dat$years[rec$eligible], est = value,
      lwr = value - z * se, upr = value + z * se, se = unname(se),
      se_scale = if (logarithmic) "log" else "reported",
      is_proj = fit$dat$is_proj[rec$eligible])
    if (logarithmic) d <- trans_est(d, exp)
    d
  }
  out <- list(mean = table("rec_log_mean", TRUE),
    prediction = table("rec_log_prediction", TRUE),
    residual = table("rec_residual"), innovation = table("rec_innovation"))
  if (!is.null(rec$curve)) {
    parent <- table("rec_log_parent", TRUE)
    R <- tidy_pop(fit, interval)$recruitment
    R <- R[match(parent$year, R$year), ]
    out$pairs <- data.frame(year = parent$year,
      parent_year = parent$year - rec$curve$lag,
      ssb = parent$est, ssb_lwr = parent$lwr, ssb_upr = parent$upr,
      recruitment = R$est, recruitment_lwr = R$lwr, recruitment_upr = R$upr,
      mean = out$mean$est, prediction = out$prediction$est,
      residual = out$residual$est, is_proj = parent$is_proj)
    curve <- .rec_curve_table(fit, interval)
    out$curve <- curve$curve
    out$reference <- curve$reference
  }
  attr(out, "interval") <- interval
  out
}

.rec_curve_table <- function(fit, interval) {
  rec <- fit$dat$rec
  p <- .tam_parameter_summary(fit, "Estimate")
  historic <- rec$data[!fit$dat$is_proj, , drop = FALSE]
  reference <- historic[1L, , drop = FALSE]
  variables <- all.vars(rec$fixed_form)
  for (name in variables) {
    value <- historic[[name]]
    reference[[name]] <- if (is.numeric(value)) stats::median(value) else {
      levels <- if (is.factor(value)) levels(value) else sort(unique(value))
      factor(levels[1L], levels = levels, ordered = is.ordered(value))
    }
  }
  X <- stats::model.matrix(rec$fixed_form, reference)
  X <- X[, colnames(rec$matrix), drop = FALSE]
  offset <- if (ncol(X)) drop(X %*% p$rec_beta) else 0
  parent <- exp(fit$rep$rec_log_parent[!fit$dat$is_proj[rec$eligible]])
  S <- seq(.05 * min(parent), 1.2 * max(parent), length.out = 100L)
  value <- as.numeric(.rec_log_curve(log(S), p, rec$curve$type) + offset)
  se <- rep(NA_real_, length(S))
  sdr <- fit[["sdrep"]]
  if (is.list(sdr) && isTRUE(sdr$pdHess)) {
    covariance <- sdr$cov.fixed
    indices <- c(which(names(sdr$par.fixed) == "log_sr_alpha"),
      which(names(sdr$par.fixed) == "log_sr_beta"), which(names(sdr$par.fixed) == "rec_beta"))
    density <- exp(p$log_sr_beta) * S
    gradient <- cbind(1, if (rec$curve$type == "bh") -density / (1 + density) else -density,
                      matrix(rep(X, each = length(S)), length(S), ncol(X)))
    if (length(indices) == ncol(gradient) && all(is.finite(covariance[indices, indices, drop = FALSE]))) {
      se <- sqrt(pmax(0, rowSums((gradient %*% covariance[indices, indices, drop = FALSE]) * gradient)))
    }
  }
  z <- stats::qnorm(.5 + interval / 2)
  list(curve = data.frame(ssb = S, est = exp(value), lwr = exp(value - z * se),
    upr = exp(value + z * se), se = se, se_scale = "log"),
    reference = data.frame(covariate = variables,
      value = vapply(variables, function(name) as.character(reference[[name]]), character(1))))
}
