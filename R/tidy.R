
#' Tidy a (named) matrix/array into long format
#'
#' @description
#' Converts a matrix or array `x` into a tidy data frame with one row per cell,
#' one column per dimension (using the dimension names if available), and a
#' value column.
#'
#' @param x A matrix or array. Works with 2D or higher dimensions.
#' @param value_name Character scalar: name of the value column. Default `"x"`.
#' @param require_dimnames Logical; if `TRUE` (default), error if any dimension
#'   of `x` lacks names. If `FALSE`, fallback names like `Var1`, `Var2`, … are used.
#'
#' @details
#' Internally uses [base::as.data.frame.table()] to reshape `x`. After reshaping,
#' all **dimension columns** are passed through [utils::type.convert()] with
#' `as.is = TRUE`, so numeric-like labels (e.g., `"1"`, `"2.5"`) become numeric.
#'
#' @return
#' A data frame with `length(x)` rows, one column per dimension, and a value
#' column named `value_name`.
#'
#' @examples
#' m <- matrix(1:6, nrow = 2, dimnames = list(age = c("2","3"), year = c("2001","2002","2003")))
#' tidy_array(m)
#'
#' a <- array(1:8,
#'            dim = c(2,2,2),
#'            dimnames = list(age = c("2","3"), year = c("2001","2002"), area = c("N","S")))
#' tidy_array(a, value_name = "val")
#'
#' # Permissive mode if dimnames are missing:
#' dimnames(m) <- NULL
#' tidy_array(m, require_dimnames = FALSE)
#'
#' @aliases tidy_mat
#' @export tidy_array tidy_mat
tidy_array <- function(x, value_name = "x", require_dimnames = TRUE) {
  if (!is.matrix(x) && !is.array(x)) {
    cli::cli_abort("{.arg x} must be a matrix or array.")
  }

  nm <- dimnames(x)

  if (require_dimnames) {
    if (is.null(nm) || any(vapply(nm, is.null, logical(1)))) {
      cli::cli_abort("All dimensions must have names (dimnames). Set {.code require_dimnames = FALSE} to allow defaults.")
    }
  }

  # Melt to data.frame; dimension columns come first
  df <- as.data.frame.table(x, responseName = value_name, stringsAsFactors = FALSE)

  k <- length(dim(x))
  if (k > 0) {
    df[seq_len(k)] <- lapply(df[seq_len(k)], type.convert, as.is = TRUE)
  }

  df
}

#' @rdname tidy_array
#' @export
tidy_mat <- tidy_array


#' Tidy observed, predicted, and residual diagnostics
#'
#' @description
#' Extracts observations, predictions, and observation SDs for `catch` and `index`
#' from a fitted TAM object, and adds standardized residuals on the log scale.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level for combined catchability intervals.
#' @param add_osa_res Logical; add one-step-ahead residuals? Hard wired to
#'                    apply the `"oneStepGaussianOffMode"` method.
#'                    See [RTMB::oneStepPredict()] for details.
#' @param ... Arguments to pass to [RTMB::oneStepPredict()].
#'
#' @details
#' `pred` is the conditional median on the natural scale. `sd` is the fitted
#' SD of log observations, not a standard error of the prediction. Standardized
#' residuals are `(log(obs) - log(pred)) / sd`; zero and missing observations
#' have missing residuals. These conditional residuals do not account for the
#' uncertainty in fitted states. Optional one-step-ahead residuals use a
#' different predictive calculation and are not currently supported with filled
#' missing observations. See [tinyAM-model] for the observation equations.
#' When stockassessment is loaded, one-step residuals run sequentially: the
#' upstream TMB parallel helper cannot select between SAM and RTMB libraries.
#'
#' @return
#' A named list with two data frames:
#'
#' - **catch**: original columns plus `pred`, `sd`, and `std_res`.
#' - **index**: original columns plus `pred`, `sd`, `q`, and `std_res`.
#'   Combined catchability intervals are `q_lwr`/`q_upr`; `q_se` is on the
#'   selected link scale named in `q_se_scale`. They include uncertainty across
#'   the full formula, rather than combining separate coefficient SEs.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' obs_pred <- tidy_obs_pred(fit)
#' head(obs_pred$catch)
#' head(obs_pred$index)
#'
#' @seealso [fit_tam()], [tidy_rep()], [tidy_sdrep()], [tidy_pop()]
#' @export

tidy_obs_pred <- function(fit, add_osa_res = FALSE, interval = .95, ...) {
  fit <- .require_tam_fit(fit, arg = "fit")
  if (!is.numeric(interval) || length(interval) != 1L || !is.finite(interval) ||
      interval <= 0 || interval >= 1) {
    cli::cli_abort("{.arg interval} must be one number strictly between zero and one.")
  }

  obs_pred <- fit$dat$obs[c("catch", "index")]
  pred <- split(exp(fit$rep$log_pred), fit$dat$obs_map$type)
  sd <- split(fit$rep$sd_obs, fit$dat$obs_map$type)

  obs_pred$catch$pred <- pred$catch
  obs_pred$catch$sd <- sd$catch
  obs_pred$catch$std_res <- with(obs_pred$catch, ifelse(obs == 0, NA, (log(obs) - log(pred)) / sd))

  obs_pred$index$pred <- pred$index
  obs_pred$index$sd <- sd$index
  obs_pred$index$q <- exp(fit$rep$log_q_obs)
  obs_pred$index$q_lwr <- obs_pred$index$q_upr <- obs_pred$index$q_se <- NA_real_
  obs_pred$index$q_se_scale <- fit$dat$index_settings$q_link
  if (!is.null(fit[["sdrep"]])) {
    estimate <- as.list(fit[["sdrep"]], "Estimate", report = TRUE)$q_link_prediction
    error <- as.list(fit[["sdrep"]], "Std. Error", report = TRUE)$q_link_prediction
    if (!is.null(estimate)) {
      estimate <- as.numeric(estimate)
      error <- as.numeric(error)
      z <- stats::qnorm(.5 + interval / 2)
      inverse <- if (identical(fit$dat$index_settings$q_link, "logit")) stats::plogis else exp
      obs_pred$index$q_lwr <- inverse(estimate - z * error)
      obs_pred$index$q_upr <- inverse(estimate + z * error)
      obs_pred$index$q_se <- error
    }
  }
  obs_pred$index$std_res <- with(obs_pred$index, ifelse(obs == 0, NA, (log(obs) - log(pred)) / sd))

  if (add_osa_res) {
    args <- list(...)
    if (isTRUE(args$parallel) && "stockassessment" %in% loadedNamespaces()) {
      args$parallel <- FALSE
      cli::cli_inform("Calculating one-step residuals sequentially while SAM and RTMB are loaded.")
    }
    osa_res <- do.call(RTMB::oneStepPredict,
      c(list(obj = fit$obj, method = "oneStepGaussianOffMode"), args))
    split_osa_res <- split(osa_res$residual, fit$dat$obs_map$type[fit$dat$obs_map$is_observed])
    split_is_observed <- split(fit$dat$obs_map$is_observed, fit$dat$obs_map$type)
    obs_pred$catch$osa_res <- obs_pred$index$osa_res <- NA
    obs_pred$catch$osa_res[split_is_observed$catch] <- split_osa_res$catch
    obs_pred$index$osa_res[split_is_observed$index] <- split_osa_res$index
  }

  attr(obs_pred, "interval") <- interval
  obs_pred
}


#' Tidy reported trends
#'
#' @description
#' Converts year and age x year objects in `fit$rep` to long (tidy) data frames.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#'
#' @return
#' A named list of data frames (e.g., `N`, `abundance`, `ssb`, `F`, `F_bar`),
#' where each data frame has one column per dimension (e.g., `year`, `age`),
#' a value column `est`, and an `is_proj` column.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' trends <- tidy_rep(fit)
#' names(trends)
#' head(trends$N)
#'
#' @seealso [tidy_obs_pred()], [tidy_sdrep()], [tidy_pop()], [tidy_mat()]
#' @export

tidy_rep <- function(fit) {
  if (.is_tam_fit(fit)) {
    dat <- fit$dat
    rep <- fit$rep
  } else if (is.list(fit) && !is.null(fit$dat) && !is.null(fit$rep)) {
    dat <- fit$dat
    rep <- fit$rep
  } else {
    cli::cli_abort(c(
      "`{.arg fit}` must be either a {.cls tam_fit} or a list with",
      "an {.field dat} element and an {.field rep} element."
    ))
  }

  if (is.null(dat$years) || is.null(dat$is_proj)) {
    cli::cli_abort("`{.arg fit}$dat` must provide `years` and `is_proj` vectors.")
  }

  if (length(dat$years) != length(dat$is_proj)) {
    cli::cli_abort("`{.arg fit}$dat$is_proj` must align with `{.arg fit}$dat$years`.")
  }

  keep <- vapply(rep, function(x) is.matrix(x) || length(x) == length(dat$years),
                 logical(1))
  keep[grepl("^eta_(q|mu_[FM])_increments$", names(keep))] <- FALSE
  rep_items <- rep[keep]

  trends <- lapply(rep_items, function(x) {
    if (is.matrix(x)) {
      d <- tidy_mat(x, value_name = "est")
      d$is_proj <- d$year %in% dat$years[dat$is_proj]
      return(d)
    }

    data.frame(
      year = dat$years,
      est = unname(x),
      is_proj = dat$is_proj,
      stringsAsFactors = FALSE
    )
  })

  trends
}


#' Transform and rescale estimate columns
#'
#' @description
#' Applies a transformation (e.g., `exp`) and optional rescaling to the columns
#' `est`, `lwr`, and `upr` of a data frame.
#' SE columns and their scale labels are left unchanged; this helper does not
#' calculate transformed standard errors.
#'
#' @param data A data frame containing columns `est`, `lwr`, and `upr`.
#' @param transform A function applied to `est`, `lwr`, and `upr`
#'   (set to `NULL` to skip). Default `exp`.
#' @param scale Numeric scale factor by which transformed columns are divided.
#'   Default `1`.
#'
#' @return
#' The input data frame with `est`, `lwr`, and `upr` transformed and rescaled.
#'
#' @examples
#' d <- data.frame(est = log(100), lwr = log(80), upr = log(120))
#' trans_est(d)                # exp + no rescale
#' trans_est(d, transform = NULL, scale = 1000)  # no transform, rescale
#'
#' @keywords internal
#' @export
trans_est <- function(data, transform = exp, scale = 1) {
  if (!is.null(transform)) {
    data[, c("est", "lwr", "upr")] <-  apply(data[, c("est", "lwr", "upr")], 2, transform)
  }
  data[, c("est", "lwr", "upr")] <- data[, c("est", "lwr", "upr")] / scale
  data
}



#' Population time series with confidence intervals
#'
#' @description
#' Returns annual recruitment, abundance, biomass, spawning biomass, and
#' average F and M, with confidence limits and clearly labelled SE scales.
#' Estimates and limits are reported in the units of each population quantity.
#'
#' @details
#' Assumptions:
#'
#' - All ADREPORTED series used here have length equal to `length(fit$dat$years)`.
#' - Estimates are on the log scale and are transformed with `exp` via [trans_est()].
#' - List element names are cleaned by removing a leading `"log_"` prefix.
#'
#' The interval is constructed as `est ± z * se` with
#' `z = qnorm(1 - (1 - interval) / 2)`.
#' Estimates and confidence limits are then exponentiated. The `se` column is
#' the standard error on the log scale, labelled by `se_scale = "log"`.
#' A small log-scale SE approximates the coefficient of variation (CV):
#' for example, `0.10` is approximately 10% relative uncertainty. This
#' approximation becomes less accurate for large SEs; use the confidence limits
#' to describe uncertainty on the reported scale.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level in `(0, 1)`; default `0.95`.
#'
#' @return
#' A named list of data frames (one per series), each with columns:
#'
#' - `year`, `est`, `lwr`, `upr`, `se`, `se_scale`, `is_proj`.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' trends <- tidy_sdrep(fit, interval = 0.9)
#' names(trends)
#' head(trends$ssb)
#'
#' @seealso [trans_est()], [tidy_rep()], [tidy_pop()]
#' @export
tidy_sdrep <- function(fit, interval = 0.95) {
  fit <- .require_tam_fit(fit, arg = "fit")

  if (is.null(fit[["sdrep"]])) {
    series <- intersect(c("recruitment", "abundance", "biomass", "ssb", "F_bar", "M_bar"),
                        names(fit$rep))
    return(stats::setNames(lapply(series, function(nm) {
      data.frame(year = fit$dat$years, est = unname(fit$rep[[nm]]),
        lwr = NA_real_, upr = NA_real_, se = NA_real_, se_scale = "log",
        is_proj = fit$dat$is_proj)
    }), series))
  }

  ## assumes all ADREPORTED objects are equal length to years and are in log space
  vals <- as.list(fit[["sdrep"]], "Estimate", report = TRUE)
  ses <- as.list(fit[["sdrep"]], "Std. Error", report = TRUE)
  vals <- vals[!grepl("^(q_link_prediction|eta_(q|mu_[FM])_increments)$", names(vals))]
  ses <- ses[names(vals)]
  df <- lapply(seq_along(vals), function(i) {
    d <- data.frame(year = fit$dat$years,
                    est = vals[[i]],
                    lwr = vals[[i]] - qnorm(1 - ((1 - interval) / 2)) * ses[[i]],
                    upr = vals[[i]] + qnorm(1 - ((1 - interval) / 2)) * ses[[i]],
                    se = ses[[i]],
                    se_scale = "log",
                    is_proj = fit$dat$is_proj) |>
      trans_est(transform = exp)
  })
  names(df) <- gsub("log_", "", names(vals))
  df
}

#' Collect population trends and age-specific estimates
#'
#' @description
#' Combines population time series with uncertainty from [tidy_sdrep()] and
#' age-specific estimates and other reported quantities from [tidy_rep()]
#' into a single named list for downstream plotting and summaries.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level for intervals passed to [tidy_sdrep()]; default `0.95`.
#'
#' @return
#' A named list containing the elements returned by [tidy_sdrep()] and
#' [tidy_rep()] (names preserved).
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' pop <- tidy_pop(fit)
#' names(pop)
#'
#' @seealso [tidy_sdrep()], [tidy_rep()], [tidy_obs_pred()]
#' @export
tidy_pop <- function(fit, interval = 0.95) {
  fit <- .require_tam_fit(fit, arg = "fit")

  sdrep_trends <- tidy_sdrep(fit, interval = interval)
  rep_trends <- tidy_rep(fit)
  not_in_sdrep <- setdiff(names(rep_trends), names(sdrep_trends))
  out <- c(sdrep_trends, rep_trends[not_in_sdrep])
  attr(out, "interval") <- interval
  out
}


#' Tidy parameter estimates (fixed & random) with CIs and back-transforms
#'
#' @description
#' Summarizes fitted parameters with estimates, confidence limits, and clearly
#' labelled SE scales. Use this to inspect process variability, catchability
#' effects, and latent population states. Actual catchability for each survey
#' observation is available from [tidy_obs_pred()].
#'
#' @details
#' Combines estimates (`Estimate`) and standard errors (`Std. Error`) from
#' `fit$sdrep`. Parameters whose names begin with `log_` or `logit_` are
#' back-transformed to the natural scale:
#'
#' - `log_`  → `exp()` (and the `log_` prefix is dropped, e.g. `log_sd_r` → `sd_r`)
#' - `logit_` → `plogis()` (and the `logit_` prefix is dropped, e.g. `logit_phi_f` → `phi_f`)
#'
#' An exception is `logit_q`: formula coefficients, confidence limits and SEs
#' remain on the fitted logit scale (`se_scale = "logit"`). A slope or contrast
#' is not an absolute q and must not be inverse-logit transformed on its own.
#' Actual observation-specific q is reported by [tidy_obs_pred()].
#'
#' Estimates and confidence limits are shown on the reported scale. SEs remain
#' on the fitted scale, identified by `se_scale`: `"log"`, `"logit"`, or
#' `"reported"` (the same scale as the displayed estimate). A small log-scale
#' SE approximates a CV: `0.10` is approximately 10% relative uncertainty.
#' This interpretation does not apply to logit or reported-scale SEs, and is
#' unreliable for large log-scale SEs. Confidence limits are the preferred
#' summary of uncertainty on the reported scale.
#'
#' - For [mono()] catchability terms, `dq` is a non-negative step on the
#'   selected q-link scale (log-q or logit-q). Estimates and SEs are reported
#'   directly on this same scale,
#'   with untransformed Wald intervals (which may cross zero at a boundary).
#'   These local curvature SEs are not boundary-adjusted inference.
#'   Coefficient names identify
#'   transitions and groups. Actual observation-specific q is in [tidy_obs_pred()].
#'
#' Formula coefficients named `log_mu_f`, `log_sd_*`, or `log_q` are
#' exponentiated in summaries. An exponentiated slope is a multiplicative
#' change per unit covariate, not the fitted F, SD, or q surface itself.
#' An exception is the mean-\eqn{M} formula
#' coefficients (`mu_m`), which operate on the log scale but are named without a
#' `log_` prefix to reflect that their values may be positive or negative; they
#' therefore print on the fitted log scale.
#'
#' Fixed-effect parameters are returned in a single data frame (`$fixed`);
#' random-effect parameters are returned as a named list of data frames
#' (`$random`), one per random block (e.g. `log_f`, `log_r`, `missing`, …).
#'
#' Labels are added where applicable:
#' - For parameters specified using a formula in [prepare_tam()] (e.g., `log_q`, `logit_q`,
#'   `log_sd_catch`, `log_sd_index`), a `coef` column is added.
#' - For `log_r`, `year` contains years 2:Y; full recruitment is in [tidy_pop()].
#' - For `log_n0`, an `age` column identifies the initial older-age state.
#'   `log_r0` and `log_n0` are exponentiated to abundance levels; `log_sd_n0`
#'   is exponentiated to the initial-age residual SD.
#' - For matrices (e.g., `log_f`, `log_n`), `year` and/or `age` columns are added
#'   via [tidy_mat()].
#'
#' @param fit A fitted TAM object (from [fit_tam()]) containing an `sdrep`
#'   (an [RTMB::sdreport()] object) and `obj$env$.random` (names of random effects).
#' @param interval Confidence level for Wald intervals; default `0.95`.
#'
#' @return
#' A list with two elements:
#' - `fixed`: a data frame stacking all fixed-effect parameters with columns
#'   parameter/index columns (`par`, `coef`, `year`, `age`, as applicable), then
#'   `est`, `lwr`, `upr`, `se`, `se_scale`, and `is_proj` when applicable.
#' - `random`: a named list of data frames (one per random block) with the same
#'   columns as `fixed` (indices appropriate to each random effect).
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' par_tab <- tidy_par(fit)
#' names(par_tab)
#' par_tab$fixed
#' names(par_tab$random)        # e.g., "log_f", "log_r", "missing", ...
#' head(par_tab$random$log_f)
#'
#' @seealso [tidy_mat()], [fit_tam()]
#' @export
tidy_par <- function(fit, interval = 0.95) {
  fit <- .require_tam_fit(fit, arg = "fit")

  est <- .tam_parameter_summary(fit, "Estimate")
  se  <- .tam_parameter_summary(fit, "Std. Error")
  nms <- intersect(names(est), names(se))
  nms <- nms[lengths(est[nms]) > 0L]

  ran_nms <- intersect(fit$obj$env$.random, nms)
  fix_nms <- setdiff(nms, ran_nms)
  z       <- stats::qnorm(0.5 + interval / 2)

  .par2df <- function(nm) {
    e <- est[[nm]]; s <- se[[nm]]
    if (is.matrix(e)) {
      df <- tidy_mat(e, value_name = "est")
      df$is_proj <- df$year %in% fit$dat$years[fit$dat$is_proj]
      df$se <- as.vector(s)
    } else {
      if (is.null(names(e))) {
        df <- data.frame(coef = NA, est = e, se = s)
      } else {
        if (nm == "log_r") {
          df <- data.frame(year = fit$dat$years[-1], est = e, se = s, is_proj = fit$dat$is_proj[-1])
        } else if (nm == "log_n0") {
          df <- data.frame(coef = names(e), age = as.integer(names(e)), est = e, se = s)
        } else {
          df <- data.frame(coef = names(e), est = e, se = s)
        }
      }
    }
    df <- cbind(data.frame(par = nm), df)
    df$lwr <- df$est - z * df$se
    df$upr <- df$est + z * df$se

    df$se_scale <- if (startsWith(nm, "logit_")) {
      "logit"
    } else if (startsWith(nm, "log_")) {
      "log"
    } else {
      "reported"
    }

    if (nm == "logit_q") {
      df <- trans_est(df, transform = NULL, scale = 1)
    } else if (startsWith(nm, "logit_")) {
      df <- trans_est(df, transform = plogis, scale = 1)
      df$par <- sub("^logit_", "", df$par)
    } else if (startsWith(nm, "log_")) {
      df <- trans_est(df, transform = exp, scale = 1)
      df$par <- sub("^log_",   "", df$par)
    } else {
      df <- trans_est(df, transform = NULL, scale = 1)
    }
    values <- c("est", "lwr", "upr", "se", "se_scale")
    df[, c(setdiff(names(df), c(values, "is_proj")), values,
           intersect("is_proj", names(df))), drop = FALSE]
  }

  fixed  <- if (length(fix_nms)) stack_list(lapply(fix_nms, .par2df), label = NULL) else
    data.frame(par = character(), est = numeric(), lwr = numeric(), upr = numeric(),
               se = numeric(), se_scale = character(), check.names = FALSE)
  rownames(fixed) <- NULL
  values <- c("est", "lwr", "upr", "se", "se_scale")
  fixed <- fixed[, c(setdiff(names(fixed), c(values, "is_proj")), values,
                     intersect("is_proj", names(fixed))), drop = FALSE]

  random <- stats::setNames(lapply(ran_nms, .par2df), ran_nms)

  attr(fixed, "interval") <- interval
  attr(random, "interval") <- interval

  list(fixed = fixed, random = random)
}

.tam_parameter_summary <- function(fit, what) {
  if (!is.null(fit[["sdrep"]])) return(as.list(fit[["sdrep"]], what))
  estimates <- fit$obj$env$parList(par = fit$parameter_values)
  attr(estimates, "check.passed") <- NULL
  attr(estimates, "what") <- what
  if (what == "Estimate") return(estimates)
  lapply(estimates, function(x) { x[] <- NA_real_; x })
}


#' Stack a list of tables with an identifier column
#'
#' @description
#' Convenience wrapper around [base::rbind()] for stacking a (named) list of
#' data frames while recording the source list element in a label column. When
#' the list is named, the names are used as labels; otherwise, integer indices
#' (`"1"`, `"2"`, …) are used. Columns absent in some elements are filled with
#' `NA` so inputs need only share column *names* in common.
#'
#' @param x A list of data frames (or objects coercible to data frames).
#' @param label A character scalar giving the identifier column name to add. Set
#'   to `NULL` to omit the identifier column. Default is `"model"`.
#' @param label_type Desired type for the identifier column: automatic
#'   conversion via [utils::type.convert()] (`"auto"`, the default), or
#'   explicit coercion to `"numeric"`, `"character"`, or `"factor"`.
#'
#' @return A single data frame produced by row-binding the list elements. If
#'   `label` is not `NULL`, the identifier column is the first column in the
#'   output.
#'
#' @examples
#' lst <- list(
#'   retro_1 = data.frame(age = 2:4, rho = runif(3)),
#'   retro_2 = data.frame(age = 2:4, rho = runif(3))
#' )
#' stack_list(lst, label = "retro")
#'
#' @export
stack_list <- function(x, label = "model",
                       label_type = c("auto", "numeric", "character", "factor")) {
  label_type <- match.arg(label_type)

  if (!length(x)) {
    cli::cli_abort("{.arg x} must contain at least one element.")
  }

  ids <- names(x)
  if (is.null(ids)) ids <- as.character(seq_along(x))

  pieces <- lapply(seq_along(x), function(i) {
    df <- x[[i]]
    if (is.null(df)) return(NULL)
    if (!is.data.frame(df)) df <- as.data.frame(df)
    if (!nrow(df)) return(NULL)
    if (!is.null(label)) df[[label]] <- ids[[i]]
    df
  })

  keep <- !vapply(pieces, is.null, logical(1))
  if (!any(keep)) {
    cli::cli_abort("No data frames to stack.")
  }
  pieces <- pieces[keep]
  ids <- ids[keep]

  all_cols <- Reduce(union, lapply(pieces, names))
  pieces <- lapply(pieces, function(df) {
    missing <- setdiff(all_cols, names(df))
    if (length(missing)) {
      for (nm in missing) df[[nm]] <- NA
    }
    df <- df[, all_cols, drop = FALSE]
    df
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL

  if (!is.null(label)) {
    out[[label]] <- switch(label_type,
                           auto      = utils::type.convert(out[[label]], as.is = TRUE),
                           numeric   = suppressWarnings(as.numeric(out[[label]])),
                           character = as.character(out[[label]]),
                           factor    = factor(out[[label]], levels = ids)
    )
    ord <- c(label, setdiff(names(out), label))
    out <- out[, ord, drop = FALSE]
  }

  out
}


#' Stack identically named subtables from a nested list
#'
#' @param x A named list of results (e.g. sims or models), each containing
#'   a named list of data.frames (e.g. "ssb", "N", "recruitment", ...).
#'   Shape: list(<id> = list(<subtable> = data.frame, ...), ...)
#' @param label Name of the column to add with the outer id. Default "model".
#'   Set to NULL to omit the id column.
#' @param label_type Controls how the label column is coerced. One of:
#'   - `"auto"` (default): numeric when possible, else character.
#'   - `"numeric"`: force numeric conversion (with `NA` for non-numeric).
#'   - `"character"`: keep as character.
#'   - `"factor"`: convert to factor.
#'
#' @return A named list of data.frames. One element per subtable name. Each
#'   data.frame is the row-bound stack across outer ids, with an added column
#'   named by `label` (unless `label = NULL`).
#' @examples
#' res <- list(
#'   sim1 = list(ssb = data.frame(year=1:3, est=1:3),
#'               N   = data.frame(year=1:2, age=2:3, est=5:6)),
#'   sim2 = list(ssb = data.frame(year=1:3, est=11:13),
#'               N   = data.frame(year=1:2, age=2:3, est=15:16))
#' )
#' stacked <- stack_nested(res, label = "sim")
#' str(stacked$ssb)  # has column 'sim'
#'
#' @importFrom stats setNames
#'
#' @export
stack_nested <- function(x, label = "model",
                         label_type = c("auto", "numeric", "character", "factor")) {
  label_type <- match.arg(label_type)

  outer_ids <- names(x)
  if (is.null(outer_ids)) outer_ids <- as.character(seq_along(x))
  sub_names <- Reduce(union, lapply(x, names))

  out <- stats::setNames(vector("list", length(sub_names)), sub_names)

  for (nm in sub_names) {
    pieces <- lapply(seq_along(x), function(i) {
      if (!is.null(x[[i]][[nm]])) x[[i]][[nm]] else NULL
    })
    names(pieces) <- outer_ids
    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    if (!length(pieces)) {
      out[[nm]] <- NULL
      next
    }
    out[[nm]] <- stack_list(pieces, label = label, label_type = label_type)
  }
  out
}


#' Stack TAM outputs across models (retro folds, scenarios, etc.)
#'
#' @description
#' Builds tidy, stacked tables from one or more fitted TAM models. For each model it collects:
#'
#' - observation diagnostics from the model's `$obs_pred` component (falling back to [tidy_obs_pred()] when absent),
#' - population summaries from `$pop` (falling back to [tidy_pop()] when absent), and
#' - parameter summaries from `$fixed_par`/`$random_par` (falling back to [tidy_par()] when absent),
#'
#' then stacks **per component** across models (e.g., all `"catch"` tables together; all `"N"`
#' tables together; all fixed parameters together; each random-effect block together).
#'
#' @details
#' **Inputs:** Pass models through `...` or via `model_list =`.
#' Precomputed assessment reference objects of class `tam_ref` are also accepted. Their
#' existing tables and uncertainty are retained without constructing a tinyAM
#' optimizer. A `tam_ref` may contain a named `comparison_scales` vector;
#' those factors are applied to the matching tinyAM population estimates and
#' uncertainty so the combined tables use the reference's units.
#'
#' - If `...` supplies **one** model, **no label** column is added.
#' - If `...` supplies **>1** model, a label column is added using the object/expression names from `...`.
#' - If `model_list` is used, it **must be a named list**; its names are always used as labels (even when length 1).
#'
#' Names (from `...` or `model_list`) are passed through [utils::type.convert()] with `as.is = TRUE`,
#' so numeric-like labels (e.g. `"2010"`, `"2011"`) become numeric.
#'
#' Components present in any model are stacked; columns missing from particular models are filled with
#' `NA` values so diagnostics like one-step-ahead residuals survive stacking. Parameter summaries are
#' stacked separately for **fixed** and **random** effects: fixed effects in a single data frame; random
#' effects as a named list of data frames (one per random-effect block, e.g. `"log_f"`, `"log_r"`,
#' `"missing"`, …).
#'
#' @param ... One or more fitted TAM objects (as returned by [fit_tam()]). Supply these or `model_list`, not both.
#' @param model_list A **named list** of fitted TAM objects. Required to be named; the names are used as label values.
#' @param interval Confidence level passed to [tidy_pop()] and [tidy_par()] for interval construction. Default `0.95`.
#' @param label Character scalar giving the label column name to add when stacking across multiple/named models. Default `"model"`.
#' @inheritParams stack_nested
#'
#' @return
#' A named list with four basic elements, and `formula_effects` when present:
#'
#' - **obs_pred** — a named list of stacked data frames (e.g., `catch`, `index`);
#' - **pop** — a named list of stacked data frames (e.g., `ssb`, `N`, `M`, `mu_M`, `F`, `mu_F`, `Z`, …);
#' - **fixed_par** — a single stacked data frame of fixed-effect parameters with columns like `par`, `est`, `se`, `lwr`, `upr`, plus indices (e.g., `coef`, `year`, `age`) and the label column when applicable;
#' - **random_par** — a named list of stacked data frames, one per random-effect block, each with the same schema as `fixed_par` plus block-appropriate indices.
#' - **formula_effects** — signed Gaussian effect levels, RW increments and
#'   numeric-by contributions, grouped by term. Intervals use their joint fitted
#'   uncertainty; the zero RW anchor is explicitly marked as fixed.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' # Single model: no label column added
#' fit1 <- fit
#' tabs1 <- tidy_tam(fit1)
#' names(tabs1)
#' head(tabs1$fixed_par)
#'
#' # Two models via ...: label uses object names
#' fit2 <- update(fit1, years = 1983:2023)
#' tabs2 <- tidy_tam(fit1, fit2)
#' head(tabs2$obs_pred$catch)      # contains column "model"
#' head(tabs2$fixed_par)           # fixed effects stacked with "model"
#' names(tabs2$random_par)         # e.g. "log_f", "log_r", "missing", ...
#' head(tabs2$random_par$log_f)    # random block stacked with "model"
#'
#' # Named list: must be named; names used as labels (even length 1)
#' fits <- list(`2023` = fit2, `2024` = fit1)
#' tabs3 <- tidy_tam(model_list = fits, label = "retro_year")
#' head(tabs3$pop$N)               # column "retro_year" has 2023/2024
#'
#' @importFrom utils type.convert
#' @seealso [tidy_obs_pred()], [tidy_pop()], [tidy_par()], [fit_tam()], [fit_retro()]
#' @export
tidy_tam <- function(..., model_list = NULL, interval = 0.95, label = "model", label_type = "auto") {
  dots <- list(...)
  dot_expr <- as.list(substitute(list(...)))[-1]
  fits_info <- .dots_or_list(dots, dot_expr, model_list = model_list, list_arg_name = "model_list")
  model_list <- fits_info$fits
  using_dots <- fits_info$using_dots

  # add label if: multiple via ... OR any named model_list usage
  add_label <- (using_dots && length(model_list) > 1L) || (!using_dots)
  id_col    <- if (add_label) label else NULL

  cached_interval_values <- unlist(lapply(model_list, function(fit) {
    intervals <- list(
      attr(fit$pop, "interval", exact = TRUE),
      attr(fit$fixed_par, "interval", exact = TRUE),
      attr(fit$random_par, "interval", exact = TRUE)
    )
    intervals <- Filter(Negate(is.null), intervals)
    unlist(intervals, use.names = FALSE)
  }), use.names = FALSE)
  cached_interval_values <- as.numeric(cached_interval_values)
  cached_interval_values <- cached_interval_values[!is.na(cached_interval_values)]
  cached_intervals <- unique(cached_interval_values)

  if (length(model_list) > 1L && length(cached_intervals) > 1L) {
    sorted <- sort(cached_intervals)
    interval_msg <- cli::format_inline("{.val {sorted}}")
    cli::cli_inform(c(
      "i" = "tidy_tam detected cached summaries built with different intervals ({interval_msg}).",
      " " = "Recomputing summaries at interval {.val {interval}} to align across models."
    ))
  }

  interval_matches <- function(x) {
    stored <- attr(x, "interval", exact = TRUE)
    !is.null(stored) && isTRUE(all.equal(stored, interval))
  }
  obs_list <- lapply(model_list, function(fit) {
    if (inherits(fit, "tam_ref")) return(fit$obs_pred)
    if (!is.null(fit$obs_pred) && interval_matches(fit$obs_pred)) return(fit$obs_pred)
    tables <- tidy_obs_pred(fit, interval = interval)
    for (nm in names(tables)) {
      extras <- setdiff(names(fit$obs_pred[[nm]]), names(tables[[nm]]))
      tables[[nm]][extras] <- fit$obs_pred[[nm]][extras]
    }
    tables
  })

  pop_list <- lapply(model_list, function(fit) {
    if (inherits(fit, "tam_ref")) return(fit$pop)
    if (!is.null(fit$pop) && interval_matches(fit$pop)) {
      fit$pop
    } else {
      tidy_pop(fit, interval = interval)
    }
  })
  references <- Filter(function(fit) inherits(fit, "tam_ref"), model_list)
  comparison_scales <- lapply(references, `[[`, "comparison_scales")
  comparison_scales <- Filter(function(x) length(x) > 0L, comparison_scales)
  if (length(comparison_scales)) {
    comparison_scales <- lapply(comparison_scales, function(x) {
      if (is.list(x)) x <- unlist(x, use.names = TRUE)
      if (!is.numeric(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
          anyDuplicated(names(x)) || any(!is.finite(x)) || any(x <= 0)) {
        cli::cli_abort("A {.cls tam_ref} {.field comparison_scales} value must be a named vector of positive finite numbers.")
      }
      x
    })
    if (length(comparison_scales) > 1L &&
        !all(vapply(comparison_scales[-1L], identical, logical(1), comparison_scales[[1L]]))) {
      cli::cli_abort("Assessment references use different {.field comparison_scales}; the fitted outputs cannot be shown on one scale.")
    }
    comparison_scales <- comparison_scales[[1L]]
    reference_pop <- references[[1L]]$pop
    for (i in seq_along(model_list)) {
      if (inherits(model_list[[i]], "tam_ref")) next
      for (metric in intersect(names(comparison_scales), names(pop_list[[i]]))) {
        tab <- pop_list[[i]][[metric]]
        scale <- comparison_scales[[metric]]
        for (field in intersect(c("est", "se", "lwr", "upr"), names(tab))) {
          tab[[field]] <- tab[[field]] * scale
        }
        reference <- reference_pop[[metric]]
        if (is.data.frame(reference) && "unit" %in% names(reference)) {
          units <- unique(as.character(reference$unit[!is.na(reference$unit) &
                                                        nzchar(as.character(reference$unit))]))
          if (length(units) == 1L) tab$unit <- units[[1L]]
        }
        pop_list[[i]][[metric]] <- tab
      }
    }
  }
  par_list <- lapply(model_list, function(fit) {
    if (inherits(fit, "tam_ref")) {
      return(list(fixed = fit$fixed_par, random = fit$random_par))
    }
    has_fixed  <- !is.null(fit$fixed_par)
    has_random <- !is.null(fit$random_par)
    if (has_fixed && has_random && interval_matches(fit$fixed_par) && interval_matches(fit$random_par)) {
      list(fixed = fit$fixed_par, random = fit$random_par)
    } else {
      tidy_par(fit, interval = interval)
    }
  })

  obs_pred <- stack_nested(obs_list, label = id_col, label_type = label_type)
  pop      <- stack_nested(pop_list, label = id_col, label_type = label_type)
  fixed_par  <- stack_nested(lapply(par_list, `[`, "fixed"), label = id_col, label_type = label_type)
  random_par <- stack_nested(lapply(par_list, `[[`, "random"), label = id_col, label_type = label_type)

  fixed_tbl <- fixed_par$fixed
  attr(fixed_tbl, "interval") <- interval

  random_tbls <- lapply(random_par, function(tab) {
    attr(tab, "interval") <- interval
    tab
  })
  attr(random_tbls, "interval") <- interval

  out <- list(
    obs_pred   = obs_pred,
    pop        = pop,
    fixed_par  = fixed_tbl,
    random_par = random_tbls
  )
  effects <- lapply(model_list, function(fit) {
    if (inherits(fit, "tam_ref")) return(list())
    if (!is.null(fit$formula_effects) && interval_matches(fit$formula_effects)) {
      fit$formula_effects
    } else .tidy_formula_effects(fit, interval = interval)
  })
  if (any(lengths(effects))) {
    components <- c("levels", "increments", "contributions")
    out$formula_effects <- stats::setNames(lapply(components, function(nm) {
      stack_nested(lapply(effects, `[[`, nm), label = id_col, label_type = label_type)
    }), components)
  }

  attr(out, "interval") <- interval
  attr(out, "label") <- id_col
  attr(out, "label_type") <- label_type

  out
}

