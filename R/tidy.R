
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
    stop("`x` must be a matrix or array.", call. = FALSE)
  }

  nm <- dimnames(x)

  if (require_dimnames) {
    if (is.null(nm) || any(vapply(nm, is.null, logical(1)))) {
      stop("All dimensions must have names (dimnames). Set `require_dimnames = FALSE` to allow defaults.",
           call. = FALSE)
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
#' Extracts observations, predictions, and standard errors for `catch` and `index`
#' from a fitted TAM object, and adds standardized residuals on the log scale.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#'
#' @return
#' A named list with two data frames:
#'
#' - **catch**: original columns plus `pred`, `sd`, and `std_res`.
#' - **index**: original columns plus `pred`, `sd`, `q`, and `std_res`.
#'
#' @examples
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' obs_pred <- tidy_obs_pred(fit)
#' head(obs_pred$catch)
#' head(obs_pred$index)
#'
#' @seealso [fit_tam()], [tidy_rep_mats()], [tidy_sdrep()], [tidy_pop()]
#' @export
tidy_obs_pred <- function(fit) {
  obs_pred <- fit$dat$obs[c("catch", "index")]
  pred <- split(exp(fit$rep$log_pred), fit$dat$obs_map$type)
  sd <- split(fit$rep$sd_obs, fit$dat$obs_map$type)

  obs_pred$catch$pred <- pred$catch
  obs_pred$catch$sd <- sd$catch
  obs_pred$catch$std_res <- with(obs_pred$catch, ifelse(obs == 0, NA, (log(obs) - log(pred)) / sd))

  obs_pred$index$pred <- pred$index
  obs_pred$index$sd <- sd$index
  obs_pred$index$q <- exp(fit$rep$log_q_obs)
  obs_pred$index$std_res <- with(obs_pred$index, ifelse(obs == 0, NA, (log(obs) - log(pred)) / sd))

  obs_pred
}


#' Tidy key report matrices
#'
#' @description
#' Converts selected matrices in `fit$rep` to long (tidy) data frames
#' using [tidy_mat()].
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#'
#' @return
#' A named list of data frames (elements: `N`, `M`, `mu_M`, `F`, `mu_F`, `Z`, `ssb_mat`),
#' where each data frame has one column per dimension (e.g., `year`, `age`) and
#' a value column `est`.
#'
#' @examples
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' mats <- tidy_rep_mats(fit)
#' names(mats)
#' head(mats$N)
#'
#' @seealso [tidy_obs_pred()], [tidy_sdrep()], [tidy_pop()], [tidy_mat()]
#' @export
tidy_rep_mats <- function(fit) {
  mat_nms <- c("N", "M", "mu_M", "F", "mu_F", "Z", "ssb_mat")
  rep_mats <- lapply(mat_nms, function(nm) {
    tidy_mat(fit$rep[[nm]], value_name = "est")
  })
  names(rep_mats) <- mat_nms
  rep_mats
}

#' Tidy `sdreport` time series with confidence intervals
#'
#' @description
#' Extracts ADREPORTED time-series from `fit$sdrep`, computes normal-approximation
#' intervals, applies a transformation (default `exp`), and returns a list of
#' tidy data frames.
#'
#' @details
#' Assumptions:
#'
#' - All ADREPORTED series used here have length equal to `length(fit$dat$years)`.
#' - Estimates are on the log scale and are transformed with `exp` via [trans_est()].
#' - List element names are cleaned by removing a leading `"log_"` prefix.
#'
#' The interval is constructed as `est ± z * sd` with
#' `z = qnorm(1 - (1 - interval) / 2)`.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level in `(0, 1)`; default `0.95`.
#'
#' @return
#' A named list of data frames (one per series), each with columns:
#'
#' - `year`, `est`, `sd`, `lwr`, `upr` — after applying the chosen transform.
#'
#' @examples
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' trends <- tidy_sdrep(fit, interval = 0.9)
#' names(trends)
#' head(trends$ssb)
#'
#' @seealso [trans_est()], [tidy_rep_mats()], [tidy_pop()]
#' @export
tidy_sdrep <- function(fit, interval = 0.95) {
  ## assumes all ADREPORTED objects are equal length to years and are in log space
  vals <- as.list(fit$sdrep, "Estimate", report = TRUE)
  sds <- as.list(fit$sdrep, "Std. Error", report = TRUE)
  df <- lapply(seq_along(vals), function(i) {
    d <- data.frame(year = fit$dat$years,
                    est = vals[[i]],
                    sd = sds[[i]],
                    lwr = vals[[i]] - qnorm(1 - ((1 - interval) / 2)) * sds[[i]],
                    upr = vals[[i]] + qnorm(1 - ((1 - interval) / 2)) * sds[[i]]) |>
      trans_est(transform = exp)
  })
  names(df) <- gsub("log_", "", names(vals))
  df
}

#' Collect population summaries (sdreport + report matrices)
#'
#' @description
#' Convenience wrapper that combines [tidy_sdrep()] and [tidy_rep_mats()]
#' into a single named list for downstream plotting and summaries.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level for intervals passed to [tidy_sdrep()]; default `0.95`.
#'
#' @return
#' A named list containing the elements returned by [tidy_sdrep()] and
#' [tidy_rep_mats()] (names preserved).
#'
#' @examples
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' pop <- tidy_pop(fit)
#' names(pop)
#'
#' @seealso [tidy_sdrep()], [tidy_rep_mats()], [tidy_obs_pred()]
#' @export
tidy_pop <- function(fit, interval = 0.95) {
  c(tidy_sdrep(fit, interval = interval),
    tidy_rep_mats(fit))
}


#' Stack TAM outputs across models (retro folds, scenarios, etc.)
#'
#' @description
#' Builds tidy, stacked tables from one or more fitted TAM models. For each
#' model, it collects:
#'
#' - observation diagnostics via [tidy_obs_pred()], and
#' - population summaries via [tidy_pop()] (with confidence intervals),
#'
#' then stacks **per component** across models (e.g., all `"catch"` tables
#' together; all `"N"` tables together).
#'
#' @details
#' **Inputs:** Pass models through `...` or via `model_list =`.
#'
#' - If `...` supplies **one** model, **no label** column is added.
#' - If `...` supplies **more than one** model, a label column is added using
#'   the object/expression names from `...`.
#' - If `model_list` is used, it **must be a named list**; its names are always
#'   used as labels (even when length 1).
#'
#' Only components present in **all** models are stacked (intersection of
#' component names), ensuring consistent columns for base-R `rbind()`.
#'
#' @param ... One or more fitted TAM objects (as returned by [fit_tam()]).
#'   Ignored if `model_list` is provided.
#' @param model_list A **named list** of fitted TAM objects. Required to be named;
#'   the names are used as label values (even for length-1 lists).
#' @param interval Confidence level passed to [tidy_pop()] for interval construction.
#'   Default `0.95`.
#' @param label Character scalar giving the label column name to add when stacking
#'   across multiple/named models. Default `"model"`.
#'
#' @return
#' A named list with two elements:
#'
#' - **obs_pred** — a named list of stacked data frames (e.g., `catch`, `index`);
#' - **pop** — a named list of stacked data frames (e.g., `ssb`, `N`, `M`,
#'   `mu_M`, `F`, `mu_F`, `Z`, `ssb_mat`).
#'
#' Each stacked data frame contains all rows across models, with a label column
#' added according to the rules above.
#'
#' @examples
#' # Single model via ...: no label column added
#' fit1 <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' tabs1 <- tidy_tam(fit1)
#' head(tabs1$obs_pred$index)
#'
#' # Multiple via ...: label uses the object names
#' fit2 <- update(fit1, years = 1983:2023)
#' tabs2 <- tidy_tam(fit1, fit2)
#' head(tabs2$obs_pred$catch)  # contains column "model"
#'
#' # Named list: must be named; names always used as label (even length 1)
#' fits <- list(`2023` = fit2, `2024` = fit1)
#' tabs3 <- tidy_tam(model_list = fits, label = "retro_year")
#' head(tabs3$pop$N)  # column "retro_year" has 2023/2024
#'
#' @seealso [tidy_obs_pred()], [tidy_pop()], [fit_tam()], [fit_retro()]
#' @export
tidy_tam <- function(..., model_list = NULL, interval = 0.95, label = "model") {
  using_dots <- is.null(model_list)

  if (using_dots) {
    dots <- list(...)
    if (!length(dots)) stop("No models supplied.", call. = FALSE)
    # Simple labels from the ... expressions
    names(dots) <- sapply(substitute(list(...))[-1], deparse1)
    model_list <- dots
  } else {
    if (!length(model_list)) stop("No models supplied.", call. = FALSE)
    if (is.null(names(model_list)) || any(!nzchar(names(model_list)))) {
      stop("`model_list` must be a named list (all names non-empty).", call. = FALSE)
    }
  }

  # Add label iff multiple via ... OR any use of model_list (which must be named)
  add_label <- (using_dots && length(model_list) > 1L) || (!using_dots)

  .common <- function(x) Reduce(intersect, lapply(x, names))
  .stack  <- function(models, comp) {
    parts <- lapply(models, `[[`, comp)
    out <- do.call(rbind, unname(parts))
    if (add_label) out[[label]] <- rep(names(models), sapply(parts, nrow))
    rownames(out) <- NULL
    out
  }

  obs_list <- lapply(model_list, tidy_obs_pred)
  pop_list <- lapply(model_list, tidy_pop, interval = interval)

  obs_nms <- .common(obs_list)
  pop_nms <- .common(pop_list)

  obs_pred <- setNames(lapply(obs_nms, function(nm) .stack(obs_list, nm)), obs_nms)
  pop      <- setNames(lapply(pop_nms, function(nm) .stack(pop_list, nm)), pop_nms)

  list(obs_pred = obs_pred, pop = pop)
}




#' Transform and rescale estimate columns
#'
#' @description
#' Applies a transformation (e.g., `exp`) and optional rescaling to the columns
#' `est`, `lwr`, and `upr` of a data frame.
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
    data[, c("est", "lwr", "upr")] <-  transform(data[, c("est", "lwr", "upr")])
  }
  data[, c("est", "lwr", "upr")] <- data[, c("est", "lwr", "upr")] / scale
  data
}
