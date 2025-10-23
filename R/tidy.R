
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
#' @param add_osa_res Logical; add one-step-ahead residuals? Hard wired to
#'                    apply the `"oneStepGaussianOffMode"` method.
#'                    See [RTMB::oneStepPredict()] for details.
#' @param ... Arguments to pass to [RTMB::oneStepPredict()].
#'
#' @return
#' A named list with two data frames:
#'
#' - **catch**: original columns plus `pred`, `sd`, and `std_res`.
#' - **index**: original columns plus `pred`, `sd`, `q`, and `std_res`.
#'
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' obs_pred <- tidy_obs_pred(fit)
#' head(obs_pred$catch)
#' head(obs_pred$index)
#'
#' @seealso [fit_tam()], [tidy_rep()], [tidy_sdrep()], [tidy_pop()]
#' @export
tidy_obs_pred <- function(fit, add_osa_res = FALSE, ...) {
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

  if (add_osa_res) {
    osa_res <- RTMB::oneStepPredict(fit$obj, method = "oneStepGaussianOffMode", ...)
    split_osa_res <- split(osa_res$residual, fit$dat$obs_map$type[fit$dat$obs_map$is_observed])
    split_is_observed <- split(fit$dat$obs_map$is_observed, fit$dat$obs_map$type)
    obs_pred$catch$osa_res <- obs_pred$index$osa_res <- NA
    obs_pred$catch$osa_res[split_is_observed$catch] <- split_osa_res$catch
    obs_pred$index$osa_res[split_is_observed$index] <- split_osa_res$index
  }

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
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' trends <- tidy_rep(fit)
#' names(trends)
#' head(trends$N)
#'
#' @seealso [tidy_obs_pred()], [tidy_sdrep()], [tidy_pop()], [tidy_mat()]
#' @export
tidy_rep <- function(fit) {
  which_mat <- which(sapply(fit$rep, is.matrix))
  which_vec <- which(sapply(fit$rep, length) == length(fit$dat$years))
  nms <- c(names(which_mat), names(which_vec))
  trends <- lapply(nms, function(nm) {
    if (is.matrix(fit$rep[[nm]])) {
      d <- tidy_mat(fit$rep[[nm]], value_name = "est")
      d$is_proj <- d$year %in% fit$dat$years[fit$dat$is_proj]
    } else {
      d <- data.frame(year = fit$dat$years, est = fit$rep[[nm]], is_proj = fit$dat$is_proj)
    }
    d
  })
  names(trends) <- nms
  trends
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
    data[, c("est", "lwr", "upr")] <-  apply(data[, c("est", "lwr", "upr")], 2, transform)
  }
  data[, c("est", "lwr", "upr")] <- data[, c("est", "lwr", "upr")] / scale
  data
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
#' - `year`, `est`, `sd`, `lwr`, `upr`, `is_proj` — after applying the chosen transform.
#'
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' trends <- tidy_sdrep(fit, interval = 0.9)
#' names(trends)
#' head(trends$ssb)
#'
#' @seealso [trans_est()], [tidy_rep()], [tidy_pop()]
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
                    upr = vals[[i]] + qnorm(1 - ((1 - interval) / 2)) * sds[[i]],
                    is_proj = fit$dat$is_proj) |>
      trans_est(transform = exp)
  })
  names(df) <- gsub("log_", "", names(vals))
  df
}

#' Collect population summaries (sdreport + report trends)
#'
#' @description
#' Convenience wrapper that combines [tidy_sdrep()] and [tidy_rep()]
#' into a single named list for downstream plotting and summaries.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param interval Confidence level for intervals passed to [tidy_sdrep()]; default `0.95`.
#'
#' @return
#' A named list containing the elements returned by [tidy_sdrep()] and
#' [tidy_rep()] (names preserved).
#'
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' pop <- tidy_pop(fit)
#' names(pop)
#'
#' @seealso [tidy_sdrep()], [tidy_rep()], [tidy_obs_pred()]
#' @export
tidy_pop <- function(fit, interval = 0.95) {
  sdrep_trends <- tidy_sdrep(fit, interval = interval)
  rep_trends <- tidy_rep(fit)
  not_in_sdrep <- setdiff(names(rep_trends), names(sdrep_trends))
  c(sdrep_trends, rep_trends[not_in_sdrep])
}


#' Tidy parameter estimates (fixed & random) with CIs and back-transforms
#'
#' @description
#' Creates a tidy summary of parameter estimates from a fitted TAM object,
#' combining estimates (`Estimate`) and standard errors (`Std. Error`) from
#' `fit$sdrep`. Parameters whose names begin with `log_` or `logit_` are
#' back-transformed to the natural scale:
#'
#' - `log_`  → `exp()` (and the `log_` prefix is dropped, e.g. `log_sd_r` → `sd_r`)
#' - `logit_` → `plogis()` (and the `logit_` prefix is dropped, e.g. `logit_phi_f` → `phi_f`)
#'
#' Fixed-effect parameters are returned in a single data frame (`$fixed`);
#' random-effect parameters are returned as a named list of data frames
#' (`$random`), one per random block (e.g. `log_f`, `log_r`, `missing`, …).
#'
#' Labels are added where applicable:
#' - For parameters specified using a formula in [make_dat()] (e.g., `log_q`, `log_sd_obs`),
#;   a `coef` column is added.
#' - For `log_r` (recruitment path), a `year` column is used.
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
#'   `par`, `est`, `se`, `lwr`, `upr`, plus any index columns such as
#'   `coef`, `year`, `age`.
#' - `random`: a named list of data frames (one per random block) with the same
#'   columns as `fixed` (indices appropriate to each random effect).
#'
#' @examples
#' fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
#' par_tab <- tidy_par(fit)
#' names(par_tab)
#' par_tab$fixed
#' names(par_tab$random)        # e.g., "log_f", "log_r", "missing", ...
#' head(par_tab$random$log_f)
#'
#' @seealso [tidy_mat()], [fit_tam()]
#' @export
tidy_par <- function(fit, interval = 0.95) {

  est <- as.list(fit$sdrep, "Estimate")
  se  <- as.list(fit$sdrep, "Std. Error")
  nms <- intersect(names(est), names(se))

  ran_nms <- fit$obj$env$.random
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
          df <- data.frame(year = fit$dat$years, est = e, se = s, is_proj = fit$dat$is_proj)
        } else {
          df <- data.frame(coef = names(e), est = e, se = s)
        }
      }
    }
    df <- cbind(data.frame(par = nm), df)
    df$lwr <- df$est - z * df$se
    df$upr <- df$est + z * df$se

    if (startsWith(nm, "logit_")) {
      df <- trans_est(df, transform = plogis, scale = 1)
      df$par <- sub("^logit_", "", df$par)
    } else if (startsWith(nm, "log_")) {
      df <- trans_est(df, transform = exp, scale = 1)
      df$par <- sub("^log_",   "", df$par)
    } else {
      df <- trans_est(df, transform = NULL, scale = 1)
    }
    df
  }

  fixed  <- if (length(fix_nms)) do.call(rbind, lapply(fix_nms, .par2df)) else
    data.frame(par = character(), est = numeric(), se = numeric(),
               lwr = numeric(), upr = numeric(), check.names = FALSE)
  rownames(fixed) <- NULL

  random <- stats::setNames(lapply(ran_nms, .par2df), ran_nms)

  list(fixed = fixed, random = random)
}


#' Stack a list of tables with an identifier column
#'
#' @description
#' Convenience wrapper around [base::rbind()] for stacking a (named) list of
#' data frames while recording the source list element in a label column. When
#' the list is named, the names are used as labels; otherwise, integer indices
#' (`"1"`, `"2"`, …) are used.
#'
#' @param x A list of data frames (or objects coercible to data frames) with the
#'   same column structure.
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
    stop("`x` must contain at least one element.", call. = FALSE)
  }

  ids <- names(x)
  if (is.null(ids)) ids <- as.character(seq_along(x))

  pieces <- lapply(seq_along(x), function(i) {
    df <- x[[i]]
    if (!is.data.frame(df)) df <- as.data.frame(df)
    if (!is.null(label)) df[[label]] <- ids[[i]]
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
#'   data.frame is the row-bound stack across outer ids, with an added `id_col`
#'   (if not NULL).
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

  id_col <- label
  outer_ids <- names(x)
  if (is.null(outer_ids)) outer_ids <- as.character(seq_along(x))
  sub_names <- Reduce(union, lapply(x, names))

  out <- stats::setNames(vector("list", length(sub_names)), sub_names)

  for (nm in sub_names) {
    pieces <- lapply(seq_along(x), function(i) {
      if (nm %in% names(x[[i]])) {
        df <- x[[i]][[nm]]
        if (!is.data.frame(df)) df <- as.data.frame(df)
        if (!is.null(id_col)) df[[id_col]] <- outer_ids[i]
        df
      }
    })
    stk <- do.call(rbind, pieces)
    rownames(stk) <- NULL

    if (!is.null(id_col)) {
      # Apply label_type coercion
      stk[[id_col]] <- switch(label_type,
                              auto      = utils::type.convert(stk[[id_col]], as.is = TRUE),
                              numeric   = suppressWarnings(as.numeric(stk[[id_col]])),
                              character = as.character(stk[[id_col]]),
                              factor    = factor(stk[[id_col]], levels = names(x))
      )
      ord <- c(id_col, setdiff(names(stk), id_col))
      stk <- stk[, ord, drop = FALSE]
    }
    out[[nm]] <- stk
  }
  out
}


#' Stack TAM outputs across models (retro folds, scenarios, etc.)
#'
#' @description
#' Builds tidy, stacked tables from one or more fitted TAM models. For each model it collects:
#'
#' - observation diagnostics via [tidy_obs_pred()],
#' - population summaries via [tidy_pop()] (with confidence intervals), and
#' - parameter summaries via [tidy_par()] (fixed and random effects),
#'
#' then stacks **per component** across models (e.g., all `"catch"` tables together; all `"N"`
#' tables together; all fixed parameters together; each random-effect block together).
#'
#' @details
#' **Inputs:** Pass models through `...` or via `model_list =`.
#'
#' - If `...` supplies **one** model, **no label** column is added.
#' - If `...` supplies **>1** model, a label column is added using the object/expression names from `...`.
#' - If `model_list` is used, it **must be a named list**; its names are always used as labels (even when length 1).
#'
#' Names (from `...` or `model_list`) are passed through [utils::type.convert()] with `as.is = TRUE`,
#' so numeric-like labels (e.g. `"2010"`, `"2011"`) become numeric.
#'
#' Only components present in **all** models are stacked (intersection of names), ensuring column
#' compatibility for base `rbind()`. Parameter summaries are stacked separately for **fixed** and
#' **random** effects: fixed effects in a single data frame; random effects as a named list of data
#' frames (one per random-effect block, e.g. `"log_f"`, `"log_r"`, `"missing"`, …).
#'
#' @param ... One or more fitted TAM objects (as returned by [fit_tam()]). Ignored if `model_list` is provided.
#' @param model_list A **named list** of fitted TAM objects. Required to be named; the names are used as label values.
#' @param interval Confidence level passed to [tidy_pop()] and [tidy_par()] for interval construction. Default `0.95`.
#' @param label Character scalar giving the label column name to add when stacking across multiple/named models. Default `"model"`.
#' @inheritParams stack_nested
#'
#' @return
#' A named list with four elements:
#'
#' - **obs_pred** — a named list of stacked data frames (e.g., `catch`, `index`);
#' - **pop** — a named list of stacked data frames (e.g., `ssb`, `N`, `M`, `mu_M`, `F`, `mu_F`, `Z`, …);
#' - **fixed_par** — a single stacked data frame of fixed-effect parameters with columns like `par`, `est`, `se`, `lwr`, `upr`, plus indices (e.g., `coef`, `year`, `age`) and the label column when applicable;
#' - **random_par** — a named list of stacked data frames, one per random-effect block, each with the same schema as `fixed_par` plus block-appropriate indices.
#'
#' @examples
#' # Single model: no label column added
#' fit1 <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14)
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
  using_dots <- is.null(model_list)

  if (using_dots) {
    dots <- list(...)
    if (!length(dots)) stop("No models supplied.", call. = FALSE)
    # name from expressions in ...
    names(dots) <- sapply(substitute(list(...))[-1], deparse1)
    model_list <- dots
  } else {
    if (!length(model_list)) stop("No models supplied.", call. = FALSE)
    if (is.null(names(model_list)) || any(!nzchar(names(model_list)))) {
      stop("`model_list` must be a named list (all names non-empty).", call. = FALSE)
    }
  }

  # add label if: multiple via ... OR any named model_list usage
  add_label <- (using_dots && length(model_list) > 1L) || (!using_dots)
  id_col    <- if (add_label) label else NULL

  obs_list <- lapply(model_list, tidy_obs_pred)
  pop_list <- lapply(model_list, tidy_pop, interval = interval)
  par_list <- lapply(model_list, tidy_par, interval = interval)

  obs_pred <- stack_nested(obs_list, label = id_col, label_type = label_type)
  pop      <- stack_nested(pop_list, label = id_col, label_type = label_type)
  fixed_par  <- stack_nested(lapply(par_list, `[`, "fixed"), label = id_col, label_type = label_type)
  random_par <- stack_nested(lapply(par_list, `[[`, "random"), label = id_col, label_type = label_type)

  list(
    obs_pred   = obs_pred,
    pop        = pop,
    fixed_par  = fixed_par$fixed,
    random_par = random_par
  )
}


