#' Extract observation residuals
#'
#' Return keyed catch and survey rows for plotting or further checks. This is a
#' convenience extractor: standardized residuals are the existing fitted-state
#' residuals, affected by shrinkage, not independent prediction errors.
#' @param object A fitted `tam_fit`.
#' @param type `"standardized"` (default), `"log"` (log observation minus log
#'   prediction), `"response"` (observation minus prediction), or `"osa"`.
#' @param method OSA algorithm, as in [add_osa_residuals()]. Only used with
#'   `type = "osa"`. Explicitly supplying it requests a fresh calculation.
#' @param ... Computation controls passed to [RTMB::oneStepPredict()], such as
#'   `trace = FALSE` or `parallel = FALSE`.
#'   `subset` and `conditional` index included observations in likelihood order;
#'   unrequested rows retain missing residuals.
#'   Explicit controls request a fresh calculation; other residual types reject them.
#' @details `type = "osa"` reuses cached OSA columns when no computation controls
#' or method are supplied. Otherwise both catch and index residuals are calculated
#' together, using one algorithm at fitted parameter values. It can be expensive,
#' but does not rerun [fit_tam()] or re-estimate fixed parameters. Temporary
#' objectives protect the fitted object's state. Other residual types remain
#' simple extractors. Excluded observations, including zeros, retain rows with
#' missing residuals. The supplied fit is not modified or implicitly cached.
#' To retain results for later extraction and dashboard use, assign
#' `fit_osa <- add_osa_residuals(fit, method = "oneStepGeneric", trace = FALSE)`.
#' @return A data frame with year, age, type (catch/index), survey, observations,
#'   predictions, residuals, and available unit/fleet columns. OSA output also
#'   identifies `osa_method`. The `residual_type` attribute names the requested
#'   scale, and `osa` retains available method/settings/coverage metadata.
#' @examples
#' \dontrun{
#' residuals(fit)
#' residuals(fit, type = "log")
#' residuals(fit, type = "osa", method = "oneStepGeneric", trace = FALSE)
#' }
#' @seealso [tidy_obs_pred()], [plot_bubbles()]
#' @importFrom stats residuals
#' @export
residuals.tam_fit <- function(object, type = c("standardized", "log", "response", "osa"),
                             method = "oneStepGaussianOffMode", ...) {
  object <- .require_tam_fit(object, arg = "object")
  type <- match.arg(type)
  controls <- list(...)
  explicit <- !missing(method) || length(controls) > 0L
  .osa_method(method, controls)
  if (type != "osa" && explicit) cli::cli_abort("Additional computation controls are only used with type = 'osa'.")
  tables <- .observation_tables(object)
  metadata <- attr(tables, "osa", exact = TRUE)
  if (is.null(metadata)) metadata <- object$osa
  missing_osa <- vapply(tables, function(d) is.data.frame(d) && nrow(d) > 0L &&
    !"osa_res" %in% names(d) && any(is.finite(d$obs) & d$obs > 0 & if (is.null(d$is_proj)) TRUE else !d$is_proj), logical(1))
  if (type == "osa" && (explicit || any(missing_osa))) {
    computed <- .compute_osa_residuals(object, method, controls)
    metadata <- computed$metadata
    tables <- .attach_osa_tables(tables, computed$tables, metadata)
  }

  out <- .bind_rows(lapply(names(tables), function(name) {
    d <- tables[[name]]
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    residual <- switch(type, standardized = d$std_res, osa = d$osa_res,
      response = d$obs - d$pred, log = {
        value <- suppressWarnings(log(d$obs) - log(d$pred))
        if (!is.null(d$sd) && !is.null(d$std_res)) {
          underflow <- which(d$pred == 0 & d$obs > 0 & is.finite(d$std_res) & is.finite(d$sd))
          value[underflow] <- d$std_res[underflow] * d$sd[underflow]
        }
        value
      })
    if (is.null(residual)) residual <- rep(NA_real_, nrow(d))
    keep <- is.finite(d$obs) & d$obs > 0 & is.finite(residual)
    if ("is_proj" %in% names(d)) keep <- keep & !d$is_proj
    residual[!keep] <- NA_real_
    if (type == "osa" && !"osa_method" %in% names(d)) {
      d$osa_method <- rep(if (is.null(metadata$method)) NA_character_ else metadata$method, nrow(d))
    }
    columns <- c("year", "age", "survey", "fleet", "obs", "pred", "unit", "is_proj")
    if (type == "osa") columns <- c(columns, "osa_method")
    d <- d[intersect(columns, names(d))]
    d$type <- name; d$residual <- residual
    if (!"survey" %in% names(d)) d$survey <- NA_character_
    d
  }))
  attr(out, "residual_type") <- type
  if (type == "osa" && !is.null(metadata)) attr(out, "osa") <- metadata
  out
}

