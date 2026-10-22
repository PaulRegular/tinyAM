#' Extract observation residuals
#'
#' Return keyed catch and survey rows for plotting or further checks. This is a
#' convenience extractor: standardized residuals are the existing fitted-state
#' residuals, affected by shrinkage, not independent prediction errors.
#' @param object A fitted `tam_fit`.
#' @param type `"standardized"` (default), `"log"` (log observation minus log
#'   prediction), `"response"` (observation minus prediction), or `"osa"`.
#' @param ... Computation controls passed to [tidy_obs_pred()] and
#'   [RTMB::oneStepPredict()] when uncached OSA residuals are requested, such as
#'   `trace = FALSE` or `parallel = FALSE`. Not used for cached OSA values;
#'   other residual types reject these arguments.
#' @details `type = "osa"` reuses cached OSA columns, or explicitly requests
#' their calculation at fitted parameter values when absent. It can be expensive,
#' but does not rerun [fit_tam()] or re-estimate fixed parameters. Temporary
#' objectives protect the fitted object's state. Other residual types remain
#' simple extractors. Excluded observations, including zeros, retain rows with
#' missing residuals. The supplied fit is not modified or implicitly cached.
#' To retain OSA results for later extraction and dashboard use, assign
#' `fit$obs_pred <- tidy_obs_pred(fit, add_osa_res = TRUE, trace = FALSE)`.
#' @return A data frame with year, age, type (catch/index), survey, observations,
#'   predictions, residuals, and available unit/fleet columns. The `residual_type`
#'   attribute names the requested scale.
#' @examples
#' \dontrun{
#' residuals(fit)
#' residuals(fit, type = "log")
#' residuals(fit, type = "osa", trace = FALSE)
#' }
#' @seealso [tidy_obs_pred()], [plot_bubbles()]
#' @importFrom stats residuals
#' @export
residuals.tam_fit <- function(object, type = c("standardized", "log", "response", "osa"), ...) {
  object <- .require_tam_fit(object, arg = "object")
  type <- match.arg(type)
  tables <- object$obs_pred[c("catch", "index")]
  if (type != "osa" && length(list(...))) cli::cli_abort("Additional computation controls are only used with type = 'osa'.")
  missing_osa <- vapply(tables, function(d) is.data.frame(d) && nrow(d) > 0L &&
    !"osa_res" %in% names(d) && any(is.finite(d$obs) & d$obs > 0 & if (is.null(d$is_proj)) TRUE else !d$is_proj), logical(1))
  if (type == "osa" && any(missing_osa)) {
    computed <- tryCatch(tidy_obs_pred(object, add_osa_res = TRUE, ...),
      error = function(e) cli::cli_abort("Could not calculate OSA residuals: {conditionMessage(e)}", parent = e))
    for (name in names(tables)[missing_osa]) {
      d <- tables[[name]]; full <- computed[[name]]
      keys <- intersect(c("year", "age", "survey", "fleet"), names(d))
      key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
      i <- match(key(d), key(full))
      if (anyNA(i)) cli::cli_abort("OSA observation rows do not match the fitted inputs.")
      tables[[name]]$osa_res <- full$osa_res[i]
    }
  }

  out <- .bind_rows(lapply(names(tables), function(name) {
    d <- tables[[name]]
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    residual <- switch(type, standardized = d$std_res, osa = d$osa_res,
      response = d$obs - d$pred, log = suppressWarnings(log(d$obs) - log(d$pred)))
    if (is.null(residual)) residual <- rep(NA_real_, nrow(d))
    keep <- is.finite(d$obs) & d$obs > 0 & is.finite(residual)
    if ("is_proj" %in% names(d)) keep <- keep & !d$is_proj
    residual[!keep] <- NA_real_
    d <- d[intersect(c("year", "age", "survey", "fleet", "obs", "pred", "unit", "is_proj"), names(d))]
    d$type <- name; d$residual <- residual
    if (!"survey" %in% names(d)) d$survey <- NA_character_
    d
  }))
  attr(out, "residual_type") <- type
  out
}

