#' Extract observation residuals
#'
#' Return keyed catch and survey rows for plotting or further checks. This is a
#' convenience extractor: standardized residuals are the existing fitted-state
#' residuals, affected by shrinkage, not independent prediction errors.
#' @param object A fitted `tam_fit`.
#' @param type `"standardized"` (default), `"log"` (log observation minus log
#'   prediction), `"response"` (observation minus prediction), or `"osa"`.
#' @param ... Reserved for future arguments; none are currently used.
#' @details OSA residuals must already have been computed, for example with
#'   `fit$obs_pred <- tidy_obs_pred(fit, add_osa_res = TRUE)`. This method never
#'   refits the model or computes OSA residuals implicitly. Excluded observations,
#'   including zeros, retain rows with missing residuals.
#' @return A data frame with year, age, type (catch/index), survey, observations,
#'   predictions, residuals, and available unit/fleet columns. The `residual_type`
#'   attribute names the requested scale.
#' @examples
#' \dontrun{
#' residuals(fit)
#' residuals(fit, type = "log")
#' }
#' @seealso [tidy_obs_pred()], [plot_bubbles()]
#' @importFrom stats residuals
#' @export
residuals.tam_fit <- function(object, type = c("standardized", "log", "response", "osa"), ...) {
  object <- .require_tam_fit(object, arg = "object")
  type <- match.arg(type)
  tables <- object$obs_pred[c("catch", "index")]
  if (type == "osa" && !any(vapply(tables, function(x) "osa_res" %in% names(x), logical(1)))) {
    cli::cli_abort("OSA residuals have not been computed. Use fit$obs_pred <- tidy_obs_pred(fit, add_osa_res = TRUE) explicitly.")
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

