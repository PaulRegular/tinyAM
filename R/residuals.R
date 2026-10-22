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

# Use the chronological population constructor rather than duplicating survival
# equations. No fitted objective or parameter mode is changed.
.tidy_process_deviations <- function(fit) {
  if (is.null(fit$rep$Z) || (is.null(fit$parameter_values) && is.null(fit$sdrep))) return(list())
  dat <- fit$dat; par <- .start_parameters(fit, dat)
  out <- list()
  table <- function(x, process, block = FALSE) {
    if (process %in% c("rw", "cor_rw")) {
      if (nrow(x) < 2L) return(NULL)
      x <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
    }
    names(dimnames(x)) <- c("year", if (block) "age_block" else "age")
    d <- tidy_array(x, value_name = "est")
    d$is_proj <- d$year %in% dat$years[dat$is_proj]
    d$kind <- if (process %in% c("rw", "cor_rw")) "Annual log increment" else "Log deviation from mean"
    d
  }
  if (!is.null(par[["log_f"]])) {
    historical <- !dat$is_proj
    out$F <- table(log(fit$rep$F[historical, , drop = FALSE]) -
      log(fit$rep$mu_F[historical, , drop = FALSE]), dat$F_settings$process)
  }
  if (!is.null(par[["log_m"]])) {
    out$M <- table(par[["log_m"]] - log(fit$rep$mu_M[rownames(par[["log_m"]]),
      dat$M_settings$age_block_start, drop = FALSE]), dat$M_settings$process, TRUE)
  }
  if (!is.null(par[["log_n"]])) {
    population <- .population_states(par, dat, fit$rep$Z, F = fit$rep$F, M = fit$rep$M)
    out$N <- table(log(fit$rep$N[-1L, -1L, drop = FALSE]) -
      population$pred_log_N[-1L, -1L, drop = FALSE], dat$N_settings$process)
    if (!is.null(out$N) && dat$N_settings$process != "rw") out$N$kind <- "Log deviation from cohort prediction"
  }
  Filter(function(x) is.data.frame(x) && nrow(x), out)
}
