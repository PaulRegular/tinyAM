# Ordinary proportional cohort gains/losses, evaluated at destination year/age.
.parse_cohort_mean <- function(dat) {
  formula <- dat$N_settings$mu_form
  if (is.null(formula)) return(matrix(numeric(), nrow(dat$pop_data), 0L))
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    cli::cli_abort("N_settings$mu_form must be a one-sided formula or NULL.")
  }
  forbidden <- function(x) {
    if (!is.call(x)) return(FALSE)
    .formula_call_name(x) %in% c("iid", "rw", "ar1", "mono", "logistic", "bh", "ricker", "|", "offset") ||
      any(vapply(as.list(x)[-1L], forbidden, logical(1)))
  }
  if (forbidden(formula[[2L]])) {
    cli::cli_abort("N_settings$mu_form supports ordinary fixed covariates only; specify residual dependence with N_settings$process.")
  }
  d <- .population_data(dat, formula, "N_settings$mu_form")
  active <- d$year > min(dat$years) & d$age > min(dat$ages)
  historical <- active & !dat$is_proj[match(d$year, dat$years)]
  frame <- stats::model.frame(formula, d[historical, , drop = FALSE], na.action = stats::na.pass)
  design <- stats::model.matrix(attr(frame, "terms"), frame)
  if (nrow(design) != sum(historical) || any(!is.finite(design))) {
    cli::cli_abort("N_settings$mu_form needs finite covariates at every historical destination year and older modeled age in data$pop_covariates.")
  }
  # Reuse training terms (including poly() bases and factor contrasts) in projections.
  projected <- stats::model.frame(attr(frame, "terms"), d[active, , drop = FALSE], na.action = stats::na.pass,
    xlev = stats::.getXlevels(attr(frame, "terms"), frame))
  full <- stats::model.matrix(attr(frame, "terms"), projected, contrasts.arg = attr(design, "contrasts"))
  if (nrow(full) != sum(active) || any(!is.finite(full))) {
    cli::cli_abort("N_settings$mu_form needs finite destination-year covariates, including projections.")
  }
  informative <- design
  if (dat$N_settings$process == "rw" && ncol(design)) {
    ny <- sum(!dat$is_proj) - 1L
    informative <- vapply(seq_len(ncol(design)), function(j) {
      z <- matrix(design[, j], ny, length(dat$ages) - 1L)
      as.vector(z[-1L, , drop = FALSE] - z[-ny, , drop = FALSE])
    }, numeric((ny - 1L) * (length(dat$ages) - 1L)))
  }
  if (!ncol(design) || nrow(informative) < ncol(design) || qr(informative)$rank < ncol(design)) {
    cli::cli_abort(c("N_settings$mu_form has redundant or unsupported mean directions.",
      i = if (dat$N_settings$process == "rw") "Constant mean directions cancel under an intrinsic RW. Use ~ 0 + time-varying covariates with independent temporal contrasts." else "Remove redundant columns and provide variation over destination years and ages."))
  }
  out <- matrix(0, nrow(d), ncol(full), dimnames = list(NULL, colnames(full)))
  out[active, ] <- full
  out
}
