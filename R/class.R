#' TAM object constructors and validators
#'
#' @keywords internal
#' @noRd
.new_tam_fit <- function(x) {
  required <- c(
    "call", "dat", "obj", "opt", "rep", "sdrep",
    "fixed_par", "random_par", "obs_pred", "pop", "is_converged"
  )

  missing <- setdiff(required, names(x))
  if (length(missing)) {
    cli::cli_abort(c(
      "Internal error: constructed TAM fit is missing components.",
      "x" = "Missing element{?s}: {cli::format_inline('{.val {missing}}')}"
    ))
  }

  class(x) <- c("tam_fit", "list")
  x
}


#' @keywords internal
#' @noRd
.is_tam_fit <- function(x) inherits(x, "tam_fit")


#' @keywords internal
#' @noRd
.require_tam_fit <- function(x, arg = "object") {
  if (!.is_tam_fit(x)) {
    cli::cli_abort(c(
      "`{.arg {arg}}` must be a {.cls tam_fit} produced by {.fn fit_tam}.",
      "i" = "Did you call {.fn fit_tam} on your data first?"
    ))
  }
  x
}


.terminal_year <- function(fit) {
  yrs <- fit$dat$years
  proj <- fit$dat$is_proj
  if (length(yrs[!proj])) max(yrs[!proj]) else max(yrs)
}

.terminal_estimate <- function(tab, terminal_year) {
  if (is.null(tab)) {
    return(NA_real_)
  }

  if ("is_proj" %in% names(tab)) {
    tab <- tab[is.na(tab$is_proj) | !tab$is_proj, , drop = FALSE]
  }
  if ("year" %in% names(tab)) {
    tab <- tab[tab$year == terminal_year, , drop = FALSE]
  }
  if (!"est" %in% names(tab)) {
    return(NA_real_)
  }
  if (!nrow(tab)) {
    return(NA_real_)
  }

  est <- tab$est
  if ("age" %in% names(tab)) {
    est <- sum(est, na.rm = TRUE)
  } else {
    est <- est[1]
  }
  est
}


.coef_table <- function(fixed_par) {
  if (is.null(fixed_par) || !nrow(fixed_par)) {
    return(matrix(numeric(0), nrow = 0, ncol = 4,
                  dimnames = list(character(), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))))
  }

  labels <- fixed_par$par
  if ("coef" %in% names(fixed_par) && any(!is.na(fixed_par$coef))) {
    labels <- paste0(labels, ifelse(is.na(fixed_par$coef), "", paste0(" (", fixed_par$coef, ")")))
  }
  if ("year" %in% names(fixed_par) && any(!is.na(fixed_par$year))) {
    labels <- paste0(labels, ifelse(is.na(fixed_par$year), "", paste0(" [year=", fixed_par$year, "]")))
  }
  if ("age" %in% names(fixed_par) && any(!is.na(fixed_par$age))) {
    labels <- paste0(labels, ifelse(is.na(fixed_par$age), "", paste0(" [age=", fixed_par$age, "]")))
  }

  est <- fixed_par$est
  se  <- fixed_par$se
  z   <- est / se
  p   <- 2 * stats::pnorm(-abs(z))

  out <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = p)
  rownames(out) <- labels
  out
}


#' @export
print.tam_fit <- function(x, ...) {
  x <- .require_tam_fit(x, arg = "x")

  cat("Call:\n")
  print(x$call)

  sumry <- summary(x)
  conv <- sumry$convergence

  cat(sprintf("\nObjective: %s\n", format(sumry$objective, digits = 6)))
  cat(sprintf("Converged: %s\n", if (conv$converged) "Yes" else "No"))
  if (!is.null(conv$max_gradient)) {
    cat(sprintf("Max |grad|: %s\n", format(conv$max_gradient, digits = 3, scientific = TRUE)))
  }
  if (!is.null(conv$pd_hessian)) {
    cat(sprintf("pdHess: %s\n", if (conv$pd_hessian) "Yes" else "No"))
  }

  if (nrow(sumry$coefficients)) {
    cat("\nCoefficients:\n")
    stats::printCoefmat(sumry$coefficients, na.print = "NA", ...)
  } else {
    cat("\nCoefficients: (none)\n")
  }

  derived <- sumry$terminal
  if (length(derived)) {
    cat(sprintf("\nTerminal year (%s) estimates:\n", sumry$terminal_year))
    for (nm in names(derived)) {
      cat(sprintf("  %-11s %s\n", paste0(nm, ":"), format(derived[[nm]], digits = 4, scientific = FALSE)))
    }
  }

  invisible(x)
}


#' @export
summary.tam_fit <- function(object, ...) {
  object <- .require_tam_fit(object, arg = "object")

  years <- object$dat$years
  ages  <- object$dat$ages
  obs_n <- sum(object$dat$obs_map$is_observed)

  grad <- object$sdrep$gradient.fixed
  max_grad <- if (!is.null(grad)) max(abs(grad)) else NULL

  terminal_year <- .terminal_year(object)
  metrics <- c("abundance", "recruitment", "ssb", "F_bar", "M_bar")
  terminal_vals <- vapply(metrics, function(nm) {
    .terminal_estimate(object$pop[[nm]], terminal_year)
  }, numeric(1), USE.NAMES = TRUE)
  terminal_vals <- terminal_vals[!is.na(terminal_vals)]
  terminal_vals <- as.list(terminal_vals)

  coeff <- .coef_table(object$fixed_par)

  res <- list(
    call = object$call,
    objective = object$opt$objective,
    convergence = list(
      converged = isTRUE(object$is_converged),
      max_gradient = max_grad,
      pd_hessian = object$sdrep$pdHess
    ),
    data = list(
      years = range(years),
      n_years = length(unique(years)),
      ages = range(ages),
      n_ages = length(unique(ages)),
      n_observations = obs_n
    ),
    coefficients = coeff,
    terminal_year = terminal_year,
    terminal = terminal_vals
  )

  class(res) <- "summary_tam_fit"
  res
}


#' @export
print.summary_tam_fit <- function(x, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) {
    print(x$call)
  }

  cat(sprintf("\nObjective: %s\n", format(x$objective, digits = 6)))
  conv <- x$convergence
  cat(sprintf("Converged: %s\n", if (conv$converged) "Yes" else "No"))
  if (!is.null(conv$max_gradient)) {
    cat(sprintf("Max |grad|: %s\n", format(conv$max_gradient, digits = 3, scientific = TRUE)))
  }
  cat(sprintf("pdHess: %s\n", if (conv$pd_hessian) "Yes" else "No"))

  cat(sprintf("\nObservations: %d  Years: %s-%s  Ages: %s-%s\n",
              x$data$n_observations,
              x$data$years[1], x$data$years[2],
              x$data$ages[1], x$data$ages[2]))

  if (nrow(x$coefficients)) {
    cat("\nCoefficients:\n")
    stats::printCoefmat(x$coefficients, na.print = "NA", ...)
  } else {
    cat("\nCoefficients: (none)\n")
  }

  if (length(x$terminal)) {
    cat(sprintf("\nTerminal year (%s) estimates:\n", x$terminal_year))
    for (nm in names(x$terminal)) {
      cat(sprintf("  %-11s %s\n", paste0(nm, ":"), format(x$terminal[[nm]], digits = 4, scientific = FALSE)))
    }
  }

  invisible(x)
}

