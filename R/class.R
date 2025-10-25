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

.coef_table <- function(fixed_par) {
  interval <- attr(fixed_par, "interval", exact = TRUE)
  if (is.null(interval) || !is.finite(interval)) {
    interval <- 0.95
  }
  ci_level <- formatC(interval * 100, format = "fg")
  ci_names <- c(sprintf("Lower %s%%", ci_level), sprintf("Upper %s%%", ci_level))

  if (is.null(fixed_par) || !nrow(fixed_par)) {
    return(matrix(numeric(0), nrow = 0, ncol = 4,
                  dimnames = list(character(), c("Estimate", "Std. Error", ci_names))))
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
  lower <- if ("lwr" %in% names(fixed_par)) fixed_par$lwr else rep(NA_real_, length(est))
  upper <- if ("upr" %in% names(fixed_par)) fixed_par$upr else rep(NA_real_, length(est))

  out <- cbind(Estimate = est, `Std. Error` = se, lower, upper)
  colnames(out)[3:4] <- ci_names
  rownames(out) <- labels
  out
}



.terminal_table <- function(pop, terminal_year) {
  default_interval <- attr(pop, "interval", exact = TRUE)
  if (is.null(default_interval) || !is.finite(default_interval)) {
    default_interval <- 0.95
  }

  default_ci_level <- formatC(default_interval * 100, format = "fg")
  default_ci_names <- c(sprintf("Lower %s%%", default_ci_level),
                        sprintf("Upper %s%%", default_ci_level))

  empty <- matrix(numeric(0), nrow = 0, ncol = 4,
                  dimnames = list(character(),
                                  c("Estimate", "Std. Error", default_ci_names)))

  if (is.null(pop)) {
    return(empty)
  }

  metrics <- c("abundance", "recruitment", "ssb", "F_bar", "M_bar")
  rows <- lapply(metrics, function(nm) {
    tab <- pop[[nm]]
    if (is.null(tab)) {
      return(NULL)
    }

    interval <- attr(tab, "interval", exact = TRUE)
    if (is.null(interval) || !is.finite(interval)) {
      interval <- default_interval
    }
    ci_level <- formatC(interval * 100, format = "fg")
    ci_names <- c(sprintf("Lower %s%%", ci_level), sprintf("Upper %s%%", ci_level))

    if ("is_proj" %in% names(tab)) {
      tab <- tab[is.na(tab$is_proj) | !tab$is_proj, , drop = FALSE]
    }
    if ("year" %in% names(tab)) {
      tab <- tab[tab$year == terminal_year, , drop = FALSE]
    }

    if (!nrow(tab) || !"est" %in% names(tab)) {
      return(NULL)
    }

    if (nrow(tab) > 1 && "age" %in% names(tab)) {
      totals <- tab[is.na(tab$age), , drop = FALSE]
      if (nrow(totals) >= 1) {
        tab <- totals[1, , drop = FALSE]
      } else {
        tab <- tab[1, , drop = FALSE]
      }
    } else {
      tab <- tab[1, , drop = FALSE]
    }

    se <- if ("sd" %in% names(tab)) tab$sd else if ("se" %in% names(tab)) tab$se else NA_real_
    lower <- if ("lwr" %in% names(tab)) tab$lwr else NA_real_
    upper <- if ("upr" %in% names(tab)) tab$upr else NA_real_

    matrix(c(tab$est, se, lower, upper), nrow = 1,
           dimnames = list(nm, c("Estimate", "Std. Error", ci_names)))
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(empty)
  }

  mat <- do.call(rbind, rows)
  mat
}


.format_numbers <- function(x, digits, big_mark = NULL) {
  n <- length(x)
  if (!n) {
    return(character())
  }

  digits <- rep_len(as.integer(digits), n)

  if (is.null(big_mark)) {
    big_mark <- rep_len("", n)
  } else {
    big_mark <- rep_len(as.character(big_mark), n)
    big_mark[is.na(big_mark)] <- ""
  }

  out <- rep(NA_character_, n)
  is_na <- is.na(x)
  is_inf <- !is_na & !is.finite(x)
  out[is_inf] <- as.character(x[is_inf])

  keep_idx <- which(!is_na & !is_inf)
  if (length(keep_idx)) {
    out[keep_idx] <- mapply(
      function(val, dig, mark) {
        formatC(val, format = "f", digits = dig, big.mark = mark)
      },
      x[keep_idx], digits[keep_idx], big_mark[keep_idx],
      USE.NAMES = FALSE
    )
  }

  out
}


.print_formatted <- function(tab, ...) {
  dots <- list(...)
  if (!"quote" %in% names(dots)) {
    dots$quote <- FALSE
  }
  do.call(print, c(list(tab), dots))
}


.format_terminal_display <- function(tab) {
  if (!is.matrix(tab) || !nrow(tab)) {
    return(tab)
  }

  metric_digits0 <- c("abundance", "recruitment", "ssb")
  row_names <- rownames(tab)
  col_names <- colnames(tab)
  is_sd_col <- col_names == "Std. Error"

  formatted <- matrix(NA_character_, nrow = nrow(tab), ncol = ncol(tab),
                      dimnames = dimnames(tab))

  for (i in seq_len(nrow(tab))) {
    row_nm <- row_names[i]
    row_key <- tolower(row_nm)
    digits_base <- if (row_key %in% metric_digits0) 0L else 3L
    big_mark <- if (digits_base == 0L) "," else NULL

    digits <- ifelse(is_sd_col, 3L, digits_base)
    big <- ifelse(is_sd_col, "", if (is.null(big_mark)) "" else big_mark)
    formatted[i, ] <- .format_numbers(tab[i, ], digits = digits, big_mark = big)
  }

  formatted
}


.format_coef_display <- function(tab) {
  if (!is.matrix(tab) || !nrow(tab)) {
    return(tab)
  }

  row_names <- rownames(tab)
  col_names <- colnames(tab)
  is_sd_col <- col_names == "Std. Error"

  formatted <- matrix(NA_character_, nrow = nrow(tab), ncol = ncol(tab),
                      dimnames = dimnames(tab))

  for (i in seq_len(nrow(tab))) {
    row_nm <- row_names[i]
    is_r0 <- grepl("\\br0\\b", row_nm)
    digits_base <- if (is_r0) 0L else 3L
    big_mark <- if (is_r0 && digits_base == 0L) "," else NULL

    digits <- ifelse(is_sd_col, 3L, digits_base)
    big <- ifelse(is_sd_col, "", if (is.null(big_mark)) "" else big_mark)
    formatted[i, ] <- .format_numbers(tab[i, ], digits = digits, big_mark = big)
  }

  formatted
}


#' @export
print.tam_fit <- function(x, ...) {
  x <- .require_tam_fit(x, arg = "x")

  cat("Call:\n")
  print(x$call)

  sumry <- summary(x)
  conv <- sumry$convergence

  cat(sprintf("\nObjective: %s\n", format(sumry$objective, digits = 6)))
  if (!is.null(conv$max_gradient)) {
    cat(sprintf("Max |grad|: %s\n", format(conv$max_gradient, digits = 3, scientific = TRUE)))
  }
  if (!is.null(conv$pd_hessian)) {
    cat(sprintf("pdHess: %s\n", if (conv$pd_hessian) "Yes" else "No"))
  }

  if (nrow(sumry$coefficients)) {
    cat("\nCoefficients:\n")
    coef_tab <- .format_coef_display(sumry$coefficients)
    .print_formatted(coef_tab, ...)
  } else {
    cat("\nCoefficients: (none)\n")
  }

  derived <- sumry$terminal_vals
  if (nrow(derived) > 0) {
    cat(sprintf("\nTerminal year (%s) estimates:\n", sumry$terminal_year))
    terminal_tab <- .format_terminal_display(derived)
    .print_formatted(terminal_tab, ...)
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
  terminal_vals <- .terminal_table(object$pop, terminal_year)
  terminal_list <- if (length(terminal_vals) > 0) {
    stats::setNames(as.list(terminal_vals[, "Estimate"]), rownames(terminal_vals))
  } else {
    list()
  }

  coeff <- .coef_table(object$fixed_par)

  res <- list(
    call = object$call,
    objective = object$opt$objective,
    convergence = list(
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
    terminal = terminal_list,
    terminal_vals = terminal_vals
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
    coef_tab <- .format_coef_display(x$coefficients)
    .print_formatted(coef_tab, ...)
  } else {
    cat("\nCoefficients: (none)\n")
  }

  if (nrow(x$terminal_vals) > 0) {
    cat(sprintf("\nTerminal year (%s) estimates:\n", x$terminal_year))
    terminal_tab <- .format_terminal_display(x$terminal_vals)
    .print_formatted(terminal_tab, ...)
  }

  invisible(x)
}

