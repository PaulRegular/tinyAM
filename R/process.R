#' Separable AR1 process: density and simulation
#'
#' @title Separable AR1 process error: density and simulation
#'
#' @description
#' Helpers for a two-dimensional (2D) AR1 process, with separable correlation
#' across years and ages. TAM uses these helpers only for `process = "ar1"`;
#' IID processes use Normal densities and draws directly.
#'
#' - [dprocess_ar1()] returns the log-density contribution for a given matrix `x`.
#' - [rprocess_ar1()] generates a matrix draw with the requested dependence structure.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n_y \times n_a}} denote the process on
#' year (\eqn{y}) and age (\eqn{a}) indices.
#' It follows a stationary separable AR(1) in both dimensions with correlations
#' \eqn{\phi_\text{age}} (columns) and \eqn{\phi_\text{year}} (rows), satisfying \eqn{|\phi|<1}.
#'
#' The implied covariance is
#' \deqn{\mathrm{Cov}\{\mathrm{vec}(X)\} =
#'       \frac{\sigma^2}{(1-\phi_\text{age}^2)(1-\phi_\text{year}^2)}
#'       \; \Sigma_\text{age} \otimes \Sigma_\text{year},}
#'
#' where \eqn{\Sigma_\cdot} are AR(1) correlation matrices with entries
#' \eqn{\phi^{|i-j|}}.
#'
#' @param x A numeric matrix (\eqn{n_y \times n_a}) of process residuals
#'   for [dprocess_ar1()].
#' @param ny,na Positive integers: numbers of years and ages for
#'   [rprocess_ar1()].
#' @param phi Length-2 numeric vector \code{c(phi_age, phi_year)} with
#'   values in \eqn{(-1, 1)} for AR(1); zero gives independence.
#' @param sd Positive innovation-scale SD \eqn{\sigma}; the marginal SD is
#'   \eqn{\sigma/\sqrt{(1-\phi_\text{age}^2)(1-\phi_\text{year}^2)}}.
#'
#' @return
#' - [dprocess_ar1()]: a single numeric log-density value.
#' - [rprocess_ar1()]: a numeric matrix of dimension \eqn{n_y \times n_a}.
#'
#' @examples
#' # Simulate then calculate log-density
#' set.seed(1)
#' X_ar1 <- rprocess_ar1(ny = 10, na = 8, sd = 0.3, phi = c(0.5, 0.9))
#' dprocess_ar1(X_ar1, sd = 0.3, phi = c(0.5, 0.9))
#'
#' @seealso
#' [RTMB::dseparable()], [RTMB::dautoreg()], [MASS::mvrnorm()]
#'
#' @import RTMB
#' @rdname process_ar1
#' @export
dprocess_ar1 <- function(x, phi = c(0, 0), sd = 1) {
  phi_age <- phi[1]
  phi_year <- phi[2]
  fa <- function(z) dautoreg(z, phi = phi_age, log = TRUE)
  fy <- function(z) dautoreg(z, phi = phi_year, log = TRUE)
  var <- sd^2 / ((1 - phi_age^2) * (1 - phi_year^2))
  dseparable(fy, fa)(x, scale = sqrt(var))
}

#' @rdname process_ar1
#' @export
rprocess_ar1 <- function(ny, na, phi = c(0, 0), sd = 1) {
  phi_age <- phi[1]
  phi_year <- phi[2]

  ar1_cor <- function(n, phi) stats::toeplitz(phi^(0:(n - 1L)))
  C_age <- ar1_cor(na, phi_age) # columns
  C_year <- ar1_cor(ny, phi_year) # rows

  var <- sd^2 / ((1 - phi_age^2) * (1 - phi_year^2))
  Sigma <- var * kronecker(C_age, C_year)

  z <- MASS::mvrnorm(1L, mu = rep(0, ny * na), Sigma = Sigma)
  return(matrix(z, ny, na))
}


#' Temporal random-walk process: density and simulation
#'
#' @description
#' A random walk allows process deviations to persist and change from one year
#' to the next, independently for each age or age block. There is no tendency
#' to return to a mean and no correlation between ages.
#'
#' @details
#' For a year by age matrix `x`, increments
#' `x[y, a] - x[y - 1, a]` are independent Normal variables with mean zero
#' and standard deviation `sd`. The first row has no increment density.
#' [rprocess_rw()] keeps that starting row and simulates subsequent rows
#' conditionally. In TAM these helpers act on process deviations from the
#' existing mean surfaces (F/M) or cohort predictions (N).
#'
#' @param x A numeric matrix with years in rows and ages or age blocks in columns.
#'   For simulation, only the first row supplies starting values.
#' @param sd Positive scalar standard deviation of year-to-year increments.
#' @return [dprocess_rw()] returns the summed log density of increments;
#'   [rprocess_rw()] returns a matrix with the dimensions and names of `x`.
#' @examples
#' x <- matrix(0, 5, 2)
#' set.seed(1)
#' x <- rprocess_rw(x, sd = 0.2)
#' dprocess_rw(x, sd = 0.2)
#' @rdname process_rw
#' @export
dprocess_rw <- function(x, sd = 1) {
  if (nrow(x) < 2L) {
    return(0)
  }
  increments <- x[-1, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  sum(RTMB::dnorm(increments, 0, sd, log = TRUE))
}

#' @rdname process_rw
#' @export
rprocess_rw <- function(x, sd = 1) {
  if (nrow(x) > 1L) {
    for (y in 2:nrow(x)) {
      x[y, ] <- x[y - 1L, ] + stats::rnorm(ncol(x), 0, sd)
    }
  }
  x
}


#' Random walk with age-correlated increments: density and simulation
#'
#' @description
#' A random walk whose annual changes are correlated between neighbouring ages.
#' Changes are independent between years, with no tendency to return to a mean.
#' tinyAM uses this process for `F_settings$process = "cor_rw"`.
#'
#' @details
#' For a year by age matrix `x`, each vector of annual increments
#' \eqn{x_{t,a}-x_{t-1,a}} has a multivariate Normal distribution with mean zero
#' and covariance \eqn{\Sigma_{ab}=\sigma_a\sigma_b\rho^{|a-b|}}.
#' Distance is measured by column position, not by numerical age labels.
#' Each `sd` is a marginal increment SD; correlation does not inflate it.
#' At `rho = 0`, increments are independent between ages, as in [dprocess_rw()].
#'
#' The first row has no density contribution. [rprocess_cor_rw()] retains it
#' and simulates later rows conditionally. In tinyAM, `x` contains deviations
#' from the log-F mean, rather than the absolute latent log-F states.
#' [dprocess_cor_rw()] supports RTMB automatic differentiation.
#'
#' @param x Numeric matrix with years in rows and ordered ages in columns.
#'   For simulation, only the first row supplies starting values.
#' @param sd Positive marginal SD of annual increments: a scalar shared by all
#'   ages, or a vector with one value per column of `x`.
#' @param rho Correlation between adjacent columns, strictly between -1 and 1.
#'   Negative values allow alternating changes between neighbouring ages.
#' @return [dprocess_cor_rw()] returns the summed log density of annual
#'   increment vectors (zero for one row). [rprocess_cor_rw()] returns a matrix
#'   with the dimensions and names of `x`, leaving its first row unchanged.
#' @examples
#' x <- matrix(0, 6, 3, dimnames = list(2000:2005, 2:4))
#' set.seed(1)
#' draw <- rprocess_cor_rw(x, sd = c(0.2, 0.3, 0.4), rho = 0.6)
#' dprocess_cor_rw(draw, sd = c(0.2, 0.3, 0.4), rho = 0.6)
#' @seealso [dprocess_rw()], [dprocess_ar1()], [tinyAM-model]
#' @rdname process_cor_rw
#' @export
dprocess_cor_rw <- function(x, sd = 1, rho = 0) {
  if (nrow(x) < 2L) {
    return(0)
  }
  increment <- x[-1, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  sd <- rep(sd, length.out = ncol(x))
  standardized <- sweep(increment, 2, sd, "/")
  dprocess_ar1(standardized, phi = c(rho, 0), sd = sqrt(1 - rho^2)) -
    nrow(increment) * sum(log(sd))
}

#' @rdname process_cor_rw
#' @export
rprocess_cor_rw <- function(x, sd = 1, rho = 0) {
  if (nrow(x) < 2L) {
    return(x)
  }
  sd <- rep(sd, length.out = ncol(x))
  for (y in 2:nrow(x)) {
    z <- stats::rnorm(ncol(x))
    if (ncol(x) > 1L) {
      for (a in 2:ncol(x)) {
        z[a] <- rho * z[a - 1L] + sqrt(1 - rho^2) * z[a]
      }
    }
    x[y, ] <- x[y - 1L, ] + sd * z
  }
  x
}

.cor_rw_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  if (!identical(fit$dat$F_settings$process, "cor_rw")) {
    return(out)
  }
  sdr <- fit[["sdrep"]]
  if (!is.list(sdr)) {
    return(out)
  }
  i <- match("atanh_rho_f", names(sdr$par.fixed))
  if (is.na(i)) {
    return(out)
  }
  rho <- tanh(sdr$par.fixed[i])
  if (is.finite(rho) && abs(rho) > .95) {
    out <- rbind(out, data.frame(
      issue = "F_increment_correlation_boundary",
      detail = "F increment age correlation is within 0.05 of -1 or 1. Examine uncertainty and compare independent increments."
    ))
  }
  variance <- if (is.matrix(sdr$cov.fixed)) sdr$cov.fixed[i, i] else NA_real_
  if (is.finite(variance) && variance >= 0 &&
    diff(tanh(sdr$par.fixed[i] + c(-1, 1) * stats::qnorm(.975) * sqrt(variance))) > .8) {
    out <- rbind(out, data.frame(
      issue = "F_increment_correlation_uncertainty",
      detail = "The 95% F increment age-correlation interval spans more than 0.8. Age dependence is imprecise; compare independent increments."
    ))
  }
  out
}


.parse_process_sd <- function(dat, component) {
  settings <- dat[[paste0(component, "_settings")]]
  formula <- settings$sd_form
  if (is.null(formula)) formula <- ~1
  label <- paste0(component, "_settings$sd_form")
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    cli::cli_abort("{label} must be a one-sided formula.")
  }
  default <- identical(formula[[2L]], 1) || identical(formula[[2L]], 1L)
  if (settings$process == "off") {
    if (!default) cli::cli_abort("{label} has no process to describe when process = 'off'.")
    return(NULL)
  }
  walk <- function(x) {
    if (!is.call(x)) {
      return(FALSE)
    }
    .formula_call_name(x) %in% c("iid", "rw", "ar1", "mono", "logistic", "bh", "ricker", "|") ||
      any(vapply(as.list(x)[-1L], walk, logical(1)))
  }
  if (walk(formula[[2L]]) || "year" %in% all.vars(formula)) {
    cli::cli_abort("{label} supports ordinary age-based formulas, not temporal or random-effect SD terms.")
  }
  d <- .population_data(dat, formula, label)
  d$is_proj <- dat$is_proj[match(d$year, dat$years)]
  active <- if (component == "N") {
    dat$ages[-1L]
  } else if (component == "M") {
    as.numeric(names(settings$age_blocks))
  } else {
    dat$ages
  }
  d <- d[d$age %in% active & !d$is_proj, , drop = FALSE]
  frame <- stats::model.frame(formula, d, na.action = stats::na.pass)
  if (anyNA(frame)) cli::cli_abort("{label} needs complete covariates at every active age and historical year.")
  design <- stats::model.matrix(formula, frame)
  if (!ncol(design) || any(!is.finite(design))) cli::cli_abort("{label} must produce a nonempty finite design.")
  for (a in active) {
    rows <- design[d$age == a, , drop = FALSE]
    if (any(rows != rep(rows[1L, ], each = nrow(rows)))) {
      cli::cli_abort("{label} varies across years at age {a}. Temporal process SDs are not enabled.")
    }
  }
  design <- design[match(active, d$age), , drop = FALSE]
  rownames(design) <- as.character(active)
  if (component == "M") {
    for (block in levels(settings$age_blocks)) {
      rows <- design[settings$age_blocks == block, , drop = FALSE]
      if (any(rows != rep(rows[1L, ], each = nrow(rows)))) {
        cli::cli_abort("{label} varies within an M age block. Use finer age_breaks or share the SD within each block.")
      }
    }
    design <- design[as.character(settings$age_block_start), , drop = FALSE]
  }
  if (qr(design)$rank < ncol(design)) cli::cli_abort("{label} is rank deficient on the active process states. Remove redundant terms.")
  list(
    matrix = design, default = default,
    parameter = if (default) paste0("log_sd_", tolower(component)) else paste0("sd_beta_", tolower(component))
  )
}

.process_sd <- function(par, dat, component) {
  design <- dat$process_sd[[component]]
  if (design$default) {
    return(exp(par[[design$parameter]]))
  }
  exp(drop(design$matrix %*% par[[design$parameter]]))
}

.process_sd_surface <- function(par, dat, component) {
  design <- dat$process_sd[[component]]
  matrix(rep(.process_sd(par, dat, component), each = length(dat$years)),
    length(dat$years), nrow(design$matrix),
    dimnames = list(year = dat$years, age = rownames(design$matrix))
  )
}

# Standardization retains the existing AR1 innovation interpretation.
.dprocess_scaled <- function(x, sd, process, phi = c(0, 0)) {
  if (length(sd) == 1L) sd <- rep(sd, ncol(x))
  if (process == "rw") {
    if (nrow(x) < 2L) {
      return(0)
    }
    x <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  }
  z <- sweep(x, 2, sd, "/")
  density <- if (process == "ar1") dprocess_ar1(z, phi = phi, sd = 1) else sum(RTMB::dnorm(z, 0, 1, log = TRUE))
  density - nrow(x) * sum(log(sd))
}

.rprocess_scaled <- function(x, sd, process, phi = c(0, 0)) {
  if (process == "rw") {
    if (nrow(x) > 1L) for (y in 2:nrow(x)) x[y, ] <- x[y - 1L, ] + stats::rnorm(ncol(x), 0, sd)
    return(x)
  }
  z <- if (process == "ar1") {
    rprocess_ar1(nrow(x), ncol(x), phi = phi, sd = 1)
  } else {
    matrix(stats::rnorm(length(x)), nrow(x), ncol(x))
  }
  out <- sweep(z, 2, sd, "*")
  dimnames(out) <- dimnames(x)
  out
}

.process_sd_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  for (component in names(dat$process_sd)) {
    design <- dat$process_sd[[component]]
    if (is.null(design) || design$default) next
    settings <- dat[[paste0(component, "_settings")]]
    rows <- if (component == "M") length(settings$years) else sum(!dat$is_proj) - as.integer(component == "N")
    if (settings$process %in% c("rw", "cor_rw")) rows <- rows - 1L
    support <- colSums(design$matrix != 0) * rows
    if (any(support < 5L)) {
      out <- rbind(out, data.frame(
        issue = "process_SD_support",
        detail = paste0(component, " SD design has fewer than five historical process contributions supporting some coefficients. Share more ages or examine simulation recovery and uncertainty.")
      ))
    }
  }
  out
}

.process_sd_fit_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  for (component in names(fit$dat$process_sd)) {
    d <- fit$pop[[paste0("sd_", component)]]
    if (is.null(d)) next
    if (any(is.finite(d$lwr) & is.finite(d$upr) & d$lwr > 0 & d$upr / d$lwr > 10)) {
      out <- rbind(out, data.frame(
        issue = "process_SD_uncertainty",
        detail = paste0(component, " process SD intervals span more than tenfold at some ages. Examine simpler sharing groups before interpreting age differences.")
      ))
    }
  }
  out
}
