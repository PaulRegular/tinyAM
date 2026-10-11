# Annual F increments have age correlation, without temporal stationarity.
.dprocess_cor_rw <- function(x, sd = 1, rho = 0) {
  if (nrow(x) < 2L) return(0)
  increment <- x[-1, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  sd <- rep(sd, length.out = ncol(x))
  standardized <- sweep(increment, 2, sd, "/")
  dprocess_ar1(standardized, phi = c(rho, 0), sd = sqrt(1 - rho^2)) -
    nrow(increment) * sum(log(sd))
}

.rprocess_cor_rw <- function(x, sd = 1, rho = 0) {
  if (nrow(x) < 2L) return(x)
  sd <- rep(sd, length.out = ncol(x))
  for (y in 2:nrow(x)) {
    z <- stats::rnorm(ncol(x))
    if (ncol(x) > 1L) for (a in 2:ncol(x)) {
      z[a] <- rho * z[a - 1L] + sqrt(1 - rho^2) * z[a]
    }
    x[y, ] <- x[y - 1L, ] + sd * z
  }
  x
}

.cor_rw_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  if (!identical(fit$dat$F_settings$process, "cor_rw")) return(out)
  sdr <- fit[["sdrep"]]
  if (!is.list(sdr)) return(out)
  i <- match("atanh_rho_f", names(sdr$par.fixed))
  if (is.na(i)) return(out)
  rho <- tanh(sdr$par.fixed[i])
  if (is.finite(rho) && abs(rho) > .95) out <- rbind(out, data.frame(
    issue = "F_increment_correlation_boundary",
    detail = "F increment age correlation is within 0.05 of -1 or 1. Examine uncertainty and compare independent increments."))
  variance <- if (is.matrix(sdr$cov.fixed)) sdr$cov.fixed[i, i] else NA_real_
  if (is.finite(variance) && variance >= 0 &&
      diff(tanh(sdr$par.fixed[i] + c(-1, 1) * stats::qnorm(.975) * sqrt(variance))) > .8) {
    out <- rbind(out, data.frame(issue = "F_increment_correlation_uncertainty",
      detail = "The 95% F increment age-correlation interval spans more than 0.8. Age dependence is imprecise; compare independent increments."))
  }
  out
}
