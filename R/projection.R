# Projection controls retain calendar-year names after preparation.
.projection_vector <- function(x, years, label) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0)) {
    what <- if (label == "F_mult") "multipliers" else "catch biomass values"
    cli::cli_abort("proj_settings${label} must contain finite, non-negative {what}.")
  }
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), as.character(years))) {
      cli::cli_abort("Named proj_settings${label} must identify every projection year exactly once.")
    }
    x <- x[as.character(years)]
  } else if (length(x) == 1L) {
    x <- rep(x, length(years))
  }
  if (length(x) != length(years)) cli::cli_abort("length(proj_settings${label}) must equal n_proj, or be one for a constant projection.")
  stats::setNames(x, years)
}

# The continuous Baranov fraction includes its limit at zero total mortality.
.baranov_fraction <- function(Z) {
  if (inherits(Z, "advector")) return(-expm1(-Z) / Z)
  out <- rep(1, length(Z))
  positive <- Z != 0
  out[positive] <- -expm1(-Z[positive]) / Z[positive]
  out
}

.catch_biomass <- function(k, N, W, F_ref, M) {
  F <- k * F_ref
  sum(N * W * F * .baranov_fraction(F + M))
}

.solve_catch <- function(target, N, W, F_ref, M, year) {
  if (target == 0) return(0)
  biomass <- sum(N * W)
  ad <- inherits(biomass, "advector") || inherits(F_ref, "advector") || inherits(M, "advector")
  if (!ad && (!is.finite(biomass) || target >= biomass || any(F_ref <= 0))) {
    cli::cli_abort("Catch projection for {year} is infeasible: target {signif(target, 6)} must be below vulnerable biomass {signif(biomass, 6)} with positive terminal F.")
  }
  # RTMB root bounds are data, not AD parameters. Use a wide finite bracket;
  # an unresolved root is an error, never a capped fishing mortality.
  upper <- 1e8
  residual <- function(k) .catch_biomass(k, N, W, F_ref, M) - target
  if (!ad && residual(upper) <= 0) cli::cli_abort("Could not bracket catch projection for {year}; target was not changed.")
  # A log multiplier gives relative accuracy across catch/abundance scales.
  f <- function(log_k) residual(exp(log_k)) / target
  root <- RTMB::uniroot(f, lower = -700, upper = log(upper), tol = 1e-10, maxiter = 200)$root
  k <- exp(root)
  if (!ad && (!is.finite(k) || abs(residual(k)) > 1e-8 * max(1, target))) {
    cli::cli_abort("Catch projection solver did not converge for {year}; requested catch was not changed.")
  }
  k
}
