# Resolve targets using raw inputs, before any modeled plus-group pooling.
.hindcast_controls <- function(fit, years, catch, F_mult) {
  supplied <- if (is.null(F_mult)) catch else F_mult
  kind <- if (is.null(catch)) "F_mult" else "catch"
  if (!is.null(catch) && !is.null(F_mult)) cli::cli_abort("Specify catch or F_mult, not both, for hindcasts.")
  observed <- identical(catch, "observed")
  if (is.null(supplied)) supplied <- 1
  if (!observed && (!is.numeric(supplied) || !length(supplied) || any(!is.finite(supplied)) || any(supplied < 0))) {
    cli::cli_abort("Hindcast catch must be 'observed', NULL, or non-negative annual values; F_mult must be non-negative numeric values.")
  }
  if (!observed && !is.null(names(supplied))) {
    named_years <- suppressWarnings(as.numeric(names(supplied)))
    if (anyDuplicated(names(supplied)) || any(!is.finite(named_years)) || any(named_years != trunc(named_years))) {
      cli::cli_abort("Name hindcast controls with unique integer forecast years.")
    }
  } else if (!observed && length(supplied) > 1L) {
    cli::cli_abort("Name a multi-year hindcast catch or F_mult vector by forecast year.")
  }
  out <- data.frame(year = years, control = kind, target = NA_real_, target_source = NA_character_, error = NA_character_)
  source <- fit$refit_args$data
  ages <- sort(unique(c(fit$dat$ages, fit$dat$plus_ages)))
  for (i in seq_along(years)) {
    value <- if (observed) NULL else if (is.null(names(supplied))) supplied else supplied[as.character(years[i])]
    if (length(value) && !is.na(value)) {
      out$target[i] <- value
      out$target_source[i] <- if (kind == "catch") "supplied catch" else "F multiplier"
    } else if (kind == "catch") {
      d <- source$catch[source$catch$year == years[i], ]
      w <- source$weight[source$weight$year == years[i], ]
      counts <- d$obs[match(ages, d$age)]
      weights <- w$obs[match(ages, w$age)]
      if (any(!is.finite(counts)) || any(!is.finite(weights))) {
        out$error[i] <- "Observed catch biomass unavailable: incomplete catch-at-age or weights; supply an annual catch override."
      } else {
        out$target[i] <- sum(counts * weights)
        out$target_source[i] <- "observed catch-at-age x supplied weight"
      }
    } else {
      out$error[i] <- "No F multiplier supplied for this forecast year."
    }
  }
  out
}

# Pair forecasts with the original held-out indices, including when a reference
# fit failed numerical checks. Future observations never enter the fold fit.
.hindcast_index_data <- function(predictions, observations) {
  keys <- c("year", "age", "survey")
  p <- predictions[predictions$is_proj, c(keys, "pred"), drop = FALSE]
  o <- observations[!observations$is_proj, c(keys, "obs"), drop = FALSE]
  p$obs <- NA_real_; p$is_proj <- TRUE
  o$pred <- NA_real_; o$is_proj <- FALSE
  rbind(p, o[, names(p), drop = FALSE])
}
