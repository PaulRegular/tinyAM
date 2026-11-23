.obs_cor_structure <- function(x, series, setting) {
  if (!is.character(x) || !length(x) || !is.null(dim(x)) ||
      !all(x %in% c("independent", "ar1"))) {
    cli::cli_abort("{setting} must use {.val independent} or {.val ar1}.")
  }
  if (is.null(names(x))) {
    if (length(x) != 1L) cli::cli_abort("{setting} must be a scalar or name each selected series.")
    return(stats::setNames(rep(x, length(series)), series))
  }
  if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    cli::cli_abort("{setting} must have unique, nonempty series names.")
  }
  unknown <- setdiff(names(x), series)
  if (length(unknown)) cli::cli_abort("{setting} has unknown series {.val {unknown}}.")
  out <- stats::setNames(rep("independent", length(series)), series)
  out[names(x)] <- x
  out
}

.prepare_obs_correlation <- function(dat, surveys = unique(as.character(dat$obs$index$survey))) {
  catch <- .obs_cor_structure(dat$catch_settings$cor_struct, dat$fleets, "catch_settings$cor_struct")
  index <- .obs_cor_structure(dat$index_settings$cor_struct, surveys, "index_settings$cor_struct")
  if (!any(c(catch, index) == "ar1")) return(NULL)

  d <- dat$obs_map
  series <- ifelse(d$type == "catch", if (dat$has_fleet) as.character(d$fleet) else "Fishery",
                   as.character(d$survey))
  lognormal <- if (is.null(d$family)) rep(TRUE, nrow(d)) else d$family == "lognormal"
  groups <- character()
  blocks <- list()
  for (type in c("catch", "index")) {
    structures <- if (type == "catch") catch else index
    for (name in names(structures)[structures == "ar1"]) {
      rows <- which(d$type == type & series == name)
      if (any(!lognormal[rows])) {
        cli::cli_abort("Series {.val {name}} requires independent observation errors for its family. Remove its {.code cor_struct = 'ar1'} setting.")
      }
      rows <- rows[dat$is_observed[rows]]
      by_year <- split(rows, d$year[rows])
      by_year <- by_year[lengths(by_year) > 1L]
      if (!length(by_year)) {
        cli::cli_warn("{type} series {.val {name}} has no within-year age pairs. Observation errors are independent for this fit.")
        next
      }
      groups <- c(groups, paste(type, name, sep = ":"))
      blocks <- c(blocks, lapply(by_year, function(rows) {
        rows <- rows[order(d$age[rows])]
        list(rows = rows, ages = d$age[rows], group = length(groups))
      }))
    }
  }
  if (!length(groups)) return(NULL)

  # Likelihood indices refer to the compact vector of observed values.
  observed <- which(dat$is_observed & lognormal)
  positions <- lapply(blocks, function(b) match(b$rows, observed))
  from <- unlist(lapply(positions, utils::head, -1L), use.names = FALSE)
  to <- unlist(lapply(positions, function(i) i[-1L]), use.names = FALSE)
  list(groups = groups, blocks = blocks, positions = positions,
       from = from, to = to, first = setdiff(seq_along(observed), to),
       uncorrelated = setdiff(seq_along(observed), unlist(positions, use.names = FALSE)),
       gap = unlist(lapply(blocks, function(b) diff(b$ages)), use.names = FALSE),
       group = rep(vapply(blocks, function(b) b$group, integer(1)), lengths(positions) - 1L))
}

.obs_has_correlation <- function(dat) length(dat$obs_correlation$groups) > 0L

.obs_cor_nll <- function(dat, par, observed, mean, sd) {
  if (!.obs_has_correlation(dat)) {
    return(-sum(RTMB::dnorm(observed, mean, sd = sd, log = TRUE)))
  }
  cor <- dat$obs_correlation
  log_rho <- -RTMB::logspace_add(0, -par$logit_rho_obs)
  if (inherits(observed, "osa")) {
    # Preserve RTMB's observation markers for OSA's partial Gaussian likelihood.
    i <- cor$uncorrelated
    nll <- if (length(i)) -sum(RTMB::dnorm(observed[i], mean[i], sd = sd[i], log = TRUE)) else 0
    for (b in seq_along(cor$blocks)) {
      block <- cor$blocks[[b]]
      i <- cor$positions[[b]]
      size <- length(i)
      row <- rep(seq_len(size), size)
      col <- rep(seq_len(size), each = size)
      gap <- abs(block$ages[row] - block$ages[col])
      covariance <- matrix(sd[i[row]] * sd[i[col]] * exp(gap * log_rho[block$group]),
                           size, size)
      nll <- nll - RTMB::dmvnorm(observed[i], mean[i], covariance, log = TRUE)
    }
    return(nll)
  }

  # The standardized errors are Markov in age, including across missing ages.
  z <- (observed - mean) / sd
  log_phi <- log_rho[cor$group] * cor$gap
  first <- cor$first
  to <- cor$to
  -sum(RTMB::dnorm(observed[first], mean[first], sd = sd[first], log = TRUE)) -
    sum(RTMB::dnorm(z[to], exp(log_phi) * z[cor$from],
                    sd = sqrt(-expm1(2 * log_phi)), log = TRUE) - log(sd[to]))
}

.draw_log_observations <- function(dat, par, mean, sd, predictions = NULL) {
  draw <- dat$is_observed
  log_obs <- dat$log_obs
  if (!is.null(predictions)) {
    family <- dat$obs_map$family
    draw <- draw & family == "lognormal"
    delta <- which(dat$is_observed & .is_delta_family(family))
    positive <- delta[log(stats::runif(length(delta))) < predictions$log_p_positive[delta]]
    log_obs[delta] <- -Inf
    lognormal <- positive[family[positive] == "delta_lognormal"]
    log_obs[lognormal] <- stats::rnorm(length(lognormal),
      predictions$log_pred_positive[lognormal], sd[lognormal])
    gamma <- c(positive[family[positive] != "delta_lognormal"],
      which(dat$is_observed & family %in% c("gamma", "gengamma")))
    for (kind in unique(family[gamma])) {
      rows <- gamma[family[gamma] == kind]
      log_obs[rows] <- .draw_positive_log_obs(kind, predictions$log_obs_positive_mean[rows],
        sd[rows], predictions$observation_Q[rows])
    }
    tweedie <- which(dat$is_observed & family == "tweedie")
    if (length(tweedie)) {
      log_obs[tweedie] <- .draw_tweedie_log_obs(mean[tweedie], sd[tweedie],
        predictions$tweedie_power[tweedie])
    }
  }
  if (!.obs_has_correlation(dat)) {
    log_obs[draw] <- stats::rnorm(sum(draw), mean[draw], sd[draw])
    return(log_obs)
  }

  z <- numeric(length(sd))
  z[draw] <- stats::rnorm(sum(draw))
  log_rho <- stats::plogis(par$logit_rho_obs, log.p = TRUE)
  for (block in dat$obs_correlation$blocks) {
    rows <- block$rows
    log_phi <- log_rho[block$group] * diff(block$ages)
    for (i in seq_along(log_phi)) {
      z[rows[i + 1L]] <- exp(log_phi[i]) * z[rows[i]] +
        sqrt(-expm1(2 * log_phi[i])) * z[rows[i + 1L]]
    }
  }
  log_obs[draw] <- mean[draw] + z[draw] * sd[draw]
  log_obs
}

.draw_tweedie_log_obs <- function(log_mean, sd, power) {
  sd <- rep_len(sd, length(log_mean))
  power <- rep_len(power, length(log_mean))
  log_phi <- 2 * log(sd)
  log_intensity <- (2 - power) * log_mean - log_phi - log(2 - power)
  count <- stats::rpois(length(log_mean), exp(log_intensity))
  out <- rep(-Inf, length(count))
  out[is.na(count)] <- NA_real_
  positive <- which(count > 0)
  shape <- count[positive] * (2 - power[positive]) / (power[positive] - 1)
  log_scale <- log_phi[positive] + log(power[positive] - 1) +
    (power[positive] - 1) * log_mean[positive]
  out[positive] <- log_scale + log(shape) + .log_gamma_relative_draw(shape)
  out
}
