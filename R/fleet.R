# F gains fleet/age columns; N and M retain one population-covariate frame.
.prepare_fleets <- function(dat) {
  dat$has_fleet <- "fleet" %in% names(dat$obs$catch)
  fleets <- if (dat$has_fleet) {
    x <- dat$obs$catch$fleet
    if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
  } else "Fishery"
  dat$fleets <- fleets
  dat$F_data <- if (dat$has_fleet) {
    do.call(rbind, lapply(fleets, function(f) {
      x <- dat$pop_data
      x$fleet <- factor(f, levels = fleets)
      x
    }))
  } else dat$pop_data
  dat$F_columns <- expand.grid(age = dat$ages, fleet = fleets, stringsAsFactors = FALSE)
  dat$F_columns$key <- if (length(fleets) == 1L) {
    as.character(dat$ages)
  } else paste(dat$F_columns$fleet, dat$F_columns$age, sep = " / age ")
  dat
}

.F_matrix <- function(x, dat, historical = FALSE) {
  years <- if (historical) dat$years[!dat$is_proj] else dat$years
  matrix(x, length(years), nrow(dat$F_columns),
    dimnames = list(year = years, age = dat$F_columns$key))
}

.F_total <- function(x, dat) {
  A <- length(dat$ages)
  ans <- x[, seq_len(A), drop = FALSE]
  if (length(dat$fleets) > 1L) for (f in 2:length(dat$fleets)) {
    ans <- ans + x[, (f - 1L) * A + seq_len(A), drop = FALSE]
  }
  dimnames(ans) <- list(year = dat$years, age = dat$ages)
  ans
}

.F_array <- function(x, dat) {
  array(x, dim = c(length(dat$years), length(dat$ages), length(dat$fleets)),
    dimnames = list(year = dat$years, age = dat$ages, fleet = dat$fleets))
}

# Reuse each existing density on one fleet at a time: age correlation never
# crosses the boundary between the final age of one fleet and the next fleet.
.F_density <- function(states, mean, par, dat) {
  type <- dat$F_settings$process
  if (type == "off") return(0)
  result <- 0
  sd <- .process_sd(par, dat, "F")
  A <- length(dat$ages)
  for (f in seq_along(dat$fleets)) {
    cols <- (f - 1L) * A + seq_len(A)
    x <- states[, cols, drop = FALSE] - mean[!dat$is_proj, cols, drop = FALSE]
    scale <- if (length(sd) == 1L) sd else sd[cols]
    value <- if (!dat$process_sd$F$default && type != "cor_rw") {
      .dprocess_scaled(x, scale, type,
        phi = if (type == "ar1") plogis(par$logit_phi_f) else c(0, 0))
    } else if (type == "rw") {
      dprocess_rw(x, sd = scale)
    } else if (type == "cor_rw") {
      dprocess_cor_rw(x, sd = scale, rho = tanh(par$atanh_rho_f))
    } else if (type == "iid") {
      sum(RTMB::dnorm(x, 0, scale, log = TRUE))
    } else {
      dprocess_ar1(x, sd = scale, phi = plogis(par$logit_phi_f))
    }
    result <- result + value
  }
  result
}

.F_states <- function(par, dat, mean, simulate = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")
  type <- dat$F_settings$process
  states <- if (type == "off") mean[!dat$is_proj, , drop = FALSE] else par[["log_f"]]
  if (simulate && type != "off") {
    sd <- .process_sd(par, dat, "F")
    A <- length(dat$ages)
    for (f in seq_along(dat$fleets)) {
      cols <- (f - 1L) * A + seq_len(A)
      mu <- mean[!dat$is_proj, cols, drop = FALSE]
      x <- states[, cols, drop = FALSE] - mu
      scale <- if (length(sd) == 1L) sd else sd[cols]
      states[, cols] <- mu + if (!dat$process_sd$F$default && type != "cor_rw") {
        .rprocess_scaled(x, scale, type,
          phi = if (type == "ar1") plogis(par$logit_phi_f) else c(0, 0))
      } else if (type == "rw") {
        rprocess_rw(x, sd = scale)
      } else if (type == "cor_rw") {
        rprocess_cor_rw(x, sd = scale, rho = tanh(par$atanh_rho_f))
      } else if (type == "iid") {
        matrix(stats::rnorm(length(x), 0, scale), nrow(x), ncol(x))
      } else {
        rprocess_ar1(nrow(x), ncol(x), sd = scale, phi = plogis(par$logit_phi_f))
      }
    }
  }
  full <- mean
  full[!dat$is_proj, ] <- states
  if (any(dat$is_proj)) {
    terminal <- max(which(!dat$is_proj))
    k <- if (is.null(dat$proj_settings$F_mult)) rep(1, sum(dat$is_proj)) else dat$proj_settings$F_mult
    full[dat$is_proj, ] <- sweep(full[rep(terminal, sum(dat$is_proj)), , drop = FALSE], 1L, log(k), "+")
  }
  list(states = states, full = full)
}
