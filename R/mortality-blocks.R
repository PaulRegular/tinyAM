# Explicit half-open age intervals. Legacy M breaks still use cut_ages().
.explicit_age_blocks <- function(ages, breaks, component) {
  label <- paste0(component, "_settings$age_breaks")
  if (!is.numeric(breaks) || length(breaks) < 2L || any(!is.finite(breaks)) ||
      any(breaks != floor(breaks)) || any(diff(breaks) <= 0) ||
      !breaks[1L] %in% ages || utils::tail(breaks, 1L) != max(ages) + 1L) {
    cli::cli_abort("{label} must contain increasing integer boundaries, starting at a modeled age and ending at max(ages) + 1.")
  }
  active <- ages[ages >= breaks[1L]]
  starts <- utils::head(breaks, -1L); ends <- breaks[-1L] - 1L
  labels <- ifelse(starts == ends, as.character(starts), paste0(starts, "-", ends))
  stats::setNames(factor(labels[findInterval(active, starts)], levels = labels), active)
}

.prepare_f_blocks <- function(dat) {
  if (dat$has_fleet) return(.prepare_fleet_blocks(dat))
  if (is.list(dat$F_settings$age_breaks)) cli::cli_abort("Named-list F age_breaks require catch$fleet labels.")
  settings <- dat$F_settings
  blocks <- if (is.null(settings$age_breaks)) stats::setNames(factor(dat$ages, levels = dat$ages), dat$ages) else
    .explicit_age_blocks(dat$ages, settings$age_breaks, "F")
  settings$age_blocks <- blocks
  settings$age_block_start <- as.numeric(.split_cuts(levels(blocks))$start)
  settings$state_index <- match(as.character(dat$ages), names(blocks))
  settings$state_index <- as.integer(blocks[settings$state_index])
  settings$state_index[is.na(settings$state_index)] <- 0L
  dat$F_settings <- settings
  zero <- dat$obs$catch$age < min(as.numeric(names(blocks)))
  if (any(zero & is.finite(dat$obs$catch$obs) & dat$obs$catch$obs > 0)) {
    cli::cli_abort("Positive catch observations occur at structurally unfished ages below the first F age break. Correct the observations or the F age range.")
  }
  dat
}

.f_state_mean <- function(x, dat) {
  if (dat$has_fleet) {
    rows <- match(paste(dat$F_states$fleet, dat$F_states$age), paste(dat$F_columns$fleet, dat$F_columns$age))
    return(x[, rows, drop = FALSE])
  }
  starts <- dat$F_settings$age_block_start
  if (is.null(starts)) return(x) # Saved fits prepared before block support.
  x[, match(starts, dat$ages), drop = FALSE]
}

.f_blocked <- function(dat) isTRUE(dat$has_fleet) || !is.null(dat$F_settings$age_breaks)

.f_phi <- function(par) {
  "[<-" <- RTMB::ADoverload("[<-")
  value <- par$logit_phi_f
  phi <- c(age = 0, year = 0)
  if ("age" %in% names(value)) phi[1L] <- plogis(value[["age"]])
  if ("year" %in% names(value)) phi[2L] <- plogis(value[["year"]])
  phi
}

.f_rho <- function(par) if (is.null(par$atanh_rho_f)) 0 else tanh(par$atanh_rho_f)

.check_f_blocks <- function(dat) {
  if (!.f_blocked(dat)) return(invisible(NULL))
  parts <- dat$F_modmat
  labels <- paste("mu_form", colnames(parts))
  if (!is.null(dat$F_mono_modmat)) {
    parts <- cbind(parts, dat$F_mono_modmat)
    labels <- c(labels, paste("mu_form", colnames(dat$F_mono_modmat)))
  }
  for (term in dat$F_terms) {
    design <- if (term$type == "logistic") {
      z <- matrix(0, nrow(dat$F_data), 2L * length(term$groups))
      for (i in seq_along(term$groups)) {
        rows <- term$groups[[i]]$rows
        z[rows, 2L * i - 1L] <- 1
        z[rows, 2L * i] <- dat$F_data[[term$variable]][rows]
      }
      z
    } else .q_term_design(term, nrow(dat$F_data))
    colnames(design) <- paste0(term$id, ":", seq_len(ncol(design)))
    parts <- cbind(parts, design)
    labels <- c(labels, rep(paste("mu_form", term$id), ncol(design)))
  }
  if (is.null(dim(parts)) || !ncol(parts)) return(invisible(NULL))
  if (dat$has_fleet) {
    for (f in dat$fleets) {
      proxy <- dat; proxy$has_fleet <- FALSE
      proxy$F_settings$age_blocks <- dat$F_blocks[[f]]
      rows <- dat$F_data$fleet == f
      proxy$F_modmat <- parts[rows, , drop = FALSE]
      proxy$F_mono_modmat <- NULL; proxy$F_terms <- list()
      proxy$F_block_fleet_label <- f
      .check_f_blocks(proxy)
    }
    return(invisible(NULL))
  }
  blocks <- dat$F_settings$age_blocks
  for (block in levels(blocks)) {
    ages <- as.numeric(names(blocks)[blocks == block])
    if (length(ages) < 2L) next
    for (j in seq_len(ncol(parts))) {
      surface <- matrix(parts[, j], length(dat$years), length(dat$ages))
      values <- surface[, match(ages, dat$ages), drop = FALSE]
      if (any(values != values[, 1L])) cli::cli_abort(c(
        "F age block {block} (ages {paste(ages, collapse = ', ')}{if (!is.null(dat$F_block_fleet_label)) paste0('; fleet ', dat$F_block_fleet_label) else ''}) has a conflicting {labels[j]} component.",
        "i" = "Share the mean within this block or refine F_settings$age_breaks."
      ))
    }
  }
  invisible(NULL)
}

.f_state_rows <- function(dat) {
  if (dat$has_fleet) return(paste(dat$F_data$fleet, dat$F_data$age) %in% paste(dat$F_states$fleet, dat$F_states$age))
  starts <- dat$F_settings$age_block_start
  dat$F_data$age %in% if (is.null(starts)) dat$ages else starts
}
