# Catch rows are the explanatory grid for F, including unobserved catches.
.check_f_plus_covariates <- function(obs, plus_age, settings, years = unique(obs$catch$year)) {
  mean <- settings$mu_form
  parts <- list()
  if (!is.null(mean) && .contains_scope(mean)) {
    split <- .split_scopes(mean,obs$catch,"fleet")
    parts <- lapply(split$scopes,function(x) list(formula=x$formula,fleet=x$name))
    parts[[length(parts)+1L]] <- list(formula=split$formula,fleet=NULL)
  } else parts <- list(list(formula=mean,fleet=NULL))
  forms <- c(list(settings$sd_form),
    if (inherits(settings$cor_form, "formula")) list(settings$cor_form) else settings$cor_form)
  parts <- c(parts,lapply(forms,function(x) list(formula=x,fleet=NULL)))
  for (part in parts) {
    variables <- setdiff(.formula_input_variables(part$formula),c("year","age","fleet"))
    if ("obs" %in% variables) cli::cli_abort("F formulas cannot use the catch response obs as a predictor.")
    variables <- intersect(variables,names(obs$catch))
    d <- obs$catch[obs$catch$age >= plus_age & obs$catch$year %in% years,,drop=FALSE]
    if (!is.null(part$fleet)) d <- d[d$fleet == part$fleet,,drop=FALSE]
    group <- paste(d$year,if("fleet" %in% names(d)) d$fleet else "Fishery")
    for (name in variables) for (g in unique(group)) {
      if (length(unique(d[[name]][group == g])) > 1L) {
        cli::cli_abort("Catch covariate {name} conflicts across biological ages pooled into plus age {plus_age} in {g}. Supply a consistent predictor for the modeled plus group.")
      }
    }
  }
  invisible(NULL)
}

.f_data <- function(dat, formula = NULL, label = "F formula") {
  d <- dat$F_data
  if (!is.null(formula)) {
    if ("obs" %in% .formula_input_variables(formula)) cli::cli_abort("{label} cannot use the catch response obs as a predictor.")
    # The common population checker also understands formula marker arguments.
    proxy <- dat; proxy$pop_data <- d
    .population_data(proxy, formula, label, table = "data$catch")
  }
  d
}

.prepare_fleets <- function(dat) {
  dat$has_fleet <- "fleet" %in% names(dat$obs$catch)
  dat$fleets <- if (dat$has_fleet) {
    x <- dat$obs$catch$fleet
    if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
  } else "Fishery"
  d <- dat$obs$catch
  d <- d[order(if (dat$has_fleet) match(d$fleet, dat$fleets) else rep(1L, nrow(d)), d$age, d$year), , drop = FALSE]
  if (dat$has_fleet && !is.factor(d$fleet)) d$fleet <- factor(d$fleet, levels = dat$fleets)
  dat$obs$catch <- d
  dat$F_data <- d
  dat$F_columns <- expand.grid(age = dat$ages, fleet = dat$fleets, stringsAsFactors = FALSE)
  dat$F_columns$key <- if (!dat$has_fleet) as.character(dat$ages) else
    paste(dat$F_columns$fleet, dat$F_columns$age, sep = " / age ")
  dat
}

.prepare_fleet_blocks <- function(dat) {
  breaks <- dat$F_settings$age_breaks
  if (is.list(breaks) && (is.null(names(breaks)) || anyDuplicated(names(breaks)) ||
      !setequal(names(breaks), dat$fleets))) {
    cli::cli_abort("A list of F age_breaks must name every catch fleet exactly once; use NULL for unblocked fleets.")
  }
  states <- list(); indices <- integer(); offset <- 0L
  dat$F_blocks <- stats::setNames(lapply(dat$fleets, function(f) {
    b <- if (is.list(breaks)) breaks[[f]] else breaks
    if (is.null(b)) stats::setNames(factor(dat$ages, levels = dat$ages), dat$ages) else
      .explicit_age_blocks(dat$ages, b, "F")
  }), dat$fleets)
  for (f in dat$fleets) {
    blocks <- dat$F_blocks[[f]]
    starts <- as.numeric(.split_cuts(levels(blocks))$start)
    states[[f]] <- data.frame(fleet = f, age = starts, block = levels(blocks),
      key = paste(f, levels(blocks), sep = " / ages "))
    i <- as.integer(blocks[match(as.character(dat$ages), names(blocks))])
    indices <- c(indices, ifelse(is.na(i), 0L, i + offset))
    offset <- offset + length(starts)
  }
  dat$F_states <- do.call(rbind, states); rownames(dat$F_states) <- NULL
  dat$F_settings$state_index <- indices
  dat$F_settings$age_block_start <- dat$F_states$age
  dat$F_settings$age_blocks <- factor(dat$F_states$key, levels = dat$F_states$key)
  zero <- indices[match(paste(dat$obs$catch$fleet, dat$obs$catch$age),
    paste(dat$F_columns$fleet, dat$F_columns$age))] == 0L
  if (any(zero & is.finite(dat$obs$catch$obs) & dat$obs$catch$obs > 0)) {
    cli::cli_abort("Positive catch observations occur at structurally unfished ages. Correct the observations or fleet F age_breaks.")
  }
  dat
}

.F_matrix <- function(x, dat, historical = FALSE) {
  years <- if (historical) dat$years[!dat$is_proj] else dat$years
  matrix(x, length(years), nrow(dat$F_columns),
    dimnames = list(year = years, age = dat$F_columns$key))
}

.F_total <- function(x, dat) {
  A <- length(dat$ages)
  out <- x[, seq_len(A), drop = FALSE]
  if (length(dat$fleets) > 1L) for (f in 2:length(dat$fleets)) {
    out <- out + x[, (f - 1L) * A + seq_len(A), drop = FALSE]
  }
  dimnames(out) <- list(year = dat$years, age = dat$ages)
  out
}

# Log total F without exp/log round trips in the tails of selectivity curves.
# Inactive ages use a finite placeholder; their reported F remains exactly zero.
.F_log_total <- function(x, dat, active = dat$F_settings$state_index > 0L) {
  "[<-" <- RTMB::ADoverload("[<-")
  out <- matrix(0, length(dat$years), length(dat$ages),
    dimnames = list(year = dat$years, age = dat$ages))
  for (a in seq_along(dat$ages)) {
    columns <- which(dat$F_columns$age == dat$ages[a] & active)
    if (length(columns)) out[, a] <- Reduce(RTMB::logspace_add,
      lapply(columns, function(j) x[, j]))
  }
  out
}

.F_array <- function(x, dat) array(x,
  c(length(dat$years), length(dat$ages), length(dat$fleets)),
  dimnames = list(year = dat$years, age = dat$ages, fleet = dat$fleets))

# Correlation formulas select shared categorical groups, never regress rho.
.prepare_f_cor <- function(dat) {
  form <- dat$F_settings$cor_form
  if (is.null(form)) form <- ~1
  if (inherits(form, "formula")) form <- list(age = form, year = form)
  if (!is.list(form) || is.null(names(form)) || anyDuplicated(names(form)) ||
      !setequal(names(form), c("age", "year"))) {
    cli::cli_abort("F_settings$cor_form must be a formula, or list(age = ~ group, year = ~ group).")
  }
  dat$F_cor <- lapply(c("age", "year"), function(dimension) {
    formula <- form[[dimension]]
    if (!inherits(formula, "formula") || length(formula) != 2L ||
        any(c("year", "age", "obs") %in% all.vars(formula))) {
      cli::cli_abort("F correlation groups must be categorical and constant within each fleet; age and year cannot define groups.")
    }
    d <- .f_data(dat, formula, "F_settings$cor_form")
    frame <- stats::model.frame(formula, d, na.action = stats::na.pass)
    if (anyNA(frame) || any(vapply(frame, function(x) !is.factor(x) && !is.character(x) && !is.logical(x), logical(1)))) {
      cli::cli_abort("F_settings$cor_form requires complete categorical groups, not numeric predictors.")
    }
    values <- if (ncol(frame)) do.call(paste, c(frame, sep = " / ")) else rep("shared", nrow(d))
    fleet <- if (dat$has_fleet) as.character(d$fleet) else rep("Fishery", nrow(d))
    groups <- vapply(dat$fleets, function(f) {
      z <- unique(values[fleet == f])
      if (length(z) != 1L) cli::cli_abort("F correlation groups vary within fleet {f}; supply one group per fleet.")
      z
    }, character(1))
    counts <- if (dat$has_fleet) table(dat$F_states$fleet)[dat$fleets] else length(dat$F_settings$age_block_start)
    applicable <- if (dimension == "age") counts > 1L else rep(sum(!dat$is_proj) > 1L, length(dat$fleets))
    levels <- unique(unname(groups[applicable]))
    list(group = groups, levels = levels, index = match(groups, levels))
  })
  names(dat$F_cor) <- c("age", "year")
  dat
}

.f_cor_parameters <- function(dat, dimension) {
  groups <- dat$F_cor[[dimension]]$levels
  if (!length(groups)) return(NULL)
  stats::setNames(numeric(length(groups)), if (!dat$has_fleet || (length(groups) == 1L && identical(groups, "shared"))) dimension else paste(dimension, groups, sep = " / "))
}

.fleet_phi <- function(par, dat, f) {
  "[<-" <- RTMB::ADoverload("[<-")
  values <- c(0, 0)
  for (i in seq_along(values)) {
    d <- c("age", "year")[i]
    spec <- dat$F_cor[[d]]; j <- spec$index[f]
    if (is.na(j)) next
    key <- names(.f_cor_parameters(dat, d))[j]
    values[i] <- plogis(par$logit_phi_f[[key]])
  }
  values
}

# Independent densities per fleet, with the existing innovation normalization.
.fleet_F_process <- function(states, mean, par, dat, simulate = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")
  type <- dat$F_settings$process
  if (type == "off") return(if (simulate) states else 0)
  result <- 0; sd <- .process_sd(par, dat, "F")
  for (f in seq_along(dat$fleets)) {
    columns <- which(dat$F_states$fleet == dat$fleets[f])
    mu <- .f_state_mean(mean[!dat$is_proj, , drop = FALSE], dat)[, columns, drop = FALSE]
    x <- states[, columns, drop = FALSE] - mu
    scale <- if (length(sd) == 1L) sd else sd[columns]
    phi <- if (type == "ar1") .fleet_phi(par, dat, f) else c(0, 0)
    rho <- if (type == "cor_rw" && length(columns) > 1L) {
      tanh(par$atanh_rho_f[dat$F_cor$age$index[f]])
    } else 0
    if (simulate) {
      draw <- if (!dat$process_sd$F$default && type != "cor_rw") {
        .rprocess_scaled(x, scale, type, phi)
      } else switch(type, rw = rprocess_rw(x, scale),
        cor_rw = rprocess_cor_rw(x, scale, rho),
        iid = matrix(stats::rnorm(length(x), 0, scale), nrow(x)),
        ar1 = rprocess_ar1(nrow(x), ncol(x), sd = scale, phi = phi))
      states[, columns] <- mu + draw
    } else {
      result <- result + if (!dat$process_sd$F$default && type != "cor_rw") {
        .dprocess_scaled(x, scale, type, phi)
      } else switch(type, rw = dprocess_rw(x, scale),
        cor_rw = dprocess_cor_rw(x, scale, rho),
        iid = sum(RTMB::dnorm(x, 0, scale, log = TRUE)),
        ar1 = dprocess_ar1(x, sd = scale, phi = phi))
    }
  }
  if (simulate) states else result
}

.parse_fleet_sd <- function(dat, formula, default) {
  d <- .f_data(dat, formula, "F_settings$sd_form")
  keep <- !dat$is_proj[match(d$year, dat$years)]
  frame <- stats::model.frame(formula, d[keep, ], na.action = stats::na.pass)
  design <- stats::model.matrix(formula, frame)
  if (nrow(design) != sum(keep) || !ncol(design) || any(!is.finite(design))) {
    cli::cli_abort("F_settings$sd_form needs finite covariates at every active catch-grid row.")
  }
  surface <- d[keep, ]; out <- list()
  for (f in dat$fleets) {
    blocks <- dat$F_blocks[[f]]
    for (b in levels(blocks)) {
      ages <- as.numeric(names(blocks)[blocks == b])
      rows <- design[surface$fleet == f & surface$age %in% ages, , drop = FALSE]
      if (any(rows != rep(rows[1L, ], each = nrow(rows)))) {
        cli::cli_abort("F_settings$sd_form varies within fleet {f}, F age block {b} (ages {paste(ages, collapse = ', ')}) or across years. Share SD covariates or refine age_breaks.")
      }
      out[[length(out) + 1L]] <- rows[1L, , drop = FALSE]
    }
  }
  design <- do.call(rbind, out); rownames(design) <- dat$F_states$key
  if (qr(design)$rank < ncol(design)) cli::cli_abort("F_settings$sd_form is rank deficient on the active fleet states.")
  list(matrix = design, default = default,
    parameter = if (default) "log_sd_f" else "sd_beta_f")
}
