# Linear recurrences work on numeric and RTMB AD values alike.
.process_whiten <- function(x, sd, process, phi = c(0, 0), rho = 0) {
  if (process %in% c("rw", "cor_rw")) x <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  z <- sweep(x, 2L, rep(sd, length.out = ncol(x)), "/")
  if (process == "ar1") {
    y <- z
    if (nrow(z) > 1L) y[-1L, ] <- z[-1L, , drop = FALSE] - phi[2L] * z[-nrow(z), , drop = FALSE]
    y[1L, ] <- z[1L, ] * sqrt(1 - phi[2L]^2)
    z <- y
    if (ncol(y) > 1L) z[, -1L] <- y[, -1L, drop = FALSE] - phi[1L] * y[, -ncol(y), drop = FALSE]
    z[, 1L] <- y[, 1L] * sqrt(1 - phi[1L]^2)
  }
  if (process == "cor_rw" && ncol(z) > 1L) {
    y <- z
    z[, -1L] <- (y[, -1L, drop = FALSE] - rho * y[, -ncol(y), drop = FALSE]) / sqrt(1 - rho^2)
  }
  z
}

.process_innovations <- function(par, dat, population, log_F, log_mu_F, log_mu_M) {
  out <- list()
  historical <- as.character(dat$years[!dat$is_proj])
  rec <- dat$rec
  if (rec$type != "off" && length(rec$eligible)) {
    i <- rec$eligible
    u <- population$log_recruitment - population$log_mu_R
    scale <- .rec_scale(par, dat)
    z <- if (rec$type == "rw") u[i] - u[i - 1L] else u[i]
    if (rec$type == "ar1") z <- c(u[i[1L]] * sqrt(1 - scale$phi^2),
      if (length(i) > 1L) u[i[-1L]] - scale$phi * u[i[-length(i)]])
    out$recruitment <- z[!dat$is_proj[i]] / scale$sd
  }
  states <- list(
    N = if (!is.null(par[["log_n"]])) population$log_N[-1L, -1L, drop = FALSE] - population$pred_log_N[-1L, -1L, drop = FALSE],
    F = if (!is.null(par[["log_f"]])) par$log_f - .f_state_mean(log_mu_F[!dat$is_proj, , drop = FALSE], dat),
    M = if (!is.null(par[["log_m"]])) par$log_m - log_mu_M[rownames(par$log_m), dat$M_settings$age_block_start, drop = FALSE]
  )
  for (component in names(states)) {
    x <- states[[component]]
    process <- dat[[paste0(component, "_settings")]]$process
    if (is.null(x) || process == "off" || (process %in% c("rw", "cor_rw") && nrow(x) < 2L)) next
    suffix <- tolower(component)
    z <- .process_whiten(x, .process_sd(par, dat, component), process,
      phi = if (process == "ar1") {
        if (component == "F") .f_phi(par) else plogis(par[[paste0("logit_phi_", suffix)]])
      } else c(0, 0),
      rho = if (process == "cor_rw") .f_rho(par) else 0)
    out[[component]] <- z[rownames(z) %in% historical, , drop = FALSE]
  }
  do.call(c, lapply(out, as.vector))
}

.process_layout <- function(par, dat) {
  out <- list()
  if (dat$rec$type != "off") {
    i <- dat$rec$eligible
    out$recruitment <- data.frame(year = dat$years[i[!dat$is_proj[i]]])
  }
  for (component in c("N", "F", "M")) {
    x <- par[[paste0("log_", tolower(component))]]
    process <- dat[[paste0(component, "_settings")]]$process
    if (is.null(x) || process == "off") next
    if (process %in% c("rw", "cor_rw")) x <- x[-1L, , drop = FALSE]
    x <- x[rownames(x) %in% as.character(dat$years[!dat$is_proj]), , drop = FALSE]
    d <- expand.grid(year = as.numeric(rownames(x)), coordinate = colnames(x), stringsAsFactors = FALSE)
    if (component == "M" || (component == "F" && .f_blocked(dat))) {
      names(d)[2L] <- "age_block"
      d$state_block <- match(d$age_block, colnames(x))
    } else {
      d$age <- as.numeric(d$coordinate)
      d$coordinate <- NULL
    }
    out[[component]] <- d
  }
  Filter(function(d) nrow(d) > 0L, out)
}

.process_method <- function() "gaussian-joint-delta-v1"

.validate_process_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    cli::cli_abort("{.arg process_seed} must be a single non-negative integer.")
  }
  as.integer(seed)
}

# Q = P' L L' P; delta = P' L'^(-1) z has covariance Q^(-1).
.process_perturbation <- function(precision, seed) {
  # Matrix caches Cholesky factors in its input. Give it a fresh sparse object
  # so the fitted object's stored uncertainty remains byte-for-byte unchanged.
  precision <- precision * 1
  precision <- methods::as(Matrix::forceSymmetric(precision), "CsparseMatrix")
  factor <- withCallingHandlers(Matrix::Cholesky(precision, LDL = FALSE, super = FALSE),
    warning = function(w) invokeRestart("muffleWarning"))
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  z <- stats::rnorm(nrow(precision))
  as.numeric(Matrix::solve(factor, Matrix::solve(factor, z, system = "Lt"), system = "Pt"))
}

.tidy_process_residuals <- function(fit, seed, enabled = TRUE) {
  status <- function(code, detail) list(status = data.frame(status = code, detail = detail,
    seed = seed, method = .process_method(), stringsAsFactors = FALSE))
  if (!enabled) return(status("disabled", "Process residual calculation was disabled."))
  if (inherits(fit, "tam_ref")) return(status("reference", "Assessment references do not contain a tinyAM joint distribution."))
  active <- c(fit$dat$rec$type, fit$dat$N_settings$process, fit$dat$F_settings$process, fit$dat$M_settings$process)
  if (!length(active) || all(active == "off")) return(status("deterministic", "No stochastic historical processes."))
  tryCatch({
    mode <- fit$parameter_values
    Q <- fit$sdrep$jointPrecision
    if (!length(mode) || any(!is.finite(mode)) || is.null(Q) ||
        !identical(dim(Q), c(length(mode), length(mode))) ||
        !identical(names(mode), rownames(Q)) || !identical(rownames(Q), colnames(Q)) ||
        any(!is.finite(Q@x)) || !isTRUE(Matrix::isSymmetric(Q))) {
      cli::cli_abort("Stored joint precision or parameter ordering is unavailable or inconsistent.")
    }
    if (identical(fit$sdrep$pdHess, FALSE)) cli::cli_abort("The stored uncertainty calculation has a non-positive-definite Hessian.")
    random <- fit$obj$env$lrandom()
    if (length(random) != length(mode)) cli::cli_abort("Stored random/fixed parameter ordering is inconsistent.")
    for (part in c("fixed", "random")) {
      stored <- fit$sdrep[[paste0("par.", part)]]
      if (!is.null(stored) && !isTRUE(all.equal(stored, mode[if (part == "random") random else !random]))) {
        cli::cli_abort("Joint precision was evaluated at a different parameter mode or ordering.")
      }
    }
    par <- fit$obj$env$parList(par = mode)
    dat <- fit$dat
    layout <- .process_layout(par, dat)
    if (!length(layout)) return(status("deterministic", "No penalized historical process contributions."))
    dat$.likelihood_component <- "process_innovations"
    obj <- RTMB::MakeADFun(function(p) {
      innovations <- nll_fun(p, dat)
      RTMB::ADREPORT(innovations)
      0
    }, par, map = fit$parameter_map, ADreport = TRUE, silent = TRUE)
    if (!identical(names(obj$par), names(mode)) || !isTRUE(all.equal(unname(obj$par), unname(mode)))) {
      cli::cli_abort("Mapped diagnostic parameter ordering differs from the fitted mode.")
    }
    tape <- RTMB::GetTape(obj)
    fitted <- as.numeric(tape(mode))
    J <- methods::as(tape$jacfun(sparse = TRUE)(mode), "CsparseMatrix")
    if (length(fitted) != sum(vapply(layout, nrow, integer(1))) ||
        any(!is.finite(fitted)) || any(!is.finite(J@x))) cli::cli_abort("Process innovation values or derivatives are invalid.")
    delta <- .process_perturbation(Q, seed)
    sample <- fitted + as.numeric(J %*% delta)
    if (any(!is.finite(sample))) cli::cli_abort("Joint-sample residuals are non-finite.")
    offset <- 0L
    tables <- lapply(layout, function(d) {
      i <- seq_len(nrow(d)) + offset
      offset <<- offset + nrow(d)
      d$est <- sample[i]
      d$fitted <- fitted[i]
      d$is_proj <- FALSE
      d
    })
    c(tables, status("ok", "One Gaussian joint perturbation propagated through the AD innovation Jacobian."))
  }, error = function(e) {
    detail <- conditionMessage(e)
    cli::cli_warn("Joint-sample process residuals unavailable: {detail}")
    status("unavailable", detail)
  })
}
