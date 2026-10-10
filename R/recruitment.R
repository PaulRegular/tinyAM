# Recruitment residuals act on absolute log recruitment, not additional states.
.parse_recruitment <- function(dat) {
  form <- dat$N_settings$rec_form
  if (is.null(form)) form <- ~ rw(year)
  if (!inherits(form, "formula") || length(form) != 2L) {
    cli::cli_abort("N_settings$rec_form must be a one-sided formula.")
  }
  process <- list()
  strip <- function(x) {
    if (!is.call(x)) return(x)
    type <- .formula_call_name(x)
    if (type %in% c("iid", "rw", "ar1")) {
      fun <- switch(type, iid = iid, rw = rw, ar1 = ar1)
      x[[1L]] <- as.name(type)
      args <- tryCatch(as.list(match.call(fun, x))[-1L], error = function(e)
        cli::cli_abort("Invalid recruitment {type}() arguments."))
      if (!identical(args$x, quote(year)) ||
          (!is.null(args$by) && !identical(args$by, quote(NULL)))) {
        cli::cli_abort("Recruitment residuals must use year without by grouping.")
      }
      fixed <- function(x, name, upper = Inf, zero = FALSE) {
        if (is.null(x) || identical(x, quote(NULL))) return(NULL)
        value <- tryCatch(eval(x, environment(form)), error = function(e) NULL)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0 || (!zero && value == 0) || value >= upper) {
          cli::cli_abort("Recruitment {name} must be one finite {if (zero) 'non-negative' else 'positive'} number below {upper}.")
        }
        value
      }
      process[[length(process) + 1L]] <<- list(type = type,
        sd = fixed(args$sd, "sd"), phi = fixed(args$phi, "phi", 1, TRUE))
      return(NULL)
    }
    if (type %in% c("bh", "ricker", "mono", "logistic", "|")) {
      cli::cli_abort("Unsupported recruitment formula term: {type}().")
    }
    if (type == "(") return(strip(x[[2L]]))
    contains <- function(z) is.call(z) &&
      (.formula_call_name(z) %in% c("iid", "rw", "ar1") ||
       any(vapply(as.list(z)[-1L], contains, logical(1))))
    if (type == "+") {
      parts <- Filter(Negate(is.null), lapply(as.list(x)[-1L], strip))
      if (!length(parts)) return(NULL)
      return(Reduce(function(a, b) call("+", a, b), parts))
    }
    if (type == "-" && length(x) == 3L && !contains(x[[3L]])) {
      left <- strip(x[[2L]])
      return(call("-", if (is.null(left)) 1 else left, x[[3L]]))
    }
    if (contains(x)) cli::cli_abort("Recruitment process terms must be additive.")
    x
  }
  ordinary <- strip(form[[2L]])
  if (length(process) != 1L) cli::cli_abort("rec_form requires exactly one iid(year), rw(year), or ar1(year) residual process.")
  rec <- process[[1L]]
  fixed_form <- form
  fixed_form[[2L]] <- if (is.null(ordinary)) 1 else ordinary
  data <- dat$obs$maturity[dat$obs$maturity$age == min(dat$ages), , drop = FALSE]
  data <- data[match(dat$years, data$year), , drop = FALSE]
  variables <- all.vars(fixed_form)
  if (any(!variables %in% names(data))) cli::cli_abort("Recruitment covariates must be columns of obs$maturity at the youngest modeled age.")
  frame <- stats::model.frame(fixed_form, data, na.action = stats::na.pass)
  matrix <- stats::model.matrix(fixed_form, frame)
  if (nrow(matrix) != length(dat$years) || any(!is.finite(matrix))) {
    cli::cli_abort("rec_form needs finite covariate values in every modeled recruitment year.")
  }
  if (rec$type == "rw") {
    constant <- vapply(seq_len(ncol(matrix)), function(j) all(diff(matrix[, j]) == 0), logical(1))
    removed <- colnames(matrix)[constant & colnames(matrix) != "(Intercept)"]
    if (length(removed)) cli::cli_warn("Constant recruitment RW covariates cancel from increments and are removed: {paste(removed, collapse = ', ')}.")
    matrix <- matrix[, !constant, drop = FALSE]
  }
  eligible <- seq.int(2L, length(dat$years))
  historical <- eligible[!dat$is_proj[eligible]]
  n <- length(historical)
  if (rec$type == "iid" && is.null(rec$sd) && n < 2L) cli::cli_abort("Estimating recruitment IID SD needs at least two process years.")
  if (rec$type == "ar1" && is.null(rec$phi) && n < 3L) cli::cli_abort("Estimating recruitment AR1 correlation needs at least three process years; supply phi.")
  design <- if (rec$type == "rw") matrix[historical, , drop = FALSE] - matrix[historical - 1L, , drop = FALSE] else matrix[historical, , drop = FALSE]
  if (ncol(design) && (qr(design)$rank < ncol(design) ||
      (is.null(rec$sd) && ncol(design) >= n))) {
    cli::cli_abort("Recruitment fixed effects are redundant or saturated; simplify rec_form.")
  }
  held <- setdiff(variables, "year")
  if (any(dat$is_proj) && length(held)) cli::cli_inform("Recruitment projections hold terminal covariate values: {paste(held, collapse = ', ')}.")
  dat$N_settings$rec_form <- form
  dat$rec <- c(rec, list(matrix = matrix, eligible = eligible, historical = historical,
                       boundary = 1L, data = data, covariates = held))
  dat
}

.rec_mean <- function(par, dat) {
  if (!ncol(dat$rec$matrix)) return(rep(0, length(dat$years)))
  drop(dat$rec$matrix %*% par$rec_beta)
}

.rec_scale <- function(par, dat) {
  list(sd = if (is.null(dat$rec$sd)) exp(par$log_sd_r) else dat$rec$sd,
       phi = if (dat$rec$type == "ar1") {
         if (is.null(dat$rec$phi)) stats::plogis(par$logit_phi_r) else dat$rec$phi
       } else 0)
}

.rec_nll <- function(log_r, mean, par, dat) {
  i <- dat$rec$eligible
  scale <- .rec_scale(par, dat)
  u <- log_r - mean
  if (dat$rec$type == "rw") return(-sum(RTMB::dnorm(u[i] - u[i - 1L], 0, scale$sd, log = TRUE)))
  if (dat$rec$type == "iid") return(-sum(RTMB::dnorm(u[i], 0, scale$sd, log = TRUE)))
  nll <- -RTMB::dnorm(u[i[1L]], 0, scale$sd / sqrt(1 - scale$phi^2), log = TRUE)
  if (length(i) > 1L) nll <- nll - sum(RTMB::dnorm(u[i[-1L]], scale$phi * u[i[-length(i)]], scale$sd, log = TRUE))
  nll
}

.simulate_rec <- function(par, dat) {
  mean <- .rec_mean(par, dat)
  scale <- .rec_scale(par, dat)
  u <- numeric(length(dat$years))
  u[1L] <- par$log_r0 - mean[1L]
  for (i in dat$rec$eligible) {
    expected <- switch(dat$rec$type, rw = u[i - 1L], iid = 0,
      ar1 = if (i == dat$rec$eligible[1L]) 0 else scale$phi * u[i - 1L])
    sd <- if (dat$rec$type == "ar1" && i == dat$rec$eligible[1L]) scale$sd / sqrt(1 - scale$phi^2) else scale$sd
    u[i] <- stats::rnorm(1L, expected, sd)
  }
  (mean + u)[-1L]
}
