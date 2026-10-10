# Recruitment residuals act on absolute log recruitment, not additional states.
.parse_recruitment <- function(dat) {
  form <- dat$N_settings$rec_form
  if (is.null(form)) form <- ~ rw(year)
  if (!inherits(form, "formula") || length(form) != 2L) {
    cli::cli_abort("N_settings$rec_form must be a one-sided formula.")
  }
  process <- list()
  curves <- list()
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
    if (type %in% c("bh", "ricker")) {
      x[[1L]] <- as.name(type)
      args <- tryCatch(as.list(match.call(bh, x))[-1L], error = function(e)
        cli::cli_abort("Invalid {type}() arguments."))
      if (!identical(args$ssb, quote(ssb))) cli::cli_abort("{type}() must use the model's ssb.")
      lag <- if (is.null(args$lag) || identical(args$lag, quote(NULL))) min(dat$ages) else
        tryCatch(eval(args$lag, environment(form)), error = function(e) NULL)
      if (!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag < 0 || lag != floor(lag)) {
        cli::cli_abort("Stock-recruit lag must be one non-negative integer.")
      }
      curves[[length(curves) + 1L]] <<- list(type = type, lag = as.integer(lag))
      return(NULL)
    }
    if (type %in% c("mono", "logistic", "|")) {
      cli::cli_abort("Unsupported recruitment formula term: {type}().")
    }
    if (type == "(") return(strip(x[[2L]]))
    contains <- function(z) is.call(z) &&
      (.formula_call_name(z) %in% c("iid", "rw", "ar1", "bh", "ricker") ||
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
  if (length(curves) > 1L) cli::cli_abort("rec_form allows one stock-recruit curve.")
  curve <- if (length(curves)) curves[[1L]] else NULL
  if (!is.null(curve) && rec$type == "rw") cli::cli_abort("Stock-recruit curves currently support IID or AR1 residuals, not RW residuals.")
  fixed_form <- form
  fixed_form[[2L]] <- if (is.null(ordinary)) 1 else ordinary
  data <- dat$obs$maturity[dat$obs$maturity$age == min(dat$ages), , drop = FALSE]
  data <- data[match(dat$years, data$year), , drop = FALSE]
  variables <- all.vars(fixed_form)
  if ("ssb" %in% variables) cli::cli_abort("Use modeled ssb only inside bh() or ricker(); supply other annual covariates explicitly.")
  if (any(!variables %in% names(data))) cli::cli_abort("Recruitment covariates must be columns of obs$maturity at the youngest modeled age.")
  frame <- stats::model.frame(fixed_form, data, na.action = stats::na.pass)
  matrix <- stats::model.matrix(fixed_form, frame)
  if (nrow(matrix) != length(dat$years) || any(!is.finite(matrix))) {
    cli::cli_abort("rec_form needs finite covariate values in every modeled recruitment year.")
  }
  if (!is.null(curve)) matrix <- matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  if (rec$type == "rw") {
    constant <- vapply(seq_len(ncol(matrix)), function(j) all(diff(matrix[, j]) == 0), logical(1))
    removed <- colnames(matrix)[constant & colnames(matrix) != "(Intercept)"]
    if (length(removed)) cli::cli_warn("Constant recruitment RW covariates cancel from increments and are removed: {paste(removed, collapse = ', ')}.")
    matrix <- matrix[, !constant, drop = FALSE]
  }
  first <- if (is.null(curve)) 2L else max(2L, curve$lag + 1L)
  if (first > sum(!dat$is_proj)) cli::cli_abort("No historical parent-recruit pairs are available; shorten lag or extend modeled years.")
  eligible <- seq.int(first, length(dat$years))
  historical <- eligible[!dat$is_proj[eligible]]
  n <- length(historical)
  if (rec$type == "iid" && is.null(rec$sd) && n < 2L) cli::cli_abort("Estimating recruitment IID SD needs at least two process years.")
  if (rec$type == "ar1" && is.null(rec$phi) && n < 3L) cli::cli_abort("Estimating recruitment AR1 correlation needs at least three process years; supply phi.")
  if (!is.null(curve)) {
    if (curve$lag == 0L && any(dat$W[, 1L] * dat$P[, 1L] != 0)) {
      cli::cli_abort("Same-year SSB is circular when recruits contribute mature biomass. Use a positive lag.")
    }
    parent <- eligible - curve$lag
    possible <- rowSums(dat$W * dat$P) > 0
    if (!is.null(dat$plus_ages)) possible <- possible | rowSums(dat$W_plus_input * dat$P_plus_input) > 0
    if (any(!possible[parent])) cli::cli_abort("Parent SSB is zero in some stock-recruit years; supply valid weight and maturity or revise lag/years.")
  }
  design <- if (rec$type == "rw") matrix[historical, , drop = FALSE] - matrix[historical - 1L, , drop = FALSE] else matrix[historical, , drop = FALSE]
  if (!is.null(curve) && ncol(design) && qr(cbind(1, design))$rank < ncol(design) + 1L) {
    cli::cli_abort("Recruitment covariates duplicate the stock-recruit productivity baseline; remove constant or redundant columns.")
  }
  if (ncol(design) && (qr(design)$rank < ncol(design) ||
      (is.null(rec$sd) && ncol(design) >= n))) {
    cli::cli_abort("Recruitment fixed effects are redundant or saturated; simplify rec_form.")
  }
  held <- setdiff(variables, "year")
  if (any(dat$is_proj) && length(held)) cli::cli_inform("Recruitment projections hold terminal covariate values: {paste(held, collapse = ', ')}.")
  dat$N_settings$rec_form <- form
  dat$rec <- c(rec, list(matrix = matrix, eligible = eligible, historical = historical,
                       boundary = seq_len(first - 1L), data = data, covariates = held,
                       curve = curve, fixed_form = fixed_form))
  dat
}

#' Recruitment formulas and stock-recruit relationships
#'
#' @description
#' Recruitment is the number of fish entering the youngest modeled age.
#' Use `N_settings$rec_form` to describe persistent changes, independent good and
#' bad years, environmental effects, or a relationship with spawning biomass.
#' The default `~ rw(year)` retains tinyAM's original recruitment random walk.
#'
#' @details
#' Use exactly one [iid()], [rw()], or [ar1()] residual term on `year`, without
#' `by` grouping. Ordinary covariates are read from `obs$maturity` at the youngest
#' modeled age; their values at older ages may be `NA`. Covariates refer to the
#' recruitment year. Construct lagged or standardized columns explicitly.
#' Projections hold their terminal values and announce that assumption.
#'
#' `bh(ssb)` rises toward a plateau; `ricker(ssb)` can decline at high spawning
#' biomass. Either curve supports IID or AR1 residuals. These curves describe
#' median recruitment, not the arithmetic mean, and use the model's existing
#' start-of-year SSB and supplied weight/maturity units.
#' With parent biomass \eqn{S_{t-L}}, their definitions are
#' \deqn{g_{BH}(S)=\alpha S/(1+\beta S),\qquad
#' g_{Ricker}(S)=\alpha S\exp(-\beta S),\quad \alpha,\beta>0.}
#' \deqn{\log R_t=\log g(S_{t-L})+X_t\gamma+u_t.}
#' The curve supplies its baseline, so the ordinary intercept is omitted.
#' IID residuals have variance \eqn{\sigma_R^2}. AR1 residuals obey
#' \eqn{u_t=\phi u_{t-1}+\epsilon_t}, with innovation variance
#' \eqn{\sigma_R^2} and stationary first variance
#' \eqn{\sigma_R^2/(1-\phi^2)}. No lognormal mean correction is applied.
#' Alpha has recruitment-number/SSB units; beta has inverse-SSB units.
#' Both are optimized on the log scale, independently of equilibrium assumptions.
#'
#' Without a curve, IID/AR1 use \eqn{\log R_t=X_t\gamma+u_t}.
#' A RW instead obeys
#' \deqn{\log R_t-\log R_{t-1}=(X_t-X_{t-1})\gamma+\epsilon_t.}
#' Constant columns cancel and are removed. SD is estimated unless supplied;
#' estimated AR1 correlation is positive, or a fixed `phi` in `[0, 1)` can be
#' supplied. See [formula_effects] for these process conventions.
#'
#' First-year `log_r0` is always a free fixed anchor. Any additional recruitment
#' years whose parent SSB precedes the modeled period are also free fixed states
#' (`log_r_init`), without process penalties. Later `log_r` values are absolute
#' random log states; residuals are calculated internally, without an additional
#' annual random-effect block. The relationship does not reconstruct pre-model
#' SSB. AR1 starts at the first eligible year with its stationary density.
#' A zero lag is allowed only when recruits contribute no mature biomass.
#' Changing recruitment age changes the default lag and boundary states.
#'
#' A converged fit does not establish that a curve's shape is well estimated.
#' Limited SSB variation, correlated covariates, and other population processes
#' can obscure the relationship. Inspect [check_tam()] and recovery simulations.
#'
#' @param ssb The literal `ssb`, referring to modeled spawning stock biomass.
#' @param lag Non-negative integer parent-year lag. `NULL` uses the youngest
#'   modeled age. SSB is start-of-year, not adjusted to spawning time.
#' @return Formula markers; direct calls raise an error.
#' @examples
#' ~ rw(year)
#' ~ iid(year)
#' ~ bh(ssb) + iid(year)
#' ~ ricker(ssb, lag = 2) + ar1(year)
#' # Add temperature to obs$maturity at the youngest modeled age, then use:
#' ~ bh(ssb) + temperature + iid(year)
#' @seealso [prepare_tam()], [fit_tam()], [sim_tam()], [check_tam()]
#' @name recruitment_formulas
#' @export
bh <- function(ssb, lag = NULL) {
  cli::cli_abort("bh() is a recruitment formula marker; use it in N_settings$rec_form.")
}

#' @rdname recruitment_formulas
#' @export
ricker <- function(ssb, lag = NULL) {
  cli::cli_abort("ricker() is a recruitment formula marker; use it in N_settings$rec_form.")
}

.rec_log_curve <- function(log_ssb, par, type) {
  density <- par$log_sr_beta + log_ssb
  par$log_sr_alpha + log_ssb - if (type == "bh") RTMB::logspace_add(0, density) else exp(density)
}

.rec_mean <- function(par, dat) {
  if (!ncol(dat$rec$matrix)) return(rep(0, length(dat$years)))
  drop(dat$rec$matrix %*% par$rec_beta)
}

.rec_scale <- function(par, dat) {
  list(sd = if (is.null(dat$rec$sd)) exp(par$log_sd_r) else dat$rec$sd,
       phi = if (dat$rec$type == "ar1") {
         if (is.null(dat$rec$phi)) plogis(par$logit_phi_r) else dat$rec$phi
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

# One chronological population construction for estimation and simulation.
.population_states <- function(par, dat, Z, simulate = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")
  T <- length(dat$years)
  A <- length(dat$ages)
  older <- 2:A
  log_N <- pred_log_N <- matrix(0, T, A, dimnames = list(year = dat$years, age = dat$ages))
  log_recruitment <- numeric(T)
  log_recruitment[1L] <- par$log_r0
  if (length(dat$rec$boundary) > 1L) log_recruitment[dat$rec$boundary[-1L]] <- par$log_r_init
  log_recruitment[dat$rec$eligible] <- par$log_r
  log_N[, 1L] <- log_recruitment
  if (dat$N_settings$process != "off") log_N[-1L, -1L] <- par$log_n
  initial_error <- numeric(A - 1L)
  if (simulate && dat$N_settings$init == "random") initial_error[] <- stats::rnorm(A - 1L, 0, exp(par$log_sd_n0))
  for (a in older) {
    pred_log_N[1L, a] <- log_N[1L, a - 1L] - Z[1L, a - 1L]
    log_N[1L, a] <- if (dat$N_settings$init == "exp" || (simulate && dat$N_settings$init == "random")) {
      pred_log_N[1L, a] + initial_error[a - 1L]
    } else par$log_n0[a - 1L]
  }
  eta_log_n0 <- log_N[1L, -1L] - pred_log_N[1L, -1L]
  eta_log_N <- matrix(0, T - 1L, A - 1L)
  if (simulate && dat$N_settings$process == "iid") eta_log_N[] <- stats::rnorm(length(eta_log_N), 0, exp(par$log_sd_n))
  if (simulate && dat$N_settings$process == "ar1") eta_log_N <- rprocess_ar1(T - 1L, A - 1L, sd = exp(par$log_sd_n), phi = stats::plogis(par$logit_phi_n))
  W <- dat$W
  P <- dat$P
  log_ssb <- numeric(T)
  mean <- .rec_mean(par, dat)
  log_mu_R <- log_pred_R <- mean
  u <- numeric(T)
  scale <- .rec_scale(par, dat)
  if (!is.null(dat$plus_ages)) {
    B <- length(dat$plus_ages)
    log_N_plus <- matrix(0, T, B, dimnames = list(year = dat$years, age = dat$plus_ages))
  }
  for (y in seq_len(T)) {
    if (y > 1L) {
      pred_log_N[y, older] <- log_N[y - 1L, older - 1L] - Z[y - 1L, older - 1L]
      pred_log_N[y, A] <- RTMB::logspace_add(pred_log_N[y, A], log_N[y - 1L, A] - Z[y - 1L, A])
      if (simulate && dat$N_settings$process == "rw" && y == 2L) {
        eta_log_N[1L, ] <- par$log_n[1L, ] - pred_log_N[y, older]
        eta_log_N <- rprocess_rw(eta_log_N, sd = exp(par$log_sd_n))
      }
      if (dat$N_settings$process == "off" || simulate) log_N[y, older] <- pred_log_N[y, older] + eta_log_N[y - 1L, ]
    }
    if (!is.null(dat$plus_ages)) {
      if (y == 1L) {
        components <- -(seq_len(B) - 1L) * Z[1L, A]
        components[-B] <- components[-B] + log(-expm1(-Z[1L, A]))
      } else {
        components <- c(log_N[y - 1L, A - 1L] - Z[y - 1L, A - 1L], log_N_plus[y - 1L, -B] - Z[y - 1L, A])
        components[B] <- RTMB::logspace_add(components[B], log_N_plus[y - 1L, B] - Z[y - 1L, A])
      }
      shares <- components - Reduce(RTMB::logspace_add, components)
      log_N_plus[y, ] <- log_N[y, A] + shares
      W[y, A] <- sum(exp(shares) * dat$W_plus_input[y, ])
      P[y, A] <- if (all(dat$W_plus_input[y, ] == 0)) 0 else
        sum(exp(shares) * dat$W_plus_input[y, ] * dat$P_plus_input[y, ]) / W[y, A]
    }
    if (y %in% dat$rec$eligible) {
      if (!is.null(dat$rec$curve)) {
        parent <- y - dat$rec$curve$lag
        parent_ssb <- if (parent == y) {
          # Zero-lag validation guarantees no mature biomass in recruitment age.
          log(sum(exp(log_N[y, older]) * W[y, older] * P[y, older]))
        } else log_ssb[parent]
        log_mu_R[y] <- mean[y] + .rec_log_curve(parent_ssb, par, dat$rec$curve$type)
      }
      log_pred_R[y] <- switch(dat$rec$type,
        rw = log_recruitment[y - 1L] + mean[y] - mean[y - 1L],
        iid = log_mu_R[y],
        ar1 = log_mu_R[y] + if (y == dat$rec$eligible[1L]) 0 else scale$phi * u[y - 1L])
      if (simulate) {
        sd <- if (dat$rec$type == "ar1" && y == dat$rec$eligible[1L]) scale$sd / sqrt(1 - scale$phi^2) else scale$sd
        log_recruitment[y] <- stats::rnorm(1L, log_pred_R[y], sd)
        log_N[y, 1L] <- log_recruitment[y]
      }
      u[y] <- log_recruitment[y] - log_mu_R[y]
    }
    log_ssb[y] <- log(sum(exp(log_N[y, ]) * W[y, ] * P[y, ]))
  }
  names(log_recruitment) <- dat$years
  list(log_N = log_N, pred_log_N = pred_log_N, log_recruitment = log_recruitment,
    log_mu_R = log_mu_R, log_pred_R = log_pred_R,
    eta_R = log_recruitment - log_pred_R, eta_R_state = u,
    eta_log_n0 = eta_log_n0, W = W, P = P,
    N_plus = if (!is.null(dat$plus_ages)) exp(log_N_plus) else NULL)
}

.rec_process_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  rec <- dat$rec
  if (!is.null(rec) && (rec$type != "rw" || ncol(rec$matrix) || !is.null(rec$curve)) &&
      length(rec$historical) < 5L && (is.null(rec$sd) || (rec$type == "ar1" && is.null(rec$phi)))) {
    out <- rbind(out, data.frame(issue = "recruitment_support", detail =
      "Fewer than five historical recruitment process years inform estimated SD/correlation. Consider a simpler formula or externally supported sd/phi."))
  }
  out
}

.initialize_rec_curve <- function(par, dat) {
  mean_M <- matrix(dat$log_mu_supplied_m + drop(dat$M_modmat %*% if (is.null(par$mu_m)) dat$mu_m else par$mu_m),
                   length(dat$years), length(dat$ages), dimnames = list(dat$years, dat$ages))
  M <- exp(mean_M)
  if (!is.null(par[["log_m"]])) M[rownames(par[["log_m"]]), names(dat$M_settings$age_blocks)] <-
    exp(par[["log_m"]][, dat$M_settings$age_blocks, drop = FALSE])
  F <- matrix(0, length(dat$years), length(dat$ages))
  F[!dat$is_proj, ] <- exp(par$log_f)
  if (any(dat$is_proj)) F[dat$is_proj, ] <- sweep(exp(par$log_f[rep(nrow(par$log_f), sum(dat$is_proj)), , drop = FALSE]),
                                                  1L, dat$proj_settings$F_mult, `*`)
  initial <- .population_states(par, dat, F + M)
  i <- dat$rec$eligible[1L] - dat$rec$curve$lag
  S <- sum(exp(initial$log_N[i, ]) * initial$W[i, ] * initial$P[i, ])
  par$log_sr_beta <- -log(S)
  par$log_sr_alpha <- par$log_r0 - log(S) + if (dat$rec$curve$type == "bh") log(2) else 1
  par
}

.rec_fit_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  rec <- fit$dat$rec
  if (is.null(rec)) return(out)
  if (!is.null(rec$curve)) {
    S <- exp(fit$rep$rec_log_parent[!fit$dat$is_proj[rec$eligible]])
    if (length(S) && all(is.finite(S)) && max(S) / min(S) < 2) {
      out <- rbind(out, data.frame(issue = "stock_recruit_support", detail =
        "Fitted parent SSB varies by less than a factor of two. Stock-recruit curve shape may be weakly supported; examine uncertainty and covariate sensitivity."))
    }
  }
  sdr <- fit[["sdrep"]]
  if (!is.list(sdr) || !isTRUE(sdr$pdHess)) return(out)
  p <- sdr$par.fixed
  se <- sqrt(pmax(0, diag(sdr$cov.fixed)))
  if (!is.null(rec$curve)) {
    i <- which(names(p) %in% c("log_sr_alpha", "log_sr_beta"))
    if (any(2 * stats::qnorm(.975) * se[i] > log(10), na.rm = TRUE)) {
      out <- rbind(out, data.frame(issue = "stock_recruit_uncertainty", detail =
        "A stock-recruit parameter's 95% interval spans more than tenfold. Curve shape is weakly estimated; compare simpler recruitment models before interpreting density dependence."))
    }
  }
  i <- which(names(p) == "logit_phi_r")
  if (length(i) && is.finite(se[i]) && diff(stats::plogis(p[i] + c(-1, 1) * stats::qnorm(.975) * se[i])) > .5) {
    out <- rbind(out, data.frame(issue = "recruitment_AR1_uncertainty", detail =
      "The recruitment AR1 correlation interval spans more than 0.5. Persistence is imprecise; compare IID recruitment or an externally supported phi."))
  }
  out
}
