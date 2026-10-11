.parse_process_sd <- function(dat, component) {
  settings <- dat[[paste0(component, "_settings")]]
  formula <- settings$sd_form
  if (is.null(formula)) formula <- ~ 1
  label <- paste0(component, "_settings$sd_form")
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    cli::cli_abort("{label} must be a one-sided formula.")
  }
  default <- identical(formula[[2L]], 1) || identical(formula[[2L]], 1L)
  if (settings$process == "off") {
    if (!default) cli::cli_abort("{label} has no process to describe when process = 'off'.")
    return(NULL)
  }
  walk <- function(x) {
    if (!is.call(x)) return(FALSE)
    .formula_call_name(x) %in% c("iid", "rw", "ar1", "mono", "logistic", "bh", "ricker", "|") ||
      any(vapply(as.list(x)[-1L], walk, logical(1)))
  }
  if (walk(formula[[2L]]) || "year" %in% all.vars(formula)) {
    cli::cli_abort("{label} supports ordinary age-based formulas, not temporal or random-effect SD terms.")
  }
  d <- dat$obs$weight
  active <- if (component == "N") dat$ages[-1L] else if (component == "M") {
    as.numeric(names(settings$age_blocks))
  } else dat$ages
  d <- d[d$age %in% active & !d$is_proj, , drop = FALSE]
  frame <- stats::model.frame(formula, d, na.action = stats::na.pass)
  if (anyNA(frame)) cli::cli_abort("{label} needs complete covariates at every active age and historical year.")
  design <- stats::model.matrix(formula, frame)
  if (!ncol(design) || any(!is.finite(design))) cli::cli_abort("{label} must produce a nonempty finite design.")
  for (a in active) {
    rows <- design[d$age == a, , drop = FALSE]
    if (any(rows != rep(rows[1L, ], each = nrow(rows)))) {
      cli::cli_abort("{label} varies across years at age {a}. Temporal process SDs are not enabled.")
    }
  }
  design <- design[match(active, d$age), , drop = FALSE]
  rownames(design) <- as.character(active)
  if (component == "M") {
    for (block in levels(settings$age_blocks)) {
      rows <- design[settings$age_blocks == block, , drop = FALSE]
      if (any(rows != rep(rows[1L, ], each = nrow(rows)))) {
        cli::cli_abort("{label} varies within an M age block. Use finer age_breaks or share the SD within each block.")
      }
    }
    design <- design[as.character(settings$age_block_start), , drop = FALSE]
  }
  if (qr(design)$rank < ncol(design)) cli::cli_abort("{label} is rank deficient on the active process states. Remove redundant terms.")
  list(matrix = design, default = default,
    parameter = if (default) paste0("log_sd_", tolower(component)) else paste0("sd_beta_", tolower(component)))
}

.process_sd <- function(par, dat, component) {
  design <- dat$process_sd[[component]]
  if (design$default) return(exp(par[[design$parameter]]))
  exp(drop(design$matrix %*% par[[design$parameter]]))
}

.process_sd_surface <- function(par, dat, component) {
  design <- dat$process_sd[[component]]
  matrix(rep(.process_sd(par, dat, component), each = length(dat$years)),
    length(dat$years), nrow(design$matrix),
    dimnames = list(year = dat$years, age = rownames(design$matrix)))
}

# Standardization retains the existing AR1 innovation interpretation.
.dprocess_scaled <- function(x, sd, process, phi = c(0, 0)) {
  if (length(sd) == 1L) sd <- rep(sd, ncol(x))
  if (process == "rw") {
    if (nrow(x) < 2L) return(0)
    x <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  }
  z <- sweep(x, 2, sd, "/")
  density <- if (process == "ar1") dprocess_ar1(z, phi = phi, sd = 1) else sum(RTMB::dnorm(z, 0, 1, log = TRUE))
  density - nrow(x) * sum(log(sd))
}

.rprocess_scaled <- function(x, sd, process, phi = c(0, 0)) {
  if (process == "rw") {
    if (nrow(x) > 1L) for (y in 2:nrow(x)) x[y, ] <- x[y - 1L, ] + stats::rnorm(ncol(x), 0, sd)
    return(x)
  }
  z <- if (process == "ar1") rprocess_ar1(nrow(x), ncol(x), phi = phi, sd = 1) else
    matrix(stats::rnorm(length(x)), nrow(x), ncol(x))
  out <- sweep(z, 2, sd, "*")
  dimnames(out) <- dimnames(x)
  out
}

.process_sd_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  for (component in names(dat$process_sd)) {
    design <- dat$process_sd[[component]]
    if (is.null(design) || design$default) next
    settings <- dat[[paste0(component, "_settings")]]
    rows <- if (component == "M") length(settings$years) else sum(!dat$is_proj) - as.integer(component == "N")
    if (settings$process %in% c("rw", "cor_rw")) rows <- rows - 1L
    support <- colSums(design$matrix != 0) * rows
    if (any(support < 5L)) out <- rbind(out, data.frame(issue = "process_SD_support",
      detail = paste0(component, " SD design has fewer than five historical process contributions supporting some coefficients. Share more ages or examine simulation recovery and uncertainty.")))
  }
  out
}

.process_sd_fit_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  for (component in names(fit$dat$process_sd)) {
    d <- fit$pop[[paste0("sd_", component)]]
    if (is.null(d)) next
    if (any(is.finite(d$lwr) & is.finite(d$upr) & d$lwr > 0 & d$upr / d$lwr > 10)) {
      out <- rbind(out, data.frame(issue = "process_SD_uncertainty",
        detail = paste0(component, " process SD intervals span more than tenfold at some ages. Examine simpler sharing groups before interpreting age differences.")))
    }
  }
  out
}
