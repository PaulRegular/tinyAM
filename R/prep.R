
#' Cut integer sequences into labeled blocks (ages or years)
#'
#' @description
#' `cut_int()` turns an integer vector into labeled
#' blocks using an increasing vector of break points.
#'
#' Labels are of the form `"start-end"` for multi-year/age blocks and
#' `"k"` for single blocks. The last label always ends at `max(x)`
#' (e.g. `"2003-2025"` or `"14"`).
#'
#' Convenience wrappers [cut_ages()] and [cut_years()] call `cut_int()`
#' with argument names that read naturally for common assessment inputs.
#'
#' @details
#' Let `breaks = (b1, b2, ..., bm)`. Blocks are
#' `[b1, b2-1], [b2, b3-1], ..., [b_{m-1}, b_m]` on the integer line.
#' If `b_{m-1} = b_m`, the last block is the singleton `{b_m}`.
#'
#' **Input requirements (enforced):**
#'
#' - `x`, `ages`, `years` are numeric, integer-valued, and non-`NA`.
#' - `breaks` is numeric, increasing, and non-`NA`.
#' - `min(x) == breaks[1]` and `max(x) == tail(breaks, 1)`.
#'
#' @param x,ages,years Integer vector to be grouped.
#' @param breaks Integer vector of group starts (strictly increasing).
#'   Must start at `min(x)` and end at `max(x)`.
#' @param ordered Logical; should the returned factor be ordered?
#'   Default is `FALSE`.
#'
#' @return
#' A factor with the same length as `x`, whose levels enumerate
#' the blocks in increasing order.
#'
#' @examples
#' # Ages: single-year blocks
#' cut_ages(2:14, 2:14)
#'
#' # Ages: a wide block then singletons
#' cut_ages(2:14, c(2, 10, 11:14))
#'
#' # Ages: even-width blocks, last block closes at max age
#' cut_ages(2:14, seq(2, 14, by = 2))
#'
#' # Years: management-era blocks
#' cut_years(1983:2025, c(1983, 1992, 1997, 2003, 2025))
#'
#' # Direct use with ordered = TRUE (if useful for contrasts)
#' cut_int(2:10, c(2, 5, 8, 10), ordered = TRUE)
#'
#' @seealso
#' [base::findInterval()], [stats::model.matrix()]
#'
#' @export
#' @rdname cut_int
cut_int <- function(x, breaks, ordered = FALSE) {

  nm <- deparse1(substitute(x))

  stopifnot(is.numeric(x), is.numeric(breaks))
  if (anyNA(x))                       stop(sprintf("`%s` must be non-NA.", nm), call. = FALSE)
  if (any(x %% 1 != 0))               stop(sprintf("`%s` must be integer-valued.", nm), call. = FALSE)
  if (anyNA(breaks) || any(diff(breaks) <= 0)) stop("`breaks` must be strictly increasing and non-NA.", call. = FALSE)
  if (min(x) != breaks[1])            stop(sprintf("The first break must equal min(%s).", nm), call. = FALSE)
  if (max(x) != tail(breaks, 1))      stop(sprintf("The last break must equal max(%s).", nm), call. = FALSE)

  k <- length(breaks)
  open_end <- k >= 2L && (breaks[k] - breaks[k - 1L] > 1L)
  starts <- if (open_end) breaks[-k] else breaks
  ends <- c(starts[-1L] - 1L, max(breaks))
  labs <- ifelse(starts == ends, starts, paste0(starts, "-", ends))
  idx    <- findInterval(x, starts)

  factor(labs[idx], levels = labs, ordered = ordered)
}

##' @export
##' @rdname cut_int
cut_ages  <- function(ages,  breaks) cut_int(ages,  breaks, ordered = FALSE)

##' @export
##' @rdname cut_int
cut_years <- function(years, breaks) cut_int(years, breaks, ordered = FALSE)



#' Build a self-contained data list for TAM
#'
#' @description
#' `make_dat()` converts tidy observation inputs and modeling options into the
#' structured list `dat` expected by TAM’s likelihood and simulation functions.
#' It expands an age–year grid, merges observations, constructs design matrices for
#' observation SDs, catchability, and mean-\eqn{F} / mean-\eqn{M} (when used),
#' and derives helper mappings and settings.
#'
#' @details
#' **Observation handling**
#'
#' - Inputs are expected as a list with components `catch`, `index`, `weight`,
#'   and `maturity`. Each must include columns `year`, `age`, and a value column
#'   named `obs` (for `catch`/`index`) or renamed from `weight`/`mat`.
#' - Observations are merged to the full `expand.grid(year, age)`.
#' - A combined observation table is created for catch and index; `log(0)` is
#'   treated as `NA` (to be handled via random effects).
#'
#' **Design matrices**
#'
#' - `obs_settings$sd_form` is evaluated on the combined obs map to produce
#'   `sd_obs_modmat`.
#' - `obs_settings$q_form` is evaluated on the index table to produce
#'   `q_modmat`.
#' - If `M_settings$mu_form` is provided, `M_modmat <- model.matrix(mu_form,
#'   data = obs$weight)`. If an `assumption` is also supplied (i.e., fixed
#'   offsets), the intercept in `mu_form` is dropped and a warning is issued.
#' - If neither `M_settings$mu_form` nor `M_settings$assumption` is supplied,
#'   the function stops, because \eqn{M} must be identified by either a fixed
#'   assumption or a mean structure.
#'
#' **Process options and guards**
#'
#' - If `N_settings$process == "off"` and `init_N0 == FALSE`, `init_N0` is
#'   forced to `TRUE` (with a warning) so the first-year abundance is
#'   estimable.
#' - `M_settings$age_breaks` (vector of break points on ages \eqn{\ge} min modeled age + 1)
#'   defines `M_settings$age_blocks` via [cut_ages()], used
#'   to couple \eqn{M} deviations across age.
#' - The AR(1) correlation parameters are only initialized for
#'   processes whose `process == "ar1"`. Correlations are assumed to be 0
#'   when `process == "iid"`, and 0.99 when `process == "approx_rw"` to approximate
#'   a random walk across ages and years.
#'
#' @param obs A list of tidy observation data.frames: `catch`, `index`,
#'   `weight`, and `maturity`. See **Details**.
#' @param years Integer vector of model years (strictly increasing).
#'   Inferred from `obs` if `NULL`.
#' @param ages Integer vector of model ages (strictly increasing).
#'   Inferred from `obs` if `NULL`.
#' @param N_settings A list with elements:
#' - `process`: one of `"off"`, `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `init_N0`: logical; if `TRUE`, estimate an initial level for the
#'   first-year abundance. If `process == "off"` and `init_N0 == FALSE`,
#'   this is forced to `TRUE`.
#' @param F_settings A list with elements:
#' - `process`: one of `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `mu_form`: an optional formula for mean-\eqn{F}.
#' @param M_settings A list with elements:
#' - `process`: one of `"off"`, `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `mu_form`: optional formula for mean-\eqn{M} (on the log scale) built
#'   on `obs$weight`. If provided together with `assumption`, the intercept
#'   in `mu_form` is dropped (warning) so assumed levels act as fixed offsets.
#' - `assumption`: optional one-sided formula giving fixed (non-estimated)
#'   log-\eqn{M} offsets, e.g. `~ I(0.2)` or a column reference such as
#'   `~ log(m_assumption)`.
#' - `age_breaks`: optional integer break points used by [cut_ages()] to
#'   define `age_blocks` for coupling \eqn{M} deviations across ages.
#' @param obs_settings A list with elements:
#' - `sd_form`: formula for observation SD blocks, evaluated on the combined obs map (e.g. `~ sd_obs_block`).
#' - `q_form`: formula for catchability blocks, evaluated on the index table (e.g. `~ q_block`).
#'
#' @return
#' A named list `dat` containing:
#'
#' - `years`, `ages`, `obs` (per-type tables merged to full grid)
#' - `SW`, `MO` weight-at-age and maturity matrices (`year × age`)
#' - `obs_map`, `log_obs`
#' - design matrices: `sd_obs_modmat`, `q_modmat`, and optionally
#'   `F_modmat`, `M_modmat`
#' - mean-level placeholders: `log_mu_f` and/or `log_mu_m` (or `log_mu_assumed_m`)
#' - process settings: `N_settings`, `F_settings`, `M_settings`
#' - AR(1) parameter placeholders `logit_phi_*` set only for `"ar1"` processes
#'
#' @examples
#' dat <- make_dat(
#'   northern_cod_data,
#'   years = 1983:2024,
#'   ages  = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' str(dat)
#'
#' @seealso [stats::model.matrix()], [cut_ages()]
#' @export
make_dat <- function(
    obs,
    years = NULL,
    ages = NULL,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~I(0.2), age_breaks = NULL),
    obs_settings = list(sd_form = ~sd_obs_block, q_form = ~q_block)
) {

  all_years <- unique(unlist(lapply(obs, function(d) d$year)))
  all_ages  <- unique(unlist(lapply(obs, function(d) d$age)))
  years <- if (is.null(years)) seq(min(all_years, na.rm = TRUE), max(all_years, na.rm = TRUE)) else as.integer(years)
  ages  <- if (is.null(ages)) seq(min(all_ages,  na.rm = TRUE), max(all_ages,  na.rm = TRUE)) else as.integer(ages)

  dat <- mget(ls())

  if (N_settings$process == "off" && !N_settings$init_N0) {
    dat$N_settings$init_N0 <- TRUE
    warning("The first year would lack parameters with process set to 'off' and init_N0 set to FALSE in N_settings; forcing init_N0 to TRUE to estimate initial levels.")
  }
  M_ages <- ages[-1]
  if (!is.null(M_settings$age_breaks)) {
    breaks <- dat$M_settings$age_breaks[-1]
    if (!min(M_ages) %in% breaks) breaks <- c(min(M_ages), breaks)
    dat$M_settings$age_blocks <- cut_ages(M_ages, breaks)
  } else {
    dat$M_settings$age_blocks <- cut_ages(M_ages, M_ages)
  }

  grid_ay <- expand.grid(year = years, age = ages)
  dat$obs <- lapply(dat$obs, function(d) {
    d_sub <- merge(grid_ay, d, by = c("age", "year"), all.x = TRUE)
    d_sub[order(d_sub$age, d_sub$year), ] |>
      droplevels()
  })

  dat$SW <- dat$MO <- matrix(NA, nrow = length(years), ncol = length(ages),
                             dimnames = list(year = years, age = ages))
  dat$SW[] <- dat$obs$weight$obs
  dat$MO[] <- dat$obs$maturity$obs

  catch <- dat$obs$catch
  index <- dat$obs$index
  catch[setdiff(names(index), names(catch))] <- NA
  index[setdiff(names(catch), names(index))] <- NA
  catch$type <- "catch"
  index$type <- "index"
  obs_fit <- rbind(catch, index)
  obs_fit$log_obs <- log(obs_fit$obs)
  obs_fit$log_obs[is.infinite(obs_fit$log_obs)] <- NA # treat zeros as NA for simplicity; NAs filled using random effects
  obs_fit$is_na_obs <- is.na(obs_fit$log_obs)

  dat$obs_map <- obs_fit[, setdiff(names(obs_fit), c("obs", "log_obs"))]
  dat$log_obs <- obs_fit$log_obs
  dat$is_na_obs <- obs_fit$is_na_obs

  dat$sd_obs_modmat <- model.matrix(obs_settings$sd_form, data = dat$obs_map)
  dat$q_modmat <- model.matrix(obs_settings$q_form, data = dat$obs$index)
  if (!is.null(dat$F_settings$mu_form)) {
    dat$F_modmat <- model.matrix(F_settings$mu_form, data = dat$obs$catch)
  } else {
    dat$log_mu_f <- 0
    dat$F_modmat <- 0
  }

  if (!is.null(dat$M_settings$mu_form)) {
    dat$M_modmat <- model.matrix(M_settings$mu_form, data = dat$obs$weight)
    if ("(Intercept)" %in% colnames(dat$M_modmat) %in% !is.null(dat$M_settings$assumption)) {
      dat$M_modmat <- model.matrix(update(M_settings$mu_form, ~ 0 + .), data = dat$obs$weight)
      warning("Dropping intercept term in M mu_form since assumed levels are supplied. Set assumption to NULL to estimate the intercept.")
    }
  } else {
    dat$log_mu_m <- 0
    dat$M_modmat <- 0
  }
  if (!is.null(dat$M_settings$assumption)) {
    dat$log_mu_assumed_m <- log(unlist(model.frame(dat$M_settings$assumption, data = dat$obs$weight)))
  } else {
    dat$log_mu_assumed_m <- 0
  }
  if (is.null(dat$M_settings$mu_form) && is.null(dat$M_settings$assumption)) {
    stop("Please supply an assumption or mu_form for M.")
  }

  if (dat$N_settings$process == "iid") {
    dat$logit_phi_n <- c(0, 0)
  }
  if (dat$N_settings$process == "approx_rw") {
    dat$logit_phi_n <- c(0.99, 0.99)
  }

  .set_phi <- function(type = c("off", "iid", "approx_rw", "ar1")) {
    match.arg(type)
    if (type == "iid") {
      return(qlogis(c(0, 0)))
    }
    if (type == "approx_rw") {
      return(qlogis(c(0.99, 0.99)))
    }
    if (type == "ar1") {
      return(NULL)
    }
  }
  dat$logit_phi_n <- .set_phi(dat$N_settings$process)
  dat$logit_phi_f <- .set_phi(dat$F_settings$process)
  dat$logit_phi_m <- .set_phi(dat$M_settings$process)

  dat

}


#' Initialize parameter list for TAM
#'
#' @description
#' `make_par()` builds a named list of initial values and shapes for all
#' fixed and random-effect parameters used by TAM, based on the structure in
#' a previously constructed `dat` list (see [make_dat()]).
#'
#' @details
#' The function inspects `dat` to decide which parameters are required and what
#' their dimensions should be. For example, if `dat$F_settings$process == "ar1"`
#' it initializes a 2-vector `logit_phi_f`; if `dat$F_settings$mu_form` is not
#' `NULL` it creates a coefficient vector `log_mu_f` of length
#' `ncol(dat$F_modmat)`, and so on.
#'
#' All numeric parameters are initialized at `0`, and all matrices are created
#' with appropriate `dimnames` (`year × age` or `year × age_block`).
#'
#' **Created elements (when applicable) include:**
#'
#' - **Recruitment & variability**
#'   - `log_r0` (only if `dat$N_settings$init_N0`)
#'   - `log_r` (length `length(dat$years)`)
#'   - `log_sd_r`
#'
#' - **Abundance deviations (N)**
#'   - `log_sd_n` (if `dat$N_settings$process != "off"`)
#'   - `logit_phi_n` length 2 (if `process == "ar1"`)
#'   - `log_n` matrix (`year` × `age[-1]`) if `process != "off"`
#'
#' - **Fishing mortality (F)**
#'   - `log_sd_f`
#'   - `logit_phi_f` length 2 (if `process == "ar1"`)
#'   - `log_mu_f` coefficients (length `ncol(dat$F_modmat)`) if a mean structure was supplied
#'   - `log_f` matrix (`year` × `age`)
#'
#' - **Natural mortality (M)**
#'   - `log_sd_m` (if `dat$M_settings$process != "off"`)
#'   - `logit_phi_m` length 2 (if `process == "ar1"`)
#'   - `log_mu_m` coefficients (length `ncol(dat$M_modmat)`) if a mean structure was supplied
#'   - `log_m` matrix (`year[-1]` × `age_block`) if `process != "off"`, with
#'     `age_block = levels(dat$M_settings$age_blocks)`
#'
#' - **Observation model**
#'   - `log_q` (length `ncol(dat$q_modmat)`)
#'   - `log_sd_obs` (length `ncol(dat$sd_obs_modmat)`)
#'   - `missing` vector of length `sum(is.na(dat$log_obs))` (placeholders for
#'     imputed `log_obs`)
#'
#' All scalar SD parameters are on the log scale, and AR(1) parameters are on
#' the logit scale (later transformed by `plogis()` in the likelihood).
#'
#' @param dat A data list returned by [make_dat()], containing design matrices,
#'   settings, and observation mappings. The shapes and presence/absence of
#'   parameters depend on elements in `dat` (e.g., `F_modmat`, `M_modmat`,
#'   `q_modmat`, `sd_obs_modmat`, `N_settings`, `F_settings`, `M_settings`,
#'   `years`, `ages`, and `M_settings$age_blocks`).
#'
#' @return
#' A named list of initialized parameters suitable to pass to the TAM objective
#' function, with elements as described in **Details**. All numeric entries are
#' initialized to `0` (or the correct length filled with `0`), and matrices have
#' informative `dimnames`.
#'
#' @examples
#' dat <- make_dat(
#'   northern_cod_data,
#'   years = 1983:2024,
#'   ages  = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' par <- make_par(dat)
#' str(par)
#'
#' @seealso [make_dat()]
#' @export
make_par <- function(dat) {

  par <- list()
  if (dat$N_settings$init_N0) {
    par$log_r0 <- 0
  }
  par$log_sd_r <- 0
  par$log_sd_f <- 0
  if (!is.null(dat$F_settings$mu_form)) {
    par$log_mu_f <- numeric(ncol(dat$F_modmat))
  }
  if (dat$N_settings$process != "off") {
    par$log_sd_n <- 0
  }
  if (dat$M_settings$process != "off") {
    par$log_sd_m <- 0
  }
  if (!is.null(dat$M_settings$mu_form)) {
    par$log_mu_m <- numeric(ncol(dat$M_modmat))
  }
  if (dat$N_settings$process == "ar1") {
    par$logit_phi_n <- c(0, 0)
  }
  if (dat$F_settings$process == "ar1") {
    par$logit_phi_f <- c(0, 0)
  }
  if (dat$M_settings$process == "ar1") {
    par$logit_phi_m <- c(0, 0)
  }
  par$log_q <- numeric(ncol(dat$q_modmat))
  par$log_sd_obs <- numeric(ncol(dat$sd_obs_modmat))

  par$missing <- numeric(sum(is.na(dat$log_obs)))

  par$log_r <- numeric(length(dat$years))
  if (dat$N_settings$process != "off") {
    par$log_n <- matrix(0, nrow = length(dat$years), ncol = length(dat$ages) - 1,
                        dimnames = list(year = dat$years, age = dat$ages[-1]))
  }
  if (dat$M_settings$process != "off") {
    par$log_m <- matrix(0, nrow = length(dat$years) - 1, ncol = nlevels(dat$M_settings$age_blocks),
                        dimnames = list(year = dat$years[-1], age_block = levels(dat$M_settings$age_blocks)))
  }
  par$log_f <- matrix(0, nrow = length(dat$years), ncol = length(dat$ages),
                      dimnames = list(year = dat$years, age = dat$ages))

  par

}


