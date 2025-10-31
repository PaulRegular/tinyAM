
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
#' @importFrom utils tail
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
  if (max(x) != utils::tail(breaks, 1))      stop(sprintf("The last break must equal max(%s).", nm), call. = FALSE)

  k <- length(breaks)
  open_end <- k >= 2L && (breaks[k] - breaks[k - 1L] > 1L)
  starts <- if (open_end) breaks[-k] else breaks
  ends <- c(starts[-1L] - 1L, max(breaks))
  labs <- ifelse(starts == ends, starts, paste0(starts, "-", ends))
  idx    <- findInterval(x, starts)

  out <- factor(labs[idx], levels = labs, ordered = ordered)
  names(out) <- x
  out

}

##' @export
##' @rdname cut_int
cut_ages  <- function(ages,  breaks) cut_int(ages,  breaks, ordered = FALSE)

##' @export
##' @rdname cut_int
cut_years <- function(years, breaks) cut_int(years, breaks, ordered = FALSE)


#' Split cut labels into start/end columns
#'
#' @description
#' Helper for parsing labels produced by [cut_int()], which are either
#' `"start-end"` for multi-width blocks or `"k"` for singletons.
#'
#' @param x A character or factor vector of labels (e.g., from `cut_int()`).
#'
#' @return
#' A data frame with two **character** columns:
#' - `start`: the block start
#' - `end`: the block end
#'
#' @details
#' This function does no validation beyond splitting on `"-"`. Inputs not
#' produced by [cut_int()] may yield unexpected results.
#'
#' @examples
#' cuts <- cut_int(2:10, c(2, 5, 8, 10))
#' .split_cuts(cuts)
#'
#' singletons <- cut_int(2:10, 2:10)
#' .split_cuts(singletons)
#'
#' mix <- cut_int(2:10, c(2:5, 10))
#' .split_cuts(mix)
#'
#' @seealso [cut_int()]
#' @keywords internal
#' @noRd
.split_cuts <- function(x) {
  m <- do.call(rbind, strsplit(as.character(x), "-"))
  if (ncol(m) == 1) m <- cbind(m, m)
  colnames(m) <- c("start", "end")
  as.data.frame(m)
}


#' Append projection rows across core obs tables (internal)
#'
#' @description
#' Internal utility to append `n_proj` years to each of `catch`, `index`,
#' `weight`, and `maturity`.
#'
#' - For `weight` and `maturity`: `obs` in projection years are the
#'   mean over the last `n_mean` terminal years by `age`.
#' - For `catch` and `index`: `obs` in projection years are set to `NA`.
#' - All other columns are copied from the terminal year (by `age`).
#' - Adds a logical `is_proj` column (`FALSE` for historical rows; `TRUE` for projections).
#'
#' @param obs Named list with data.frames `catch`, `index`, `weight`, `maturity`.
#' @param n_proj Integer number of projection years to append (>=1).
#' @param n_mean Integer number of terminal years to average for `weight`/`maturity` (>=1).
#' @return Same `obs` structure with appended projection rows and `is_proj`.
#' @keywords internal
#'
#' @importFrom stats aggregate
#'
#' @noRd
.add_proj_rows <- function(obs, n_proj = 3, n_mean = 3) {
  .add_one <- function(x) {
    max_year   <- max(x$year, na.rm = TRUE)
    proj_years <- seq.int(max_year + 1L, max_year + n_proj)
    mean_years <- seq.int(max_year - n_mean + 1L, max_year)
    aux <- x[x$year == max_year, setdiff(names(x), c("year", "obs")), drop = FALSE]
    mean_obs <- stats::aggregate(
      obs ~ age,
      data = x[x$year %in% mean_years, ],
      FUN = mean
    ) |> merge(aux, by = "age")
    proj_grid <- expand.grid(year = proj_years, age = mean_obs$age)
    proj_rows <- merge(proj_grid, mean_obs, by = "age")[, names(x)]
    proj_rows$is_proj <- TRUE
    x$is_proj <- FALSE
    x_with_proj <- rbind(x, proj_rows)
    x_with_proj[order(x_with_proj$age, x_with_proj$year), ]
  }
  obs_with_proj <- lapply(obs, .add_one)
  obs_with_proj$catch$obs[obs_with_proj$catch$is_proj] <- NA
  obs_with_proj$index$obs[obs_with_proj$index$is_proj] <- NA
  obs_with_proj
}


#' Aggregate data for plus group (internal)
#'
#' @description
#' This internal helper aggregates observations for ages greater than or equal to
#' a specified terminal (`plus_age`) into a single plus group. The `"obs"` column is
#' summed for `catch` and `index` data, and averaged for `weight` and `maturity` data.
#' Non-aggregated columns (e.g., year) retain their values from the terminal age.
#'
#' @param obs Named list with data.frames `catch`, `index`, `weight`, `maturity`.
#' @param plus_age Integer specifying the terminal age to be modeled.
#'   All ages greater than or equal to this value are combined into a plus group.
#' @return A list with the same structure as `obs`, but with all data
#'   aggregated into the specified plus group.
#' @keywords internal
#'
#' @noRd
.plus_fun <- function(obs, plus_age) {

  .aggregate_one <- function(df, plus_age, fun) {
    plus_df <- df[!is.na(df$obs) & df$age >= plus_age, ]
    if (nrow(plus_df) == 0) {
      return(df)
    } else {
      plus_group <- stats::aggregate(obs ~ year, data = plus_df, FUN = fun)
      plus_group$age <- plus_age
      sub_df <- df[df$age <= plus_age, ]
      df_out <- merge(sub_df, plus_group, by = c("year", "age"), all.x = TRUE, suffixes = c("", "_plus"))
      df_out$obs <- ifelse(is.na(df_out$obs_plus), df_out$obs, df_out$obs_plus)
      df_out$obs_plus <- NULL
      return(df_out[order(df_out$age, df_out$year), ])
    }
  }

  list(
    catch = .aggregate_one(obs$catch, plus_age, sum),
    index = .aggregate_one(obs$index, plus_age, sum),
    weight = .aggregate_one(obs$weight, plus_age, mean),
    maturity = .aggregate_one(obs$maturity, plus_age, mean)
  )

}


#' Build a self-contained data list for TAM
#'
#' @description
#' `make_dat()` converts tidy observation inputs and modeling options into the
#' structured list `dat` expected by TAM’s likelihood and simulation functions.
#' It expands an age–year grid, merges observations, constructs design matrices for
#' observation SDs, catchability, and mean-\eqn{F} and/or mean-\eqn{M} (when used),
#' and derives helper mappings and settings.
#'
#' @details
#' **Observation handling**
#'
#' - Inputs are expected as a list with components `catch`, `index`, `weight`,
#'   and `maturity`. Each must include columns `year`, `age`, and a value column
#'   named `obs` (for `catch`/`index`) or renamed from `weight`/`mat`. See
#'   [cod_obs] for an example of the required structure.
#' - Observations are merged to the full `expand.grid(year, age)`.
#' - A combined observation table is created for catch and index; `log(0)` is
#'   treated as `NA` (to be handled via random effects).
#'
#' **Design matrices**
#'
#' - `obs_settings$sd_catch_form` is evaluated on the catch table to produce
#'   `sd_catch_modmat`.
#' - `obs_settings$sd_index_form` is evaluated on the index table to produce
#'   `sd_index_modmat`.
#' - `obs_settings$q_form` is evaluated on the index table to produce
#'   `q_modmat`.
#' - If `M_settings$mu_form` is provided, `M_modmat <- model.matrix(mu_form,
#'   data = obs$weight)`. If an `assumption` is also supplied, the intercept
#'   in `mu_form` is dropped and a warning is issued.
#' - If neither `M_settings$mu_form` nor `M_settings$assumption` is supplied,
#'   the function stops, because \eqn{M} must be identified by either a fixed
#'   assumption or a mean structure.
#'
#' **Process options and guards**
#'
#' - If `N_settings$process == "off"` and `init_N0 == FALSE`, `init_N0` is
#'   forced to `TRUE` (with a warning) so the first-year abundance is
#'   estimable.
#' - `M_settings$age_breaks` (vector of break points on ages)
#'   defines `M_settings$age_blocks` via [cut_ages()], used
#'   to couple \eqn{M} deviations across age.
#' - The AR(1) correlation parameters are only initialized for
#'   processes whose `process == "ar1"`. Correlations are assumed to be 0
#'   when `process == "iid"`, and 0.99 when `process == "approx_rw"` to approximate
#'   a random walk across ages and years.
#'
#' **Projections (optional)**
#'
#' - If `proj_settings` is supplied, the function adds:
#'   - `proj_years`: the set of projection years;
#'   - `is_proj`: logical vector identifying the projection years;
#'   - `obs$...$is_proj` rows to each `obs` table identifying the projection years.
#'
#' @param obs A list of tidy observation data.frames: `catch`, `index`,
#'   `weight`, and `maturity`. See **Details**.
#' @param years Integer vector of model years (strictly increasing).
#'   Inferred from observed data (non-projection) if `NULL`.
#' @param ages Integer vector of model ages (strictly increasing).
#'   Inferred from observed data (non-projection) if `NULL`. If the ages in
#'   the data extend beyond `max(ages)`, the `"obs"` column is summed for `catch`
#'   and `index` data, and averaged for `weight` and `maturity` data.
#' @param N_settings A list with elements:
#' - `process`: one of `"off"`, `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `init_N0`: logical; if `TRUE`, estimate an initial level for the
#'   first-year abundance. If `process == "off"` and `init_N0 == FALSE`,
#'   this is forced to `TRUE`.
#' @param F_settings A list with elements:
#' - `process`: one of `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `mu_form`: an optional formula for mean-\eqn{F}.
#' - `mean_ages`: optional vector of ages to include in population weighted
#'   average F (`F_bar`) calculations. All ages used if absent.
#' @param M_settings A list with elements:
#' - `process`: one of `"off"`, `"iid"`, `"approx_rw"`, or `"ar1"`.
#' - `mu_form`: optional formula for mean-\eqn{M} (on the log scale) built
#'   on `obs$weight`. If provided together with `assumption`, the intercept
#'   in `mu_form` is dropped (warning) so assumed levels act as fixed offsets.
#' - `assumption`: optional one-sided formula giving fixed (non-estimated)
#'   log-\eqn{M} offsets, e.g. `~ I(0.2)` or a column reference such as
#'   `~ log(M_assumption)` stored in the `obs$weight` data.frame.
#' - `age_breaks`: optional integer break points used by [cut_ages()] to
#'   define `age_blocks` for coupling \eqn{M} deviations across ages.
#'   When a narrower set of `age_breaks` than modeled `ages` is provided,
#'   deviations are only estimated for ages within the specified range; ages
#'   outside this range are fixed to their mean or assumed levels. By default,
#'   deviations are estimated for all modeled ages except the youngest, reducing
#'   confounding between \eqn{M} and recruitment.
#' - `first_dev_year`: integer year at which to start estimating \eqn{M} deviations.
#'   Defaults to the second modeled year if `NULL`. Anchoring the first year
#'   to mean or assumed levels helps minimize confounding between early \eqn{M}
#'   estimates, initial abundance, and catchability, leading to more stable estimation.
#' - `mean_ages`: optional vector of ages to include in population weighted
#'   average M (`M_bar`) calculations. All ages used if absent.
#' @param obs_settings A list with elements:
#' - `sd_catch_form`: formula for observation SD blocks for catch-at-age data.
#' - `sd_index_form`: formula for observation SD blocks for index-at-age data.
#' - `q_form`: formula for catchability blocks, evaluated on the index table (e.g. `~ q_block`).
#' - `fill_missing`: logical - fill missing values, and zeros, using random effects? Default = `TRUE`.
#'                   Note that one-step-ahead residuals are not currently working when `TRUE`.
#' @param proj_settings Optional list with elements:
#' - `n_proj`: number of years to project (default `NULL` disables projections).
#' - `n_mean`: number of terminal years to average when adding projection rows. Only the `"obs"` columns
#'             are averaged. All other columns are copied from the terminal year (by `age`).
#' - `F_mult`: multiplier to apply to terminal F to set a level to carry forward in the projection years
#'   (default = `1` to assume status quo F through the projection years). Can be a value of length 1 or
#'   length = `n_proj`. When it is a vector of length one, that multiplier is recycled across all
#'   projection years.
#'
#' @return
#' A named list `dat` containing:
#'
#' - `years`, `ages` — modeled ranges (years includes `proj_years`, if used)
#' - `is_proj` — logical vector identifying whether year is projected
#' - `proj_years` — integer vector of projection years, if used
#' - `obs` — per-type tables restricted to `years` x `ages` (including `proj_years`, if used)
#' - `W`, `P` — mean weight-at-age, and proportion mature at age matrices (`year x age`)
#' - `obs_map` — stack of `obs$catch` and `obs$index` mapping variables
#' - `log_obs`, `is_missing`, `is_observed`, `observed` - vector of log observations (with NA),
#'   logical vector indicating missing and observed values, and vector of non-missing values, respectively.
#' - design matrices: `sd_catch_modmat`, `sd_index_modmat`, `q_modmat`, and optionally `F_modmat`, `M_modmat`
#' - mean-level placeholders: `log_mu_f` and/or `log_mu_m` (or `log_mu_assumed_m`)
#' - process settings: `N_settings`, `F_settings`, `M_settings`
#' - projection settings: `proj_settings`
#' - AR(1) parameter assumptions, `logit_phi_*`, if applicable
#'
#' @example inst/examples/example_dat_default.R
#' @examples
#' names(dat)
#'
#' ## With projection settings
#' dat <- make_dat(
#'   cod_obs,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", assumption = ~ I(0.3)),
#'   obs_settings = list(sd_catch_form = ~ 1, sd_index_form = ~ 1, q_form = ~ q_block),
#'   proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)
#' )
#'
#' @importFrom stats model.frame model.matrix
#'
#' @seealso [stats::model.matrix()], [cut_ages()]
#' @export
make_dat <- function(
    obs,
    years = NULL,
    ages = NULL,
    N_settings = list(process = "iid", init_N0 = FALSE),
    F_settings = list(process = "approx_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, assumption = ~I(0.2), age_breaks = NULL, first_dev_year = NULL),
    obs_settings = list(sd_catch_form = ~1, sd_index_form = ~1, q_form = ~q_block, fill_missing = TRUE),
    proj_settings = NULL
) {

  dat <- mget(ls())

  check_obs(obs)

  ## Subset obs
  all_obs_years <- sort(unique(unlist(lapply(obs, `[[`, "year"))))
  all_obs_ages  <- sort(unique(unlist(lapply(obs, `[[`, "age"))))
  dat$years <- if (is.null(years)) seq(min(all_obs_years), max(all_obs_years)) else as.integer(years)
  dat$ages <- if (is.null(ages)) seq(min(all_obs_ages), max(all_obs_ages)) else as.integer(ages)
  if (max(all_obs_ages) > max(dat$ages)) {
    dat$obs <- .plus_fun(dat$obs, max(dat$ages))
  }
  dat$obs <- lapply(dat$obs, function(d) {
    d_sub <- d[d$year %in% dat$years & d$age %in% dat$ages, ]
    d_sub[order(d_sub$age, d_sub$year), ] |>
      droplevels()
  })
  dat$is_proj <- rep(FALSE, length(dat$years))

  ## Add projection dat
  if (!is.null(proj_settings) && proj_settings$n_proj > 0) {
    dat$obs <- .add_proj_rows(dat$obs, n_proj = proj_settings$n_proj, n_mean = proj_settings$n_mean)
    years_plus <- sort(unique(unlist(lapply(dat$obs, `[[`, "year"))))
    dat$proj_years <- setdiff(years_plus, dat$years)
    dat$is_proj <- c(dat$is_proj, rep(TRUE, proj_settings$n_proj))
    dat$years <- years_plus # update years vec to include proj_years
    if (is.null(proj_settings$F_mult) || any(is.na(proj_settings$F_mult))) {
      cli::cli_abort("{.strong Please specify proj_settings$F_mult (non-NA).}")
    }
    if (length(proj_settings$F_mult) == 1L) {
      dat$proj_settings$F_mult <- rep(proj_settings$F_mult, proj_settings$n_proj)
    } else {
      dat$proj_settings$F_mult <- proj_settings$F_mult
    }
    if (length(dat$proj_settings$F_mult) != dat$proj_settings$n_proj) {
      cli::cli_abort("{.strong length(proj_settings$F_mult) must equal proj_settings$n_proj}")
    }
    if (any(dat$proj_settings$F_mult == 0)) {
      cli::cli_warn("Zero F not supported; replacing 0 with 1e-12.")
      dat$proj_settings$F_mult[dat$proj_settings$F_mult == 0] <- 1e-12
    }
    names(dat$proj_settings$F_mult) <- dat$proj_years
  } else {
    dat$proj_settings <- list(n_proj = 0)
    for (nm in names(dat$obs)) {
      dat$obs[[nm]]$is_proj <- FALSE
    }
  }

  if (N_settings$process == "off" && !N_settings$init_N0) {
    dat$N_settings$init_N0 <- TRUE
    cli::cli_warn("The first year would lack parameters with process set to 'off' and init_N0 set to FALSE in N_settings; forcing init_N0 to TRUE to estimate initial levels.")
  }
  if (!is.null(M_settings$age_breaks)) {
    m_age_range <- range(dat$M_settings$age_breaks)
    age_range <- range(dat$ages)
    m_ages <- seq(max(age_range[1], m_age_range[1]), min(age_range[2], m_age_range[2]), by = 1)
    dat$M_settings$age_blocks <- cut_ages(m_ages, dat$M_settings$age_breaks)
  } else {
    dat$M_settings$age_blocks <- cut_ages(dat$ages[-1], range(dat$ages[-1]))
  }
  if (is.null(M_settings$first_dev_year)) {
    dat$M_settings$first_dev_year <- dat$years[2]
  }
  dat$M_settings$years <- dat$years[dat$years >= dat$M_settings$first_dev_year]
  dat$M_settings$age_block_start <- .split_cuts(levels(dat$M_settings$age_blocks))$start

  empty_mat <- matrix(NA, nrow = length(dat$years), ncol = length(dat$ages),
                      dimnames = list(year = dat$years, age = dat$ages))
  dat$W <- dat$P <- empty_mat
  dat$W[] <- dat$obs$weight$obs
  dat$P[] <- dat$obs$maturity$obs

  catch <- dat$obs$catch
  index <- dat$obs$index
  catch[setdiff(names(index), names(catch))] <- NA
  index[setdiff(names(catch), names(index))] <- NA
  catch$type <- "catch"
  index$type <- "index"
  obs_fit <- rbind(catch, index)
  obs_fit$log_obs <- log(obs_fit$obs)
  obs_fit$log_obs[is.infinite(obs_fit$log_obs)] <- NA # treat zeros as NA for simplicity; NAs filled using random effects
  obs_fit$is_missing <- is.na(obs_fit$log_obs)
  obs_fit$is_observed <- !obs_fit$is_missing

  if (is.null(obs_settings$fill_missing)) {
    dat$obs_settings$fill_missing <- TRUE
    cli::cli_warn("obs_settings$fill_missing was NULL; forcing to TRUE")
  }
  if (sum(obs_fit$is_missing) == 0) {
    dat$obs_settings$fill_missing <- FALSE # set fill_missing to FALSE if there are no missing values
  }

  dat$obs_map <- obs_fit[, setdiff(names(obs_fit), c("obs", "log_obs"))]
  dat$log_obs <- obs_fit$log_obs
  dat$is_missing <- obs_fit$is_missing
  dat$is_observed <- obs_fit$is_observed
  dat$observed <- dat$log_obs[!dat$is_missing]

  dat$sd_catch_modmat <- stats::model.matrix(obs_settings$sd_catch_form, data = dat$obs$catch)
  dat$sd_index_modmat <- stats::model.matrix(obs_settings$sd_index_form, data = dat$obs$index)
  dat$q_modmat <- stats::model.matrix(obs_settings$q_form, data = dat$obs$index)
  if (!is.null(dat$F_settings$mu_form)) {
    dat$F_modmat <- stats::model.matrix(F_settings$mu_form, data = dat$obs$catch)
  } else {
    dat$log_mu_f <- 0
    dat$F_modmat <- 0
  }

  if (!is.null(dat$M_settings$mu_form)) {
    dat$M_modmat <- stats::model.matrix(M_settings$mu_form, data = dat$obs$weight)
    if ("(Intercept)" %in% colnames(dat$M_modmat) && !is.null(dat$M_settings$assumption)) {
      dat$M_modmat <- stats::model.matrix(update(M_settings$mu_form, ~ 0 + .), data = dat$obs$weight)
      cli::cli_warn("Dropping intercept term in M mu_form since assumed levels are supplied. Set assumption to NULL to estimate the intercept.")
    }
  } else {
    dat$log_mu_m <- 0
    dat$M_modmat <- 0
  }
  if (!is.null(dat$M_settings$assumption)) {
    dat$log_mu_assumed_m <- log(unlist(stats::model.frame(dat$M_settings$assumption, data = dat$obs$weight)))
  } else {
    dat$log_mu_assumed_m <- 0
  }
  if (is.null(dat$M_settings$mu_form) && is.null(dat$M_settings$assumption)) {
    stop("Please supply an assumption or mu_form for M.")
  }

  .check_ages <- function(x, ages, label) {
    if (is.null(x)) return(ages)
    if (!all(x %in% ages)) {
      cli::cli_abort("{.strong {label}} must be a subset of ages: {paste(ages, collapse = ', ')}")
    }
    x
  }
  dat$F_settings$mean_ages <- .check_ages(dat$F_settings$mean_ages, dat$ages, "F_settings$mean_ages")
  dat$M_settings$mean_ages <- .check_ages(dat$M_settings$mean_ages, dat$ages, "M_settings$mean_ages")

  .set_phi <- function(type = c("off", "iid", "approx_rw", "ar1")) {
    match.arg(type)
    if (type == "iid") {
      return(qlogis(c("age" = 0, "year" = 0)))
    }
    if (type == "approx_rw") {
      return(qlogis(c("age" = 0.99, "year" = 0.99)))
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

