
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
#' Breaks normally mark the first value in each block. The final break must
#' equal `max(x)`. If the gap between the last two breaks exceeds one, the final
#' break closes the preceding block rather than starting a singleton. For
#' example, `c(2, 5, 8)` gives `2-4` and `5-8`, whereas `c(2, 5, 7, 8)` gives
#' `2-4`, `5-6`, `7`, and `8`. This convention also defines M age blocks.
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

  if (!is.numeric(x) || !is.numeric(breaks)) {
    cli::cli_abort("{.arg x} and {.arg breaks} must be numeric.")
  }
  if (anyNA(x)) {
    cli::cli_abort("{.arg {nm}} must be non-NA.")
  }
  if (any(x %% 1 != 0)) {
    cli::cli_abort("{.arg {nm}} must be integer-valued.")
  }
  if (anyNA(breaks) || any(diff(breaks) <= 0)) {
    cli::cli_abort("{.arg breaks} must be strictly increasing and non-NA.")
  }
  if (min(x) != breaks[1]) {
    cli::cli_abort(sprintf("The first break must equal min(%s).", nm))
  }
  if (max(x) != utils::tail(breaks, 1)) {
    cli::cli_abort(sprintf("The last break must equal max(%s).", nm))
  }

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
  max_year <- max(unlist(lapply(obs, `[[`, "year")))
  .add_one <- function(x, nm) {
    proj_years <- seq.int(max_year + 1L, max_year + n_proj)
    mean_years <- seq.int(max_year - n_mean + 1L, max_year)
    template <- x[x$year == max(x$year), , drop = FALSE]
    if (nm %in% c("weight", "maturity")) {
      mean_obs <- stats::aggregate(obs ~ age, data = x[x$year %in% mean_years, ], FUN = mean)
      template$obs <- mean_obs$obs[match(template$age, mean_obs$age)]
    } else {
      template$obs <- NA_real_
    }
    proj_rows <- do.call(rbind, lapply(proj_years, function(year) {
      rows <- template
      rows$year <- year
      rows
    }))
    proj_rows$is_proj <- TRUE
    x$is_proj <- FALSE
    x_with_proj <- rbind(x, proj_rows)
    x_with_proj[order(x_with_proj$age, x_with_proj$year), ]
  }
  obs_with_proj <- lapply(names(obs), function(nm) .add_one(obs[[nm]], nm))
  names(obs_with_proj) <- names(obs)
  obs_with_proj$catch$obs[obs_with_proj$catch$is_proj] <- NA
  obs_with_proj$index$obs[obs_with_proj$index$is_proj] <- NA
  obs_with_proj
}


#' Aggregate data for plus group (internal)
#'
#' @description
#' This internal helper aggregates observations for ages greater than or equal to
#' a specified terminal (`plus_age`) into a single plus group. The `"obs"` column is
#' summed for `catch` and `index` data. Weight and maturity tables retain all
#' ages for abundance-based aggregation inside [nll_fun()]. Other catch/index
#' columns retain their values from the plus age or first available older age.
#'
#' @param obs Named list with data.frames `catch`, `index`, `weight`, `maturity`.
#' @param plus_age Integer specifying the terminal age to be modeled.
#'   All ages greater than or equal to this value are combined into a plus group.
#' @return A list with the same structure as `obs`, with catch/index observations
#'   aggregated into the specified plus group and biological tables unchanged.
#' @keywords internal
#'
#' @noRd
.plus_fun <- function(obs, plus_age) {

  .aggregate_one <- function(df, plus_age, fun, groups = "year") {
    plus_df <- df[df$age >= plus_age, , drop = FALSE]
    if (nrow(plus_df) == 0) {
      return(df)
    } else {
      plus_group <- stats::aggregate(plus_df["obs"], plus_df[groups], FUN = function(x) {
        if (all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
      })
      # Retain metadata at the plus age, or the first available older age.
      template <- plus_df[order(plus_df$age), , drop = FALSE]
      template <- template[!duplicated(template[groups]), , drop = FALSE]
      template$age <- plus_age
      plus_rows <- merge(template[, setdiff(names(template), "obs"), drop = FALSE],
                         plus_group, by = groups)
      df_out <- rbind(df[df$age < plus_age, , drop = FALSE], plus_rows[, names(df), drop = FALSE])
      return(df_out[order(df_out$age, df_out$year), ])
    }
  }

  list(
    catch = .aggregate_one(obs$catch, plus_age, sum),
    index = .aggregate_one(obs$index, plus_age, sum, groups = c("year", "survey")),
    weight = obs$weight,
    maturity = obs$maturity
  )

}


#' Build a self-contained data list for TAM
#'
#' @description
#' Prepares catch, survey, weight, and maturity data for the chosen model years
#' and ages. Checks input coverage, combines older ages into a plus group when
#' needed, and sets up mortality, catchability, and observation-error formulas.
#' Optional projection settings append future years. [fit_tam()] calls this
#' function automatically; call it directly to inspect a model specification.
#'
#' @details
#' **Observation handling**
#'
#' - Inputs are expected as a list with components `catch`, `index`, `weight`,
#'   and `maturity`. Each must include columns `year`, `age`, and a value column
#'   named `obs`. Rename any differently named input columns before calling. See
#'   [cod_obs] for an example of the required structure.
#' - Catch, weight, and maturity must already contain exactly one row per
#'   year-age combination across the input range. Survey tables may be sparse.
#' - A combined observation table is created for catch and index; `log(0)` is
#'   treated as `NA` (to be handled via random effects).
#' - When biological inputs extend above the modeled plus age, the underlying
#'   weight and maturity values are retained in `W_plus_input` and `P_plus_input`,
#'   indexed by `years` and `plus_ages`. Their projection rows use the same
#'   recent-year averages by biological age as other weight/maturity inputs.
#'   `dat$W`, `dat$P`, and `dat$obs` retain the input values at modeled ages;
#'   their terminal biological values are placeholders. [nll_fun()] replaces
#'   terminal W by an abundance-weighted mean and terminal P by a biomass-weighted
#'   proportion using its reconstructed hidden age composition. See [tinyAM-model].
#'   Formula covariates continue to use the input rows at modeled ages, not these
#'   effective biological values; this avoids circular mortality calculations.
#'
#' **Design matrices**
#'
#' - `catch_settings$sd_form` is evaluated on the catch table to produce
#'   `sd_catch_modmat` and associated **log-scale** parameters `log_sd_catch`
#'   (intercept removed when `sd_supplied` is provided so supplied SDs act as
#'   offsets).
#' - `index_settings$sd_form` is evaluated on the index table to produce
#'   `sd_index_modmat` and **log-scale** parameters `log_sd_index`.
#' - `index_settings$q_form` is evaluated on the index table to produce
#'   `q_modmat` and **log-scale** parameters `log_q`.
#'   Additive [mono()] terms instead contribute cumulative indicators
#'   in `q_mono_modmat`, with directly fitted non-negative log-q increments `dq`.
#'   `q_mono_steps` records each transition and its optional group.
#' - If `M_settings$mu_form` is provided, `M_modmat <- model.matrix(mu_form,
#'   data = obs$weight)` and the resulting coefficients are parameters `mu_m`.
#'   These coefficients are applied on the log scale to build \eqn{M}, but
#'   retain the `mu_` prefix to emphasize they can be positive or negative.
#'   If `mu_supplied` is also provided, the intercept in `mu_form` is dropped
#'   and a warning is issued.
#' - If neither `M_settings$mu_form` nor `M_settings$mu_supplied` is supplied,
#'   the function stops, because \eqn{M} must be identified by either a supplied
#'   surface or a mean structure.
#'
#' **Process options and guards**
#'
#' - Initial abundance is controlled by `N_settings$init`, independently of
#'   the subsequent N process. First-year recruitment is the fixed state
#'   `log_r0`. All older initial ages use survivorship from recruitment under
#'   first-year total mortality Z as their baseline, with no equilibrium
#'   plus-group correction.
#' - `"exp"` is the default: a parsimonious, comparatively stable nuisance-state
#'   initialization. More flexible choices are available when first-year data
#'   support them. `"free"` estimates unpenalized fixed `log_n0` states;
#'   `"random"` estimates random `log_n0` states whose survivorship residuals
#'   `eta_log_n0` are IID normal with separately estimated `sd_n0`.
#' - `sd_r` describes temporal recruitment variation; `sd_n0` describes
#'   variation across historical cohorts on the initial age margin (potentially
#'   including mortality-history variation); `sd_n` describes subsequent
#'   cohort-process deviations. These SDs are estimated separately.
#' - Random initialization requires at least two ages. Fewer than ten ages
#'   triggers a heuristic weak-identification warning for `sd_n0`, not a
#'   prohibition or automatic fallback. The threshold may be revised after
#'   simulation testing; inspect convergence and sensitivity carefully.
#' - `M_settings$age_breaks` (vector of break points on ages)
#'   defines `M_settings$age_blocks` via [cut_ages()], used
#'   to share absolute latent \eqn{M} states across ages. All mean components
#'   must be constant across ages within each fitted block.
#' - The AR(1) correlation parameters are only initialized for
#'   processes whose `process == "ar1"`. Correlations are assumed to be 0
#'   when `process == "iid"`. A temporal random walk (`"rw"`) has independent
#'   year-to-year increments within each age or age block, with no age
#'   correlation or AR parameter. Its first process row is unpenalized.
#'   For N the walk acts on cohort residuals; for F/M it acts on deviations
#'   from their existing log mean surfaces. Initial abundance, the M process
#'   start year, and terminal-F projections retain their usual meanings.
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
#' @param years Consecutive historical model years (at least two).
#'   Inferred from observed data (non-projection) if `NULL`.
#' @param ages Consecutive model ages (at least two).
#'   Inferred from observed data (non-projection) if `NULL`. If the ages in
#'   the data extend beyond `max(ages)`, the `"obs"` column is summed for `catch`
#'   and `index` data within each year (and survey for indices). Older-age
#'   weight and maturity are retained for abundance-based aggregation inside
#'   [nll_fun()], preserving biomass and mature biomass in the plus group.
#' @param N_settings A list with elements:
#' - `process`: `"off"` for deterministic cohort survival, `"iid"` for independent
#'   cohort residuals, `"rw"` for residuals that accumulate through time, or
#'   `"ar1"` for residuals correlated between years and ages. See [tinyAM-model].
#' - `init`: `"exp"` (default), `"free"`, or `"random"`, independently of
#'   `process`. All use fixed first-year recruitment `log_r0` as the starting
#'   anchor. `"exp"` uses deterministic survivorship; `"free"` estimates fixed
#'   older-age `log_n0` states; `"random"` estimates random `log_n0` states with
#'   IID survivorship residuals and separate SD `sd_n0`. See **Details**.
#' @param F_settings A list with elements:
#' - `process`: `"iid"` for independent departures from the mean log F,
#'   `"rw"` for departures that accumulate over time, or `"ar1"` for departures
#'   correlated between years and ages.
#' - `mu_form`: an optional formula for mean-\eqn{F} (coefficients estimated as
#'   **log-scale** parameters `log_mu_f`).
#' - `mean_ages`: optional vector of ages to include in population weighted
#'   average F (`F_bar`) calculations. All ages used if absent.
#' @param M_settings A list with elements:
#' - `process`: `"off"` uses only supplied/mean M. Otherwise `"iid"`, `"rw"`,
#'   and `"ar1"` have the same meanings as for F, on selected years/age blocks.
#' - `mu_form`: optional formula for mean-\eqn{M} (applied on the log scale) built
#'   on `obs$weight`, yielding coefficients `mu_m`. These enter the log-\eqn{M}
#'   surface directly and may therefore be positive or negative; they intentionally
#'   avoid a `log_` prefix to prevent automatic back-transformation when tidied.
#'   If provided together with `mu_supplied`, the intercept in `mu_form` is
#'   dropped (warning) so supplied levels act as fixed offsets.
#' - `mu_supplied`: optional one-sided formula giving supplied (non-estimated)
#'   \eqn{M}, e.g. `~ I(0.2)` or a column reference such as
#'   `~ M_assumption` stored in the `obs$weight` data.frame.
#' - `age_breaks`: optional integer break points used by [cut_ages()] to
#'   define `age_blocks` sharing absolute latent M states across ages.
#'   When a narrower set of `age_breaks` than modeled `ages` is provided,
#'   process states are only estimated for ages within the specified range; ages
#'   outside this range are fixed to their mean or assumed levels. By default,
#'   blocks use `cut_ages(ages[-1], unique(range(ages[-1])))`, usually one block
#'   excluding the youngest age; two older ages form separate singleton blocks.
#'   This restricts
#'   how M can trade off against recruitment; it does not guarantee identifiability.
#' - `first_dev_year`: one historical modeled year at which M process states
#'   begin. Defaults to the second year if `NULL`. Earlier M stays at its mean
#'   or supplied value; earlier years are not coupled to the first latent state.
#'   Ignored when the M process is off.
#' - `mean_ages`: optional vector of ages to include in population weighted
#'   average M (`M_bar`) calculations. All ages used if absent.
#' @param catch_settings A list with elements:
#' - `sd_form`: formula for observation SD blocks for catch-at-age data.
#' - `sd_supplied`: optional one-sided formula giving supplied SDs (on the natural
#'   scale of the log-observation residuals) for catch-at-age data. When provided,
#'   the intercept is removed from `sd_form` so supplied SDs act as offsets.
#' - `fill_missing`: logical – fill missing values, and zeros, using random effects?
#'   Defaults to `TRUE`. Note that one-step-ahead residuals are not currently working when `TRUE`.
#' @param index_settings A list with elements:
#' - `sd_form`: formula for observation SD blocks for index-at-age data.
#' - `sd_supplied`: optional one-sided formula giving supplied SDs (on the natural
#'   scale of the log-observation residuals) for index-at-age data. When provided,
#'   the intercept is removed from `sd_form` so supplied SDs act as offsets.
#' - `q_form`: formula for catchability, evaluated on the index table. Ordinary
#'   `~ q_block` is unconstrained; `~ mono(q_block)` is non-decreasing across ordered
#'   blocks. See [mono()] for independent non-decreasing curves by survey.
#' - `fill_missing`: logical – fill missing values, and zeros, using random effects?
#'   Defaults to `TRUE`. Note that one-step-ahead residuals are not currently working when `TRUE`.
#' @param proj_settings Optional list with elements:
#' - `n_proj`: number of years to project (default `NULL` disables projections).
#' - `n_mean`: number of recent historical years used to average weight and
#'   maturity by age. Catch/index observations in projection years are missing.
#'   Other columns are copied from the last available year in each table.
#' - `F_mult`: multiplier to apply to terminal F to set a level to carry forward in the projection years
#'   (required when projecting; use `1` for status quo F). Can be a value of length 1 or
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
#' - `plus_ages`, `W_plus_input`, `P_plus_input` — retained biological ages and
#'   inputs, present only when biological data extend above `max(ages)`. The last
#'   retained age represents a hidden plus group. Effective modeled W/P and hidden
#'   abundance `N_plus` can then be inspected in the likelihood report (`fit$rep`).
#' - `obs_map` — stack of `obs$catch` and `obs$index` mapping variables
#' - `log_obs`, `is_missing`, `is_observed`, `observed` - vector of log observations (with NA),
#'   logical vector indicating missing and observed values, and vector of non-missing values, respectively.
#' - design matrices: `sd_catch_modmat`, `sd_index_modmat`, `q_modmat`, and optionally `F_modmat`, `M_modmat`
#' - mean-level placeholders: `log_mu_f` and/or `mu_m` (or `log_mu_supplied_m`)
#' - process settings: `N_settings`, `F_settings`, `M_settings`, `catch_settings`, `index_settings`
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
#'   N_settings = list(process = "iid", init = "exp"),
#'   F_settings = list(process = "rw", mu_form = NULL),
#'   M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
#'   catch_settings = list(sd_form = ~ 1),
#'   index_settings = list(sd_form = ~ 1, q_form = ~ q_block),
#'   proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)
#' )
#'
#' @importFrom stats model.frame model.matrix
#'
#' @seealso [fit_tam()], [tinyAM-model], [stats::model.matrix()], [cut_ages()]
#' @export
make_dat <- function(
    obs,
    years = NULL,
    ages = NULL,
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~I(0.2), age_breaks = NULL, first_dev_year = NULL),
    catch_settings = list(sd_form = ~1, sd_supplied = NULL, fill_missing = TRUE),
    index_settings = list(sd_form = ~1, sd_supplied = NULL, q_form = ~q_block, fill_missing = TRUE),
    proj_settings = NULL
) {

  dat <- mget(ls())

  dat$N_settings$process <- match.arg(dat$N_settings$process, c("off", "iid", "rw", "ar1"))
  dat$F_settings$process <- match.arg(dat$F_settings$process, c("iid", "rw", "ar1"))
  dat$M_settings$process <- match.arg(dat$M_settings$process, c("off", "iid", "rw", "ar1"))

  check_obs(obs)

  for (nm in c("years", "ages")) {
    x <- get(nm)
    if (!is.null(x) && (!is.numeric(x) || !length(x) ||
        any(!is.finite(x)) || any(x != trunc(x)))) {
      cli::cli_abort("{.arg {nm}} must be a nonempty vector of integer values.")
    }
  }

  ## Subset obs
  all_obs_years <- sort(unique(unlist(lapply(obs, `[[`, "year"))))
  all_obs_ages  <- sort(unique(unlist(lapply(obs, `[[`, "age"))))
  dat$years <- if (is.null(years)) seq(min(all_obs_years), max(all_obs_years)) else as.integer(years)
  dat$ages <- if (is.null(ages)) seq(min(all_obs_ages), max(all_obs_ages)) else as.integer(ages)
  if (min(dat$years) < min(all_obs_years) ||
      max(dat$years) > max(all_obs_years)) {
    cli::cli_abort(c(
      "{.arg years} must fall within the years available in {.arg obs}.",
      "x" = "Requested years: {min(dat$years)}-{max(dat$years)}.",
      "i" = "Available years: {min(all_obs_years)}-{max(all_obs_years)}.",
      "i" = "Use {.arg proj_settings} for years beyond the terminal data year."
    ))
  }
  if (min(dat$ages) < min(all_obs_ages) ||
      max(dat$ages) > max(all_obs_ages)) {
    cli::cli_abort(c(
      "{.arg ages} must fall within the ages available in {.arg obs}.",
      "x" = "Requested ages: {min(dat$ages)}-{max(dat$ages)}.",
      "i" = "Available ages: {min(all_obs_ages)}-{max(all_obs_ages)}."
    ))
  }
  if (any(diff(dat$years) != 1L)) {
    cli::cli_abort("{.arg years} must be consecutive years.")
  }
  if (any(diff(dat$ages) != 1L)) {
    cli::cli_abort("{.arg ages} must be consecutive ages.")
  }
  if (max(all_obs_ages) > max(dat$ages)) {
    dat$obs <- .plus_fun(dat$obs, max(dat$ages))
  }
  dat$obs <- stats::setNames(lapply(names(dat$obs), function(nm) {
    d <- dat$obs[[nm]]
    keep_age <- if (nm %in% c("weight", "maturity")) d$age >= min(dat$ages) else d$age %in% dat$ages
    d_sub <- d[d$year %in% dat$years & keep_age, ]
    d_sub[order(d_sub$age, d_sub$year), ] |>
      droplevels()
  }), names(dat$obs))
  dat$is_proj <- rep(FALSE, length(dat$years))
  if (length(dat$years) < 2L) {
    cli::cli_abort("At least two historical years are required for recruitment and cohort transitions.")
  }

  ## Add projection dat
  if (!is.null(proj_settings) && proj_settings$n_proj > 0) {
    dat$obs <- .add_proj_rows(dat$obs, n_proj = proj_settings$n_proj, n_mean = proj_settings$n_mean)
    years_plus <- sort(unique(unlist(lapply(dat$obs, `[[`, "year"))))
    dat$proj_years <- setdiff(years_plus, dat$years)
    dat$is_proj <- years_plus %in% dat$proj_years
    dat$years <- years_plus # update years vec to include proj_years
    if (is.null(proj_settings$F_mult) || any(is.na(proj_settings$F_mult))) {
      cli::cli_abort("{.strong Please specify proj_settings$F_mult (non-NA).}")
    }
    if (!is.numeric(proj_settings$F_mult) || any(!is.finite(proj_settings$F_mult)) ||
        any(proj_settings$F_mult < 0)) {
      cli::cli_abort("proj_settings$F_mult must contain finite, non-negative multipliers.")
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

  if ("init_N0" %in% names(N_settings)) {
    cli::cli_abort("N_settings$init_N0 has been retired. Use N_settings$init = \"exp\", \"free\", or \"random\" instead.")
  }
  if (is.null(N_settings$init)) dat$N_settings$init <- "exp"
  dat$N_settings$init <- match.arg(dat$N_settings$init, c("exp", "free", "random"))
  if (length(dat$ages) < 2L) {
    cli::cli_abort("N0 initialization requires at least two modeled ages.")
  }
  if (dat$N_settings$init == "random") {
    n_init_dev <- length(dat$ages) - 1L
    if (n_init_dev < 1L) {
      cli::cli_abort("Random N0 initialization requires at least one initial age-to-age deviation (two modeled ages).")
    }
    if (length(dat$ages) < 10L) {
      cli::cli_warn('"random" N0 initialization estimates sd_n0 from only {n_init_dev} initial age-to-age deviations ({length(dat$ages)} modeled ages). sd_n0 may be weakly identified. Consider N_settings$init = "exp", or inspect convergence and sensitivity carefully.')
    }
  }
  if (!is.null(M_settings$age_breaks)) {
    m_age_range <- range(dat$M_settings$age_breaks)
    age_range <- range(dat$ages)
    m_ages <- seq(max(age_range[1], m_age_range[1]), min(age_range[2], m_age_range[2]), by = 1)
    dat$M_settings$age_blocks <- cut_ages(m_ages, dat$M_settings$age_breaks)
  } else {
    dat$M_settings$age_blocks <- cut_ages(dat$ages[-1], unique(range(dat$ages[-1])))
  }
  if (is.null(M_settings$first_dev_year)) {
    dat$M_settings$first_dev_year <- dat$years[2]
  }
  if (dat$M_settings$process != "off") {
    start <- dat$M_settings$first_dev_year
    if (length(start) != 1L || !is.numeric(start) || !is.finite(start) ||
        start != as.integer(start) || !start %in% dat$years[!dat$is_proj]) {
      cli::cli_abort("M_settings$first_dev_year must be a single historical modeled year.")
    }
  }
  dat$M_settings$years <- dat$years[dat$years >= dat$M_settings$first_dev_year]
  dat$M_settings$age_block_start <- .split_cuts(levels(dat$M_settings$age_blocks))$start

  empty_mat <- matrix(NA, nrow = length(dat$years), ncol = length(dat$ages),
                      dimnames = list(year = dat$years, age = dat$ages))
  if (max(dat$obs$weight$age) > max(dat$ages)) {
    dat$plus_ages <- seq.int(max(dat$ages), max(dat$obs$weight$age))
    plus_dimnames <- list(year = dat$years, age = dat$plus_ages)
    dat$W_plus_input <- matrix(dat$obs$weight$obs[dat$obs$weight$age %in% dat$plus_ages],
                               length(dat$years), length(dat$plus_ages), dimnames = plus_dimnames)
    dat$P_plus_input <- matrix(dat$obs$maturity$obs[dat$obs$maturity$age %in% dat$plus_ages],
                               length(dat$years), length(dat$plus_ages), dimnames = plus_dimnames)
    for (nm in c("weight", "maturity")) {
      dat$obs[[nm]] <- droplevels(dat$obs[[nm]][dat$obs[[nm]]$age %in% dat$ages, ])
    }
  }
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

  dat$catch_settings <- catch_settings
  dat$index_settings <- index_settings
  if (is.null(dat$catch_settings$fill_missing)) {
    dat$catch_settings$fill_missing <- TRUE
    cli::cli_warn("catch_settings$fill_missing was NULL; forcing to TRUE")
  }
  if (is.null(dat$index_settings$fill_missing)) {
    dat$index_settings$fill_missing <- TRUE
    cli::cli_warn("index_settings$fill_missing was NULL; forcing to TRUE")
  }

  dat$obs_map <- obs_fit[, setdiff(names(obs_fit), c("obs", "log_obs"))]
  dat$log_obs <- obs_fit$log_obs
  dat$is_missing <- obs_fit$is_missing
  dat$is_observed <- obs_fit$is_observed
  dat$observed <- dat$log_obs[!dat$is_missing]

  catch_missing <- dat$is_missing & dat$obs_map$type == "catch"
  index_missing <- dat$is_missing & dat$obs_map$type == "index"
  if (sum(catch_missing) == 0) {
    dat$catch_settings$fill_missing <- FALSE
  }
  if (sum(index_missing) == 0) {
    dat$index_settings$fill_missing <- FALSE
  }
  dat$fill_missing_map <- logical(length(dat$is_missing))
  dat$fill_missing_map[dat$obs_map$type == "catch"] <- dat$catch_settings$fill_missing
  dat$fill_missing_map[dat$obs_map$type == "index"] <- dat$index_settings$fill_missing
  dat$fill_missing_map <- dat$fill_missing_map & dat$is_missing
  dat$any_fill_missing <- any(dat$fill_missing_map)

  dat$sd_catch_modmat <- stats::model.matrix(dat$catch_settings$sd_form, data = dat$obs$catch)
  if (!is.null(dat$catch_settings$sd_supplied)) {
    if ("(Intercept)" %in% colnames(dat$sd_catch_modmat)) {
      dat$sd_catch_modmat <- stats::model.matrix(update(dat$catch_settings$sd_form, ~ 0 + .), data = dat$obs$catch)
      dat$catch_settings$sd_form <- update(dat$catch_settings$sd_form, ~ 0 + .)
      cli::cli_warn("Dropping intercept term in catch sd_form since supplied SDs are provided. Set sd_supplied to NULL to estimate the intercept.")
    }
    dat$log_sd_catch_supplied <- .log_supplied(dat$catch_settings$sd_supplied, dat$obs$catch, "catch sd_supplied")
  } else {
    dat$log_sd_catch_supplied <- rep(0, nrow(dat$obs$catch))
  }

  dat$sd_index_modmat <- stats::model.matrix(dat$index_settings$sd_form, data = dat$obs$index)
  if (!is.null(dat$index_settings$sd_supplied)) {
    if ("(Intercept)" %in% colnames(dat$sd_index_modmat)) {
      dat$sd_index_modmat <- stats::model.matrix(update(dat$index_settings$sd_form, ~ 0 + .), data = dat$obs$index)
      dat$index_settings$sd_form <- update(dat$index_settings$sd_form, ~ 0 + .)
      cli::cli_warn("Dropping intercept term in index sd_form since supplied SDs are provided. Set sd_supplied to NULL to estimate the intercept.")
    }
    dat$log_sd_index_supplied <- .log_supplied(dat$index_settings$sd_supplied, dat$obs$index, "index sd_supplied")
  } else {
    dat$log_sd_index_supplied <- rep(0, nrow(dat$obs$index))
  }

  q_design <- .parse_q_formula(dat$index_settings$q_form, dat$obs$index)
  dat[names(q_design)] <- q_design
  if (!is.null(dat$F_settings$mu_form)) {
    dat$F_modmat <- stats::model.matrix(F_settings$mu_form, data = dat$obs$catch)
  } else {
    dat$log_mu_f <- 0
    dat$F_modmat <- 0
  }

  if (!is.null(dat$M_settings$mu_form)) {
    dat$M_modmat <- stats::model.matrix(M_settings$mu_form, data = dat$obs$weight)
    if ("(Intercept)" %in% colnames(dat$M_modmat) && !is.null(dat$M_settings$mu_supplied)) {
      dat$M_modmat <- stats::model.matrix(update(M_settings$mu_form, ~ 0 + .), data = dat$obs$weight)
      cli::cli_warn("Dropping intercept term in M mu_form since supplied levels are provided. Set mu_supplied to NULL to estimate the intercept.")
    }
  } else {
    dat$mu_m <- 0
    dat$M_modmat <- 0
  }
  if (!is.null(dat$M_settings$mu_supplied)) {
    dat$log_mu_supplied_m <- .log_supplied(dat$M_settings$mu_supplied, dat$obs$weight, "M mu_supplied")
  } else {
    dat$log_mu_supplied_m <- 0
  }
  if (is.null(dat$M_settings$mu_form) && is.null(dat$M_settings$mu_supplied)) {
    cli::cli_abort("Please supply mu_supplied or mu_form for M.")
  }

  # model.matrix() can silently omit rows with missing covariates.
  design_tables <- c(sd_catch_modmat = "catch", sd_index_modmat = "index",
                     q_modmat = "index", F_modmat = "catch", M_modmat = "weight")
  for (nm in names(design_tables)) {
    x <- dat[[nm]]
    if (is.matrix(x) && (nrow(x) != nrow(dat$obs[[design_tables[[nm]]]]) ||
        any(!is.finite(x)))) {
      cli::cli_abort("{nm} must have one finite row per observation. Check for missing or non-finite formula covariates.")
    }
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

  .set_phi <- function(type) {
    if (type == "iid") {
      return(qlogis(c("age" = 0, "year" = 0)))
    }
    NULL
  }
  dat$logit_phi_n <- .set_phi(dat$N_settings$process)
  dat$logit_phi_f <- .set_phi(dat$F_settings$process)
  dat$logit_phi_m <- .set_phi(dat$M_settings$process)

  dat

}

# Evaluate positive offsets without silently dropping missing observation rows.
.log_supplied <- function(formula, data, label) {
  x <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  if (ncol(x) != 1L || !is.numeric(x[[1L]]) ||
      !nrow(x) %in% c(1L, nrow(data)) || any(!is.finite(x[[1L]])) ||
      any(x[[1L]] <= 0)) {
    cli::cli_abort("{label} must supply one positive finite value, or one per observation row.")
  }
  log(x[[1L]])
}
