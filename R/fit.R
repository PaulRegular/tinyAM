
#' Fit a Tiny Assessment Model (TAM)
#'
#' @title Fit TAM
#' @description
#' Builds data with [make_dat()], initializes parameters with [make_par()],
#' constructs the RTMB objective, optimizes it, and returns a fitted object with
#' reports and standard errors.
#'
#' @details
#' This function assigns the constructed data list to a global (`dat <<- ...`)
#' to satisfy closures inside the likelihood (`nll_fun`). Random-effect blocks
#' are chosen automatically from the model settings:
#' \itemize{
#'   \item Always includes \code{log_f} and \code{log_r}.
#'   \item Includes \code{missing} if there are missing observations.
#'   \item Includes \code{log_n} if \code{N_settings$process != "off"}.
#'   \item Includes \code{log_m} if \code{M_settings$process != "off"}.
#' }
#' A warning is issued if the number of random effects exceeds 1.5 times the
#' number of observed data points (rough identifiability check).
#'
#' @param inputs A named list of tidy observation tables as expected by
#'   [make_obs()] / [make_dat()] (e.g., \code{catch}, \code{index},
#'   \code{weight}, \code{maturity}).
#' @param silent Disable tracing information?
#' @inheritDotParams make_dat
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{call}{Matched call.}
#'   \item{dat}{Data list returned by [make_dat()].}
#'   \item{obj}{RTMB ADFun object.}
#'   \item{opt}{\code{nlminb} optimization result.}
#'   \item{rep}{List from \code{obj$report()}.}
#'   \item{sdrep}{\code{RTMB::sdreport} result.}
#' }
#'
#' @examples
#'
#' fit <- fit_tam(
#'   northern_cod_data,
#'   years = 1983:2024, ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "rw"),
#'   M_settings = list(process = "off", assumption = ~I(0.3)),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#' fit$sdrep
#'
#' @seealso [make_dat()], [make_par()], [sim_tam()], [fit_retro()],
#'   RTMB::MakeADFun(), RTMB::sdreport()
#' @export
fit_tam <- function(inputs, silent = FALSE, ...) {

  call <- match.call()

  dat <<- make_dat(inputs, ...) # dat not found without global assignment
  par <- make_par(dat)

  ran <- c("log_f", "log_r")
  if (length(par$missing) > 0) {
    ran <- c(ran, "missing")
  }
  if (dat$N_settings$process != "off") {
    ran <- c(ran, "log_n")
  }
  if (dat$M_settings$process != "off") {
    ran <- c(ran, "log_m")
  }

  n_ran <- length(unlist(par[ran]))
  n_obs <- sum(!is.na(dat$log_obs))
  if (n_ran > (n_obs * 1.5)) {
    warning(sprintf("Number of random effects (%d) exceed 1.5 times the number of observations (%d). Consider simplifying your model.", n_ran, n_obs))
  }

  obj <- MakeADFun(
    nll_fun,
    par,
    random = ran,
    silent = silent
  )

  opt <- nlminb(
    obj$par, obj$fn, obj$gr,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  opt$objective
  rep <- obj$report()
  sdrep <- sdreport(obj)

  out <- list(
    call = call,
    dat = dat,
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep
  )

  out

}


#' Run a retrospective (peel) analysis
#'
#' @title Retrospective fits for TAM
#' @description
#' Creates a sequence of terminal-year peels and refits the model on each
#' truncated year range, attaching the list of refits to \code{fit$retro}.
#'
#' @details
#' Peel years are \code{(max_year - folds) : max_year}. Each refit is attempted
#' with \code{try()} so individual failures do not stop the sequence. Refits are
#' generated via \code{update(fit, years = ...)}; ensure your \code{fit} object
#' supports \code{update()} with a \code{years} argument.
#'
#' @param fit A fitted TAM object as returned by [fit_tam()].
#' @param folds Integer number of terminal peels (default \code{2}).
#'
#' @return
#' The input \code{fit} with an added \code{$retro} list whose elements are the
#' refitted objects for each peel (named by terminal year).
#'
#' @examples
#' \dontrun{
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' retro_fit <- fit_retro(fit, folds = 5)
#' lapply(retro_fit$retro, function(x) x$rep$ssb)
#' }
#'
#' @seealso [fit_tam()], [sim_tam()]
#' @export
fit_retro <- function(fit, folds = 2) {

  min_year <- min(fit$dat$years)
  max_year <- max(fit$dat$years)
  retro_years <- (max_year - folds):(max_year)

  retro <- vector("list", length(retro_years))
  names(retro) <- retro_years
  for (i in seq_along(retro_years)) {
    retro[[i]] <- try(update(fit, years = min_year:retro_years[i]))
  }

  fit$retro <- retro
  fit

}


#' Simulate from a fitted TAM
#'
#' @title Simulation using a fitted TAM
#' @description
#' Simulates random effects and observations from a fitted model by calling the
#' likelihood with \code{simulate = TRUE}, then regenerates the REPORTed
#' quantities. Simulation is performed purely in R (not via RTMB \code{simref}).
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts MLEs as a parameter list: \code{as.list(fit$sdrep, "Estimate")}.
#'   \item Calls \code{nll_fun(p, simulate = TRUE)} to draw new random effects
#'         and \code{log_obs}.
#'   \item If \code{obs_only = FALSE} (default), injects the simulated random
#'         effects back into \code{p} (\code{p[obj$env$.random]}), calls
#'         \code{nll_fun} again (to update any derived quantities), and runs
#'         \code{obj$report(unlist(p))}.
#'   \item Returns the REPORT list with \code{log_obs} replaced by the simulated
#'         values; if \code{obs_only = FALSE}, all derived objects reflect the
#'         simulated random effects.
#' }
#'
#' @param fit A fitted TAM object returned by [fit_tam()].
#' @param obs_only Logical; if \code{TRUE}, only new \code{log_obs} values are
#'   simulated and inserted into the final report; if \code{FALSE} (default),
#'   simulated random effects are also injected and all REPORTed quantities are
#'   recomputed.
#'
#' @return
#' A list of REPORTed objects (same structure as \code{fit$rep}) with
#' \code{log_obs} replaced by simulated values. If \code{obs_only = FALSE},
#' the report also reflects simulated random effects.
#'
#' @examples
#' fit <- fit_tam(northern_cod_data, years = 1983:2024, ages = 2:14)
#' rep_sim <- sim_tam(fit, obs_only = FALSE)
#' str(rep_sim$log_obs)
#'
#' @export
sim_tam <- function(fit, obs_only = FALSE) {

  obj <- fit$obj
  p <- as.list(fit$sdrep, "Estimate")
  sims <- nll_fun(p, simulate = TRUE)
  rep <- obj$report()
  if (!obs_only) {
    p[obj$env$.random] <- sims[obj$env$.random]
    sims <- nll_fun(p, simulate = TRUE)
    rep <- obj$report(unlist(p))
  }
  rep$log_obs <- sims$log_obs
  rep

}

