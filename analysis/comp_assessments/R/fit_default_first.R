# Used only by the new Pacific recipes. Existing recipes retain their starts.
fit_default_first <- function(fit_model, retry_start, attempt_file = NULL) {
  attempts <- list()
  save_attempts <- function() {
    if (!is.null(attempt_file)) {
      dir.create(dirname(attempt_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(do.call(rbind, attempts), attempt_file)
    }
  }
  attempt <- function(strategy, start = NULL) {
    cli::cli_inform(if (strategy == "default") {
      "Fitting with package starting values."
    } else {
      "Numerical checks failed; retrying once with aligned source starting values."
    })
    warnings <- character()
    begun <- Sys.time()
    number <- length(attempts) + 1L
    attempts[[number]] <<- data.frame(initialization = strategy,
      is_converged = NA, optimizer_code = NA_integer_, optimizer_message = "",
      objective = NA_real_, max_abs_gradient = NA_real_,
      sdreport_success = NA, positive_definite_hessian = NA,
      elapsed_seconds = NA_real_, notes = "Attempt started; no result returned yet.")
    save_attempts()
    fit <- tryCatch(withCallingHandlers(
      if (strategy == "default") fit_model() else fit_model(start_par = start),
      warning = function(w) warnings <<- c(warnings, conditionMessage(w))
    ), error = identity)
    checks <- if (inherits(fit, "tam_fit")) tinyAM::check_tam(fit) else NULL
    attempts[[number]] <<- data.frame(
      initialization = strategy,
      is_converged = !is.null(checks) && isTRUE(checks$is_converged),
      optimizer_code = if (is.list(fit$opt)) fit$opt$convergence else NA_integer_,
      optimizer_message = if (is.list(fit$opt) && length(fit$opt$message)) fit$opt$message else "",
      objective = if (is.list(fit$opt)) fit$opt$objective else NA_real_,
      max_abs_gradient = if (is.null(checks)) NA_real_ else checks$max_gradient,
      sdreport_success = is.list(fit$sdrep) && !inherits(fit$sdrep, "try-error"),
      positive_definite_hessian = isTRUE(fit$sdrep$pdHess),
      elapsed_seconds = as.numeric(difftime(Sys.time(), begun, units = "secs")),
      notes = paste(c(if (inherits(fit, "error")) conditionMessage(fit), warnings),
                    collapse = " | ")
    )
    save_attempts()
    fit
  }
  fit <- attempt("default")
  if (!attempts[[1L]]$is_converged) {
    starts <- tryCatch(retry_start(), error = identity)
    if (inherits(starts, "error")) {
      attempts[[1L]]$notes <- paste(attempts[[1L]]$notes,
        "Source starts unavailable:", conditionMessage(starts))
    } else {
      retry <- attempt("source_start", starts)
      if (inherits(retry, "tam_fit")) fit <- retry
    }
  }
  if (!inherits(fit, "tam_fit")) {
    cli::cli_abort(paste(vapply(attempts, function(x) x$notes, character(1)), collapse = " | "),
      initialization_attempts = do.call(rbind, attempts))
  }
  attr(fit, "initialization_attempts") <- do.call(rbind, attempts)
  fit
}
