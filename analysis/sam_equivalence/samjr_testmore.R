# SAMjr-style tests, adapted for tinyAM; not package code.
# Source: fishfollower/samjr, revision 0384774a3f2f8ec3dc3330375998bd73ef6bbe7c
#   testmore/{codIN3Ben,nscodXtraSd}/script.R: input and output workflows
#   samjr/R/tables.R, modeltable.default: comparison table below
#   samjr/R/fit.R, sam.fit: optimizer controls and three Newton steps below
#   Makefile, testmore-%: exact output comparison, ignoring trailing CR
# SAMjr authors: Anders Nielsen and Christoffer Moesgaard Albertsen.
# Upstream license: GPL-2. This adapted file: SPDX-License-Identifier: GPL-2.0-only
# See SOURCES.md for links, modifications and scope. No SAMjr likelihood is copied.

samjr_source_models <- function(case, folder, fit_function) {
  env <- new.env(parent = globalenv())
  env$sam.fit <- fit_function
  old <- setwd(folder); on.exit(setwd(old))
  last <- if (case == "nscodbiopro") "fit0" else if (case == "nscodXtraSd") "fit2" else "fit"
  for (expression in parse("script.R")) {
    if (identical(expression, quote(library(samjr)))) next
    eval(expression, env)
    if (exists(last, env, inherits = FALSE)) break
  }
  if (case == "nscodbiopro") list(baseline = env$fit0) else if (case == "nscodXtraSd") {
    list(Standard = env$fit, XtraSd = env$fit2)
  } else list(baseline = env$fit)
}

# Adapted from modeltable.default in samjr/R/tables.R. Only extraction changes:
# both native SAM and tinyAM expose opt$objective and opt$par. No fitted values
# are taken from the reference file. The LR p-value is an upstream test output,
# not a recommendation to use this test for an assessment model comparison.
samjr_modeltable <- function(fits) {
  nam <- names(fits)
  logL <- vapply(fits, function(f) -f$opt$objective, numeric(1))
  npar <- vapply(fits, function(f) length(f$opt$par), integer(1))
  res <- cbind("log(L)" = logL, "#par" = npar, "AIC" = -2 * logL + 2 * npar)
  rownames(res) <- nam
  o <- seq_along(fits)
  if (length(fits) == 2L) {
    o <- order(npar, decreasing = TRUE)
    if (npar[o[1L]] > npar[o[2L]]) {
      df <- npar[o[1L]] - npar[o[2L]]
      D <- 2 * (logL[o[1L]] - logL[o[2L]])
      P <- 1 - stats::pchisq(D, df)
      res <- cbind(res, c(NA, P)[o])
      colnames(res)[ncol(res)] <- paste0("Pval( ", nam[o[1L]], " -> ", nam[o[2L]], " )")
    }
  }
  res[o, , drop = FALSE]
}

# Adapted from samjr's sam.fit. Always applied to every protocol fit, never as
# a rescue step. The fresh objective uses tinyAM's likelihood, maps and random
# blocks and the source defpar() values, not fitted tinyAM/SAM parameter values.
samjr_optimizer <- function(fit, start) {
  dat <- fit$dat
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), start,
    random = fit$obj$env$.random, map = fit$parameter_map, silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
    control = list(eval.max = 2000, iter.max = 2000, rel.tol = 1e-10))
  steps <- 0L
  for (i in seq_len(3L)) {
    step <- try({
      g <- as.numeric(obj$gr(opt$par))
      h <- stats::optimHess(opt$par, obj$fn, obj$gr)
      solve(h, g)
    }, silent = TRUE)
    if (inherits(step, "try-error")) {
      warning("SAMjr Newton step failed: ", as.character(step), call. = FALSE)
      break
    }
    opt$par <- opt$par - step
    steps <- steps + 1L
  }
  opt$objective <- obj$fn(opt$par)
  sdrep <- RTMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  gradient <- as.numeric(obj$gr(opt$par))
  list(opt = opt, gradient = gradient, sdrep = sdrep, newton_steps = steps,
    initial_outer = obj$par, initial_parameters = start,
    converged = opt$convergence == 0L && isTRUE(sdrep$pdHess) &&
      all(is.finite(gradient)) && max(abs(gradient)) < fit$grad_tol)
}

samjr_write_output <- function(case, fits, path) {
  old <- options(digits = 7L, width = 80L, scipen = 0L)
  on.exit(options(old))
  if (case == "codIN3Ben") {
    cat(fits[[1L]]$opt$objective, "\n", file = path)
  } else if (case == "nscodXtraSd") {
    cat(capture.output(samjr_modeltable(fits)), file = path, sep = "\n")
  } else stop("No complete upstream output test is supported for ", case)
}

samjr_output_bytes <- function(path) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  # Matches diff --strip-trailing-cr; preserve spaces and the final newline.
  charToRaw(gsub("\r(?=\n|$)", "", rawToChar(bytes), perl = TRUE))
}

samjr_output_match <- function(actual, expected) {
  identical(samjr_output_bytes(actual), samjr_output_bytes(expected))
}

samjr_testmore <- function(case, folder, output) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
  })
  on.exit(writeLines(unique(warnings), file.path(output, "warnings.txt")), add = TRUE)
  # The upstream setup is executed without fitting: intercept sam.fit and keep
  # its data/conf/defpar arguments. This prevents optimized starts leaking in.
  specs <- capture(samjr_source_models(case, folder, function(data, conf, parameters, ...) {
    list(data = data, conf = conf, pl = parameters)
  }))
  fits <- protocol <- native <- list()
  for (name in names(specs)) {
    spec <- specs[[name]]
    input <- samjr_inputs(spec)
    dat <- do.call(prepare_tam, c(list(data = input$obs, years = input$years,
      ages = input$ages), input$settings))
    start <- samjr_parameters(spec, dat)
    saveRDS(list(source = spec, start = start, initialization = "source defpar()"),
      file.path(output, paste0(name, "-inputs.rds")))
    fits[[name]] <- capture(do.call(fit_tam, c(list(data = input$obs,
      years = input$years, ages = input$ages, start_par = start, silent = TRUE), input$settings)))
    saveRDS(list(source = spec, start = start, fit = fits[[name]]),
      file.path(output, paste0(name, "-public.rds")))
    protocol[[name]] <- capture(samjr_optimizer(fits[[name]], start))
    saveRDS(protocol[[name]], file.path(output, paste0(name, "-samjr-controls.rds")))
    native[[name]] <- capture(stockassessment::sam.fit(spec$data, spec$conf, spec$pl, silent = TRUE))
    saveRDS(native[[name]], file.path(output, paste0(name, "-native.rds")))
  }
  expected <- file.path(folder, "res.EXP")
  sources <- list(public = fits, samjr_controls = protocol, native = native)
  rows <- lapply(names(sources), function(engine) {
    path <- file.path(output, paste0("res.", engine, ".out"))
    samjr_write_output(case, sources[[engine]], path)
    numerical <- all(vapply(sources[[engine]], function(f) {
      if (engine == "public") isTRUE(f$is_converged) else if (engine == "samjr_controls") {
        isTRUE(f$converged)
      } else f$opt$convergence == 0L && isTRUE(f$sdrep$pdHess)
    }, logical(1)))
    data.frame(case = case, engine = engine, output_match = samjr_output_match(path, expected),
      numerical_success = numerical)
  })
  checks <- do.call(rbind, rows)
  write.csv(checks, file.path(output, "checks.csv"), row.names = FALSE)
  checks
}
